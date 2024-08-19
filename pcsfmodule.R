###########################################################################################################
##This software belongs to Biogen Inc. All right reserved.
##
##@file: pcsfmodule.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 08/19/2024
##@version 3.0
###########################################################################################################
## PCSF Network
##########################################################################################################
#pkgs:  "PCSF", "visNetwork" ,"DT",   "dplyr"

library(PCSF)
library(visNetwork)
load("db/hgnc.RData")
load("db/kegg.pathways.RData")
data("STRING")

ORAEnrichment <- function(deGenes,universe, gsets, logFC, Dir="Both"){
	deGenes = deGenes[which(deGenes %in% universe)]
	tmp = rep(NA, length(gsets))
	ora.stats = data.frame(p.value=tmp, p.adj = tmp, DeGeneNum=tmp,DE_UpGene= tmp, DE_DownGene=tmp, SetNum = tmp, N_q=tmp, Fold_Enrich=tmp, SetNumAll=tmp, Total_DEG=tmp, Total_Gene=tmp,  DeGene_in_Set=tmp)

	totalDE = length(deGenes)
	if (Dir=="Up") {
		totalDE = sum(logFC > 0)
	} else if (Dir=="Down") {
		totalDE = sum(logFC < 0)
	}
	n = length(universe) - totalDE

	for (j in 1:length(gsets)){
		gset = gsets[[j]]
		DEinS = intersect(gset, deGenes)
		logFCinS = logFC[DEinS]
		totalDEinS = length(intersect(gset, deGenes))
		totalSinUniverse = length(intersect(gset, universe))

		N_q=totalDEinS- 0.5
		if (Dir=="Up") {
			N_q=length(logFCinS[logFCinS > 0])-0.5; DEinS<-names(logFCinS[logFCinS > 0])
		} else if (Dir=="Down") {
			N_q=length(logFCinS[logFCinS < 0])-0.5; DEinS<-names(logFCinS[logFCinS < 0])
		}

		ora.stats[j, "p.value"] = phyper(q = N_q, m=totalDE, n = n, k = totalSinUniverse, lower.tail = FALSE)
		ora.stats[j, "DeGeneNum"] = totalDEinS
		ora.stats[j, "SetNum"] = totalSinUniverse  #previous versions used length(gset)
		ora.stats[j, "DE_UpGene"] = length(logFCinS[logFCinS > 0])
		ora.stats[j, "DE_DownGene"] = length(logFCinS[logFCinS < 0])
		ora.stats[j, "N_q"]=N_q
		ora.stats[j, "SetNumAll"]=length(gset)
		ora.stats[j, "Fold_Enrich"]= round( (totalDEinS/totalDE) / (totalSinUniverse/length(universe) )*100)/100
		ora.stats[j, "Total_DEG"]=totalDE
		ora.stats[j, "Total_Gene"]=length(universe)
		ora.stats[j, "DeGene_in_Set"]=paste(DEinS, collapse = ",")

	}
	ora.stats[, "p.adj"] = p.adjust(ora.stats[, "p.value"], method = "BH")
	ora.stats[, "DE_Direction"]=Dir
	ora.stats<-ora.stats%>%mutate(GeneSet= names(gsets))%>%
	arrange(p.value, dplyr::desc(N_q))%>% rownames_to_column('rank')%>%dplyr::select(-N_q)%>%relocate(GeneSet)

	return(ora.stats)
}

pcsf_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
		column(3,
			wellPanel(
				uiOutput(ns('loadedprojects')),
				selectInput(ns("test"), label="Select Comparison Group", choices=NULL),

				fluidRow(
					column(width=6, numericInput(ns("FCcut"), label= "Fold Change Cutoff", value = 1.2, min=1, step=0.1)),
					column(width=6, numericInput(ns("pvalcut"), label= "P Value Cutoff", value=0.01, min=0, step=0.001))
				),
				radioButtons(ns("psel"), label= "P value or P.adj Value?", choices= c("Pval"="Pval", "Padj"="Padj"), inline = TRUE),
				span(textOutput(ns("filteredgene")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
				span(textOutput(ns("filteredgene2")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
				uiOutput(ns("myTabUI")),
				radioButtons(ns("PCSFMSigDB"), label= "MSigDB Collections",
					choices= c("KEGG Pathway" = "KEGG",
						"c2.cgp (chemical and genetic perturbations)" = "c2.cgp.v6.1",
						"c2.cp.biocarta" = "c2.cp.biocarta.v6.1",
						"c2.cp.reactome" =  "c2.cp.reactome.v6.1",
						"c2.cp (Canonical pathways)" = "c2.cp.v6.1",
						"c3.all (motif gene sets)" = "c3.all.v6.1",
						"c3.tft (transcription factor targets)" = "c3.tft.v6.1",
						"c5.bp (GO biological process)" = "c5.bp.v6.1",
						"c5.cc (GO cellular component)" = "c5.cc.v6.1",
						"c5.mf (GO molecular function)" = "c5.mf.v6.1",
						"c6.all (oncogenic signatures)" = "c6.all.v6.1",
						"c7.all (immunologic signatures)" = "c7.all.v6.1",
					"h.all.v6.1 (hallmark gene sets)" = "h.all.v6.1"),
				selected = "KEGG")
			)
		),
		column(9,
			tabsetPanel(id="tabset",
				tabPanel(title="PCSF network", value ="PCSF network",
					actionButton(ns("plotPCSF"), "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
				visNetworkOutput(ns("PCSFnetwork"), height="800px"), style = "background-color: #eeeeee;"),
					tabPanel(title="PCSF Function Cluster", value ="PCSF Function Cluster",
						actionButton(ns("plotPCSFCluster"), "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
					visNetworkOutput(ns("PCSFCluster"), height="800px"), style = "background-color: #eeeeee;"),
						tabPanel(title="Help", htmlOutput("help_PCSFnetwork"))
					)
				)
			)
		}

pcsf_server <- function(id) {
	shiny::moduleServer(id,
		function(input, output, session) {
			ns <- session$ns

			output$loadedprojects <- renderUI({
				req(length(working_project()) > 0)
				radioButtons(ns("current_dataset"), label = "Change Working Dataset", choices=DS_names(), inline = F, selected=working_project())
			})

			observeEvent(input$current_dataset, {
				working_project(input$current_dataset)
			})

			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$tests)
				tests = DataInSets[[working_project()]]$tests
				updateSelectizeInput(session,'test',choices=tests, selected=tests[1])
			})

			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$results_long)
				results_long <- DataInSets[[working_project()]]$results_long

				req(input$psel)
				p_sel   <- input$psel
				test_sel <- input$test
				FCcut <- log2(as.numeric(input$FCcut))
				pvalcut <- as.numeric(input$pvalcut)
				sel_label <- "UniqueID"
				direction <- "UpDown"
				tmpdat <- GeneFilter(results_long, test_sel, p_sel, direction, pvalcut, FCcut, sel_label)
				output$filteredgene <- renderText({paste("Genes Pass Cutoff (DEGs):",nrow(tmpdat),sep="")})
				output$filteredgene2 <- renderText({paste("Genes Up: ", sum(tmpdat$logFC>0), "; Genes Down: ", sum(tmpdat$logFC<0),sep="")})

				if (nrow(tmpdat) < 5) {
					output$myTabUI <- renderUI({
						tags$h4("Minimal 5 nodes. Try lower cutoffs.", style="color: red;")
					})
				} else if (nrow(tmpdat) > 500) {
					output$myTabUI <- renderUI({
						tags$h4("Maximal 500 nodes. Try higher cutoffs", style="color: red;")
					})
				} else {
					output$myTabUI <- renderUI({
						tags$h4("Generate plot by clicking plot/refersh button", style="color: green;")
					})
				}

			})

			DataTerminalsReactive <- reactive({

				results_long <- DataInSets[[working_project()]]$results_long
				results_long <- results_long %>% dplyr::mutate (Gene.Name = toupper(Gene.Name))
				ProteinGeneName <- DataInSets[[working_project()]]$ProteinGeneName
				ProteinGeneName <-ProteinGeneName %>%
				dplyr::mutate (Gene.Name = toupper(Gene.Name))

				req(input$psel)
				p_sel   <- input$psel
				test_sel <- input$test

				FCcut = log2(as.numeric(input$FCcut))
				pvalcut = as.numeric(input$pvalcut)
				sel_label <- "UniqueID"
				direction <- "UpDown"

				terminals.df <- GeneFilter(results_long, test_sel, p_sel, direction, pvalcut, FCcut, sel_label)
				req(nrow(terminals.df) > 0 & nrow(terminals.df) < 500)

				terminals.df <- terminals.df %>%
				dplyr::inner_join(hgnc, ., by=c("symbol"="Gene.Name"))

				all_genes <- dplyr::filter(ProteinGeneName, !is.na(`Gene.Name`)) %>%
				dplyr::select(one_of(c("Gene.Name"))) %>%
				dplyr::inner_join(hgnc,., by=c("symbol"="Gene.Name")) %>%
				dplyr::pull(entrez_id) %>% unique()

				terminals <- as.numeric(as.data.frame(terminals.df)[,"logFC"])
				terminals <- 2^abs(terminals)
				names(terminals) <- as.data.frame(terminals.df)[,"symbol"]

				return(list("terminals.df" = terminals.df, "terminals"= terminals, "all_genes"= all_genes))
			})

			DataPCSFReactive <- reactive({
				if (!exists("ppi")) {
					ppi <- construct_interactome(STRING)
				}

				terminaldata <- DataTerminalsReactive()
				terminals.df <- terminaldata$terminals.df
				terminals <- terminaldata$terminals
				all_genes <- terminaldata$all_genes
				set.seed(123)
				isolate({
					subnet <- PCSF::PCSF_rand(ppi, terminals, n=10, r = 0.1, w = 2, b = 1, mu = 0.0005)
				})

				return(list("subnet" = subnet, "terminals.df" = terminals.df, "terminals"= terminals, "all_genes"= all_genes))
			})

			observeEvent(input$plotPCSF,{
				output$PCSFnetwork <- renderVisNetwork({
					withProgress(message = 'Making Network:', value = 0, {
						isolate({
							PCSFData <- DataPCSFReactive()
							subnet <- PCSFData$subnet
							plot.PCSF(subnet, style = 1, edge_width=5, node_size=40, node_label_cex = 40, Steiner_node_color = "lightblue", Terminal_node_color = "lightgreen")
						})
					})
				})
			})

			observeEvent(input$plotPCSFCluster,{
				output$PCSFCluster <- renderVisNetwork({
					withProgress(message = 'Making Network:', value = 0, {
						isolate({
							PCSFData <-	DataPCSFReactive()
							subnet <- PCSFData$subnet
							terminals <- PCSFData$terminals
							terminals.df <- PCSFData$terminals.df
							all_genes <-    PCSFData$all_genes
							clusters = cluster_edge_betweenness(subnet)

							enrichment_result = as.list(1:length(clusters))
							enrichment_result_complete = as.list(1:length(clusters))
							for (a in 1:length(clusters)) {
								sig_genes.df <- terminals.df %>% dplyr::filter(symbol %in% clusters[[a]])
								sig_genes <- as.numeric(as.data.frame(sig_genes.df)[,"logFC"])
								sig_genes <- 2^abs(sig_genes)
								names(sig_genes) <- as.data.frame(sig_genes.df)[,"entrez_id"]

								if (input$PCSFMSigDB == "KEGG") {
									gsets = kegg.pathways$human$kg.sets
								} else {
									if (!exists("gmtlist")) {
										load("db/gmtlist.RData")
									}

									gsets <- gmtlist[[input$PCSFMSigDB]]
								}
								gsa <- ORAEnrichment(deGenes=names(sig_genes), universe = all_genes, gsets, logFC = sig_genes, Dir = "Both")

								res <- 	gsa %>%
								rownames_to_column(var="ID") %>%
								dplyr::filter(p.value < 0.05) %>%
								dplyr::slice(1:10) %>%
								dplyr::select(GeneSet, rank, p.value, p.adj, DeGeneNum, DE_UpGene, DE_DownGene, SetNum)

								res[,sapply(res,is.numeric)] <- signif(res[,sapply(res,is.numeric)],3)

								enrich = "<!DOCTYPE html> <html>\n
								<head>\n
								<style>\n
								table {font-family: arial, sans-serif; font-size: 10px; border-collapse: collapse; width: 100%;} \n
								td,\n
								th { border: 1px solid #dddddd; text-align: center; padding: 5px;}\n
								tr:nth-child(even) {background-color: #dddddd;}\n
								</style> \n
								</head>\n
								<body>\n
								<table><tr><th>ID</th><th>Rank</th><th>p.value</th><th>p.adj</th><th>DeGeneNum</th><th>UpGene</th><th>DownGene</th><th>SetNum</th></tr>\n"
								for (i in 1:nrow(res)) {
									enrich = paste0(enrich, "<tr>")
									for (j in 1:ncol(res)) {
										enrich = paste0(enrich, "<td>", res[i,
										j], "</td>")
									}
									enrich = paste0(enrich, "</tr>\n")
								}
								enrich = paste0(enrich, "</table> </body> </html>")
								enrichment_result[[a]] = enrich
								enrichment_result_complete[[a]] = res
							}

							enrichment = enrichment_result
							V(subnet)$group = clusters$membership
							V(subnet)$title = paste0("Cluster ", clusters$membership,": Enrichment analysis")
							for (i in 1:length(V(subnet))) {
								V(subnet)$title[i] = paste0(V(subnet)$title[i], enrichment[[V(subnet)$group[i]]])
							}
							class(subnet) <- c("PCSFe", "igraph")
							plot.PCSFe(subnet)
						})
					})
				})
			})
		}
	)
}
