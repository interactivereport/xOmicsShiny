###########################################################################################################
##This software belongs to Biogen Inc. All right reserved.
##
##@file: pcsfmodule.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 05/23/2023
##@version 3.0
###########################################################################################################
## PCSF Network
##########################################################################################################
#pkgs:  "PCSF", "visNetwork" ,"DT",   "dplyr"

library(PCSF)
library(visNetwork)
load("db/hgnc.RData")
load("db/kegg.pathways.RData")

pcsf_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
		column(3,
			wellPanel(
				selectInput(ns("test"), label="Select Test", choices=NULL),
				column(width=6, selectInput(ns("fccut"), label= "Choose Fold Change Threshold", choices= c("all"=0,"1.2"=0.263,"1.5"=0.584,"2"=1,"3"=1.585,"4"=2), selected = 0.263)),
				column(width=6, selectInput(ns("pvalcut"), label= "Choose P-value Threshold", choices= c("0.0001"=0.0001,"0.001"=0.001,"0.01"=0.01,"0.05"=0.05,"all"=1),selected=0.01)),
				radioButtons(ns("psel"), label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
				span(textOutput(ns("PCSFfilteredgene")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
				#textOutput(ns("PCSFfilteredgene")),
				#tags$head(tags$style("#PCSFfilteredgene{color: red;	font-size: 20px; font-style: italic;}")),
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
				selected = "KEGG"),
				uiOutput(ns("PCSFTabUI"))
			)
		),
		column(9,
			tabsetPanel(id="tabset",
				tabPanel("PCSF network", visNetworkOutput(ns("PCSFnetwork"), height="800px"), style = "background-color: #eeeeee;"),
				tabPanel("PCSF Function Cluster", visNetworkOutput(ns("PCSFCluster"), height="800px"), style = "background-color: #eeeeee;"),
				tabPanel(title="Data Output",	DT::dataTableOutput(ns("dat_PCSFnetwork"))),
				tabPanel(title="Help", htmlOutput("help_PCSFnetwork"))
			)#tabsetPanel
		)#column
	)
}

pcsf_server <- function(id) {
	shiny::moduleServer(id,
		function(input, output, session) {
			ns <- session$ns

			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$tests)
				tests = DataInSets[[working_project()]]$tests
				updateSelectizeInput(session,'test',choices=tests, selected=tests[1])
			})


			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$results_long)
				results_long = DataInSets[[working_project()]]$results_long
				test = input$test
				fccut = as.numeric(input$fccut)
				pvalcut = as.numeric(input$pvalcut)

				req(input$psel)
				if (input$psel == "Padj") {
					filteredgene = results_long %>%
					dplyr::filter(abs(logFC) > fccut & Adj.P.Value < pvalcut) %>%
					dplyr::filter(test == test)
				} else {
					filteredgene = results_long %>%
					dplyr::filter(abs(logFC) > fccut & P.Value < pvalcut) %>%
					dplyr::filter(test == test)
				}

				output$PCSFfilteredgene <- renderText({paste("Selected Genes:",nrow(filteredgene),sep="")})

				if (nrow(filteredgene) > 0 & 	nrow(filteredgene) < 500){
					output$PCSFTabUI <- renderUI({
						actionButton(ns("runPCSF"),"Generate PCSF")
					})

				} else {
					output$PCSFTabUI <- renderUI({
						"too many genes or zero gene"
					})
				}

			})



			DataPCSFReactive <- reactive({
				#browser()
				data("STRING")
				ppi <- construct_interactome(STRING)

				#sample_group <- DataInSets[[working_project()]]$sample_group
				results_long <- DataInSets[[working_project()]]$results_long
				ProteinGeneName <- DataInSets[[working_project()]]$ProteinGeneName
				test = input$test
				FCcut = as.numeric(input$fccut)
				pvalcut = as.numeric(input$pvalcut)
			
				if (input$psel == "Padj") {
					terminals.df = results_long %>%
					dplyr::filter(abs(logFC) > FCcut & Adj.P.Value < pvalcut)
				} else {
					terminals.df = results_long %>%
					dplyr::filter(abs(logFC) > FCcut & P.Value < pvalcut)
				}

				terminals.df <- terminals.df %>%
				dplyr::filter(test == test) %>%
				dplyr::filter(!is.na(`Gene.Name`)) %>%
				dplyr::arrange(desc(abs(logFC))) %>%
				distinct(.,Gene.Name,.keep_all = TRUE) %>%
				as.data.frame()

				terminals.df <- dplyr::inner_join(hgnc,terminals.df, by=c("symbol"="Gene.Name"))


				all_genes <- dplyr::filter(ProteinGeneName, !is.na(`Gene.Name`)) %>%
				dplyr::select(one_of(c("Gene.Name"))) %>%
				dplyr::inner_join(hgnc,., by=c("symbol"="Gene.Name")) %>%
				#dplyr::select(one_of(c("entrez_id"))) %>% collect %>%	.[["entrez_id"]] %>% as.character() %>% unique()
				dplyr::pull(entrez_id) %>% unique()

				terminals <- as.numeric(as.data.frame(terminals.df)[,"logFC"])
				terminals <- 2^abs(terminals)
				names(terminals) <- as.data.frame(terminals.df)[,"symbol"]

				set.seed(123)
				subnet <- PCSF::PCSF_rand(ppi, terminals, n=10, r = 0.1, w = 2, b = 1, mu = 0.0005)
				return(list("subnet" = subnet, "terminals.df" = terminals.df, "terminals"= terminals, "all_genes"= all_genes))
			})


			#library(RMySQL)
			#mydb = dbConnect(MySQL(), user='bgao', password='sqladmin', host='camhpcprot01.biogen.com', dbname="msdata")
			#
			#statement = "SELECT PARTICIPANT_A, PARTICIPANT_B, INTERACTION_TYPE, INTERACTION_PUBMED_ID,INTERACTION_DATA_SOURCE FROM `pathwaycommonsv9` WHERE `INTERACTION_TYPE` IN ('controls-state-change-of', 'controls-phosphorylation-of', 'interacts-with','in-complex-with')"
			#rs <- dbSendQuery(mydb, statement)
			#sql.r = fetch(rs, n = -1)
			#dbClearResult(rs)
			#
			#sql.cost <-
			#dplyr::mutate(sql.r, cost = if_else (INTERACTION_TYPE %in% c('controls-state-change-of', 'controls-phosphorylation-of'), 0.6, 0.1)) %>%
			#dplyr::mutate(PubmedN = str_count(INTERACTION_PUBMED_ID,';')) %>%
			#dplyr::mutate(
			#	cost2 = case_when(
			#		PubmedN <= 1 ~ 0,
			#		PubmedN == 2 ~ 0.1,
			#		PubmedN == 3 ~ 0.2,
			#		PubmedN == 4 ~ 0.3,
			#		PubmedN > 4 ~ 0.4,
			#	)
			#) %>%
			#dplyr::mutate(cost = cost + cost2) %>%
			#dplyr::select(one_of("PARTICIPANT_A", "PARTICIPANT_B", "cost")) %>%
			#dplyr::arrange(., desc(cost)) %>%
			#dplyr::rename(from=PARTICIPANT_A, to=PARTICIPANT_B) %>%
			#dplyr::distinct(from,to,.keep_all=TRUE)
			#ppi <- construct_interactome(sql.cost)

			observeEvent(input$runPCSF,{
				output$PCSFnetwork <- renderVisNetwork({
					withProgress(message = 'Making Network:', value = 0, {
						isolate({
							PCSFData <-	DataPCSFReactive()
							subnet <- PCSFData$subnet
							plot.PCSF(subnet, style = 1, edge_width=5, node_size=40, node_label_cex = 40, Steiner_node_color = "lightblue", Terminal_node_color = "lightgreen")
						})
					})
				})
			})


			observeEvent(input$runPCSF,{
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
									gsets <- gmtlist[[input$PCSFMSigDB]]
								}
								#gsa <- ORAEnrichment (deGenes=names(sig_genes), universe=all_genes, gsets, logFC =sig_genes)

								gsa <- ORAEnrichment(deGenes=names(sig_genes), universe = all_genes, detected_genes = all_genes, gsets, logFC = sig_genes, Dir = "Both")

								res <- 	gsa %>%
								rownames_to_column(var="ID") %>%
								dplyr::filter(p.value < 0.05) %>%
								dplyr::slice(1:10)

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
			#termlist <- c()
			#for ( i in 1:length(enrichment_result_complete)) {
			#	cluster <- 	enrichment_result_complete[[i]]
			#	colnames(cluster) <- make.names(colnames(cluster))
			#	cluster <- dplyr::filter(cluster, p.value < 0.01) %>%
			#	#dplyr::select(one_of("Term", "Overlap","p.value", "Genes","cluster")) %>%
			#	dplyr::mutate(Adjusted.P.value =  round(as.numeric(p.value),6)) %>%
			#	dplyr::arrange(Adjusted.P.value) %>%
			#	dplyr::slice(1:10) %>%
			#	dplyr::mutate(cluster = i) %>%
			#	as.data.frame()
			#	termlist <- rbind(termlist,cluster)
			#}
			#termlist
			#
		}
	)
}
