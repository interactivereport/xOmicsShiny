###########################################################################################################
## Omics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: genesetmodule.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com); Xinmin Zhang (xinmin@bioinforx.com)
##@Date : 7/3/2023
##@version 3.0

##########################################################################################################
## Gene Set Enrichment
##########################################################################################################
# pkgs: "pathview", "biomaRt",  "ComplexHeatmap",  "DT", "stringr",  "dplyr", "tidyr"

# this module require results_long, tests, ProteinGeneName, Species, data_long(for heatmap)
library(tibble)
library(pathview)
library(biomaRt)
#library(enrichplot)
#library(enrichR)
library(ComplexHeatmap)


geneset_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
		column(3,
			wellPanel(
				uiOutput(ns('loadedprojects')),
				conditionalPanel("input.tabset!='KEGG Pathway View'",  ns=ns,
					radioButtons(ns("subset"),label="Use subset genes or upload your own subset?", choices=c("subset","upload genes"),inline = TRUE, selected="subset"),
					conditionalPanel(ns = ns, "input.subset=='subset'",
						selectInput(ns("test"), label="Select Comparison for Gene Set analysis", choices=NULL),
						column(width=6, numericInput(ns("FCcut"), label= "Fold Change Cutoff", value = 1.2, min=1, step=0.1)),
						column(width=6, numericInput(ns("pvalcut"), label= "P Value Cutoff", value=0.01, min=0, step=0.001)),
						radioButtons(ns("psel"), label= "P value or P.adj Value?", choices= c("Pval"="Pval", "Padj"="Padj"), inline = TRUE),
						radioButtons(ns("direction"), label= "Up- or Down-Regulated Genes?", choices= c("Both", "Up", "Down"), inline = TRUE),
						span(textOutput(ns("filteredgene1")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
						span(textOutput(ns("filteredgene2")), style = "color:red; font-size:15px; font-family:arial; font-style:italic")
					),
					conditionalPanel(ns = ns, "input.subset=='upload genes'",
						textAreaInput(ns("gene_list"), "Gene Name List", "", cols = 5, rows=6)
					),
				),
				conditionalPanel("input.tabset=='KEGG Pathway View'", ns= ns,
					radioButtons(ns("kegg_more_tests"), label= "Add more comparisons?", choices= c("Yes", "No"), selected="No", inline = TRUE),
					conditionalPanel("input.kegg_more_tests=='Yes'", ns= ns,
						selectInput(ns("test2"), label="2nd Comparison", choices=NULL),
						selectInput(ns("test3"), label="3rd Comparison", choices=NULL),
						selectInput(ns("test4"), label="4th Comparison", choices=NULL),
						selectInput(ns("test5"), label="5th Comparison", choices=NULL)
					)
				),
				conditionalPanel("input.tabset=='Gene Set Heatmap'", ns=ns,
					column(width=6, sliderInput(ns("hxfontsize_gsh"), "Column Font Size:", min = 2, max = 24, step = 1, value = 12)),
					column(width=6, sliderInput(ns("hyfontsize_gsh"), "Row Font Size:", min = 2, max = 24, step = 1, value = 12)),
					column(width=6, sliderInput(ns("htfontsize_gsh"), "Title Font Size:", min = 10, max = 32, step = 1, value = 14)),
					column(width=6, sliderInput(ns("hlfontsize_gsh"), "Legend Font Size:", min = 2, max = 20, step = 1, value = 10)),
					radioButtons(ns("gs_heatmap_label"), label="Gene Label",inline = TRUE, choices=c("UniqueID", "Gene.Name"), selected="Gene.Name")
				),
				conditionalPanel("input.tabset=='Gene Set Enrichment'",  ns=ns,
					radioButtons(ns("MSigDB"), label= "MSigDB Collections",
						choices= c("KEGG Pathway" = "KEGG",
							"c5.mf (GO molecular function)" = "c5.mf.v6.1",
							"c5.bp (GO biological process)" = "c5.bp.v6.1",
							"c5.cc (GO cellular component)" = "c5.cc.v6.1",
							"h.all.v6.1 (hallmark gene sets)" = "h.all.v6.1",
							"c2.cgp (chemical and genetic perturbations)" = "c2.cgp.v6.1",
							"c2.cp.biocarta" = "c2.cp.biocarta.v6.1",
							"c2.cp.reactome" =  "c2.cp.reactome.v6.1",
							"c2.cp (Canonical pathways)" = "c2.cp.v6.1",
							"c3.all (motif gene sets)" = "c3.all.v6.1",
							"c3.tft (transcription factor targets)" = "c3.tft.v6.1",
							"c6.all (oncogenic signatures)" = "c6.all.v6.1",
							"c7.all (immunologic signatures)" = "c7.all.v6.1",
							"Custom Geneset" = "Custom Geneset"
						),
					selected = "KEGG"),
					conditionalPanel(ns = ns, "input.MSigDB=='Custom Geneset'",
						textAreaInput(ns("CustomGeneset"), "Custom Geneset", "", cols = 5, rows=6)
					),
					numericInput(ns("orapvalcut"), label= "Gene Set Enrichment P Value Cutoff", value=0.05, min=0, step=0.001),
					radioButtons(ns("genebackground"),label="Gene Backgroud:", choices=c("All Genes" = "all genes", "Detected Genes"= "detected genes"),inline = TRUE, selected="all genes")
				),
				conditionalPanel("input.tabset=='KEGG Pathway View'", ns=ns,
					selectInput(ns("kegg_logFC"), label= "KEGG view log2FC Range:", choices= c(1, 2, 3), selected=1),
					radioButtons(ns("kegg_mapsample"), label= "Map Symbols to KEGG Nodes?", choices= c("YEs"=TRUE, "No"=FALSE), inline = TRUE),
				),
				conditionalPanel("input.tabset=='enrichR Plot'", ns=ns,
					selectInput(ns("SelectedPathway"), "Select Pathway", choices = NULL)
				)
			)
		),
		column(9,
			tabsetPanel(id=ns("tabset"),
				tabPanel(title="Gene Set Enrichment",value="Gene Set Enrichment",
					actionButton(ns("RunGSE"), "Run Gene Set Enrichment"),
					DT::dataTableOutput(ns("MSigDB"))
				),
				tabPanel(title="Gene Expression",value="Gene Expression",  textInput(ns('x1'), 'Row ID'),  DT::dataTableOutput(ns("Expression"))
				),
				tabPanel(title="Gene Set Heatmap",value="Gene Set Heatmap", textInput(ns('x2'), 'Row ID'),	actionButton(ns("genesetheatmap"), "Save to output"),	plotOutput(ns('SetHeatMap'), height=800)
				),
				tabPanel(title="KEGG Pathway View", value="KEGG Pathway View", textInput(ns('x3'), 'Row ID'), actionButton(ns("keggSave"), "Save to output"), plotOutput(ns("keggView"))
				),
				#tabPanel(title="enrichR Plot",value="enrichR Plot",
				#	actionButton(ns("RunEnrichR"), "Plot/Refresh (Gene set will be processed outside company!"),
				#	plotOutput(ns('RegPathway'), height=800)
				#),
				#tabPanel(title="enrichR Table",value="enrichR Table",
				#	DT::dataTableOutput(ns("RegPathwayTable"))
				#),
				tabPanel(title="Help", value="Help", htmlOutput("help_geneset")
				)
			)
		)
	)
}

geneset_server <- function(id) {
	shiny::moduleServer(id,
		function(input, output, session) {
			ns <- shiny::NS(id)

			output$loadedprojects <- renderUI({
				req(length(working_project()) > 0)
				radioButtons(ns("current_dataset"), label = "Change Working Dataset", choices=names(DataInSets), inline = F, selected=working_project())
			})

			observeEvent(input$current_dataset, {
				working_project(input$current_dataset)
			})

			#observe({
			#	is_enrichR_available <- require("enrichR")
			#	req(is_enrichR_available == TRUE)
			#	enrichRtab <- listEnrichrDbs()
			#	enrichRchoice <- enrichRtab[,"libraryName"]
			#	updateSelectizeInput(session,'SelectedPathway',choices=enrichRchoice, selected=enrichRchoice[145])
			#})
			#
			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$results_long)
				req(input$psel)

				results_long = DataInSets[[working_project()]]$results_long
				sel_test = input$test
				FCcut = log2(as.numeric(input$FCcut))
				pvalcut = as.numeric(input$pvalcut)

				if (input$psel == "Padj") {
					tmpdat = results_long %>% filter(test==sel_test & Adj.P.Value < pvalcut & abs(logFC) > FCcut)
				} else {
					tmpdat = results_long %>% filter(test==sel_test & P.Value < pvalcut & abs(logFC) > FCcut)
				}

				output$filteredgene1 <- renderText({paste("Genes Pass Cutoff (DEGs):",nrow(tmpdat),sep="")})
				output$filteredgene2 <- renderText({paste("Genes Up: ", sum(tmpdat$logFC>0), "; Genes Down: ", sum(tmpdat$logFC<0),sep="")})

				DEGs=tmpdat$Gene.Name
				updateTextAreaInput(session, "gene_list", value=paste(DEGs, collapse="\n"))
			})

			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$tests_order)
				tests = DataInSets[[working_project()]]$tests_order
				updateSelectizeInput(session,'test',choices=tests, selected=tests[1])
				tests_more=c("None", tests)
				updateSelectizeInput(session,'test2',choices=tests_more, selected="None")
				updateSelectizeInput(session,'test3',choices=tests_more, selected="None")
				updateSelectizeInput(session,'test4',choices=tests_more, selected="None")
				updateSelectizeInput(session,'test5',choices=tests_more, selected="None")
			})

			DataGenesetReactive <- reactive({
				req(DataInSets[[working_project()]]$results_long)
				results_long = DataInSets[[working_project()]]$results_long
				ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
				Prj_species=DataInSets[[working_project()]]$Species
				req(input$subset)

				if(!("hgnc" %in% names(LoadedData))){
					load("db/hgnc.RData")
					LoadedData[["hgnc"]]  <-  hgnc
				} else {
					hgnc <- LoadedData[["hgnc"]]
				}


				if (input$subset == "upload genes") {
					gene_list <- input$gene_list

					if(grepl("\n", gene_list)) {
						gene_list <-  stringr::str_split(gene_list, "\n")[[1]]
					} else if (grepl(",", gene_list)) {
						gene_list <-  stringr::str_split(gene_list, ",")[[1]]
					}

					gene_list <- gsub(" ", "", gene_list, fixed = TRUE)
					gene_list <- unique(gene_list[gene_list != ""])

					validate(need(length(gene_list)>2, message = "input gene list"))

					filteredgene <- results_long %>%
					dplyr::filter(!is.na(`Gene.Name`)) %>%
					dplyr::filter(Gene.Name %in% gene_list) %>%
					dplyr::select(one_of(c("Gene.Name","logFC"))) %>%
					dplyr::distinct(., Gene.Name,.keep_all = TRUE)
				}
				else {
					sel_test = input$test
					FCcut = log2(as.numeric(input$FCcut))
					pvalcut = as.numeric(input$pvalcut)

					filteredgene <- results_long %>%
					dplyr::mutate(P.stat=ifelse(rep(input$psel == "Padj", nrow(results_long)),  Adj.P.Value,  P.Value)) %>%
					dplyr::filter(test == sel_test) %>% dplyr::filter(!is.na(`Gene.Name`))%>%
					dplyr::filter(abs(logFC) >= FCcut & P.stat < pvalcut) %>%
					dplyr::select(one_of(c("Gene.Name","logFC"))) %>%
					dplyr::distinct(., Gene.Name,.keep_all = TRUE)
				}

				match_info = ifelse(Prj_species=="mouse", "(mouse genes mapped to human)", ifelse (Prj_species=="rat", "(rat genes mapped to human)", ""))

				if (!(Prj_species %in% c("human", "mouse", "rat")))  {
					match_info=str_c("(", Prj_species, " genes matches to human gene symbols)")
					Prj_species="human"
				}

				if (Prj_species=="human") {
					filteredgene <- filteredgene %>%
					dplyr::mutate_at(.vars = vars(Gene.Name), .funs = toupper)

					terminals.df <- dplyr::inner_join(hgnc, filteredgene, by=c("symbol"="Gene.Name"))

					all_genes <- dplyr::filter(ProteinGeneName, !is.na(`Gene.Name`)) %>%
					dplyr::mutate(Gene.Name=toupper(Gene.Name)) %>%
					dplyr::select(one_of(c("Gene.Name"))) %>%
					dplyr::inner_join(hgnc,., by=c("symbol"="Gene.Name")) %>%
					dplyr::pull(entrez_id) %>% unique()
					
				} else if  (Prj_species=="mouse") {
					load("db/mouse_rat_genes_map2_human.RData")
					terminals.df <- filteredgene %>%
					dplyr::transmute(mouse_Gene=Gene.Name, logFC) %>%
					dplyr::left_join(M_match%>%dplyr::transmute(mouse_Gene=mouse_symbol, symbol=human_symbol), by = join_by(mouse_Gene), relationship = "many-to-many") %>%
					dplyr::mutate(symbol=ifelse(is.na(symbol), toupper(mouse_Gene), symbol)) %>%
					dplyr::left_join(hgnc, by = join_by(symbol)) %>%
					dplyr::filter(!is.na(entrez_id))%>%
					dplyr::select(symbol, entrez_id, logFC, mouse_Gene)

					all_genes <- ProteinGeneName %>%
					dplyr::left_join(M_match%>%transmute(Gene.Name=mouse_symbol, symbol=human_symbol), by = join_by(Gene.Name))%>%
					dplyr::mutate(symbol=ifelse(is.na(symbol), toupper(Gene.Name), symbol)) %>%
					dplyr::left_join(hgnc, by = join_by(symbol))%>%
					dplyr::filter(!is.na(entrez_id))%>%
					dplyr::pull(entrez_id) %>% unique()
				} else if  (Prj_species=="rat") {
					load("db/mouse_rat_genes_map2_human.RData")
					terminals.df <- filteredgene %>%
					dplyr::transmute(rat_Gene=Gene.Name, logFC) %>%
					dplyr::left_join(R_match%>%transmute(rat_Gene=rat_symbol, symbol=human_symbol), by = join_by(rat_Gene), relationship = "many-to-many") %>%
					dplyr::mutate(symbol=ifelse(is.na(symbol), toupper(rat_Gene), symbol)) %>%
					dplyr::left_join(hgnc, by = join_by(symbol)) %>%
					dplyr::filter(!is.na(entrez_id))%>%
					dplyr::select(symbol, entrez_id, logFC, rat_Gene)

					all_genes <- ProteinGeneName %>%
					dplyr::left_join(R_match%>%dplyr::transmute(Gene.Name=rat_symbol, symbol=human_symbol), by = join_by(Gene.Name))%>%
					dplyr::mutate(symbol=ifelse(is.na(symbol), toupper(Gene.Name), symbol)) %>%
					dplyr::left_join(hgnc, by = join_by(symbol))%>%
					dplyr::filter(!is.na(entrez_id)) %>%
					dplyr::pull(entrez_id) %>% unique()
				}

				sig_genes <- as.numeric(as.data.frame(terminals.df)[,3])
				names(sig_genes) <- as.data.frame(terminals.df)[,2]
				detected_genes = all_genes

				if (input$genebackground == "all genes")
				all_genes <-  unique(hgnc$entrez_id)

				return(list("sig_genes" = sig_genes , "all_genes" = all_genes, "detected_genes" = detected_genes, "terminals.df" = terminals.df))
			})

			observeEvent(input$RunGSE, {withProgress(message = 'Processing...', value = 0, {

				getresults <- DataGenesetReactive()
				sig_genes <- 	getresults$sig_genes
				all_genes <- 	getresults$all_genes
				detected_genes <- getresults$detected_genes
				orapvalcut <- as.numeric(input$orapvalcut)

				if (input$MSigDB == "KEGG") {
					if(!("kegg.pathways" %in% names(LoadedData))){
						load("db/kegg.pathways.RData")
						LoadedData[["kegg.pathways"]]  <-  kegg.pathways
					} else{
						kegg.pathways <- LoadedData[["kegg.pathways"]]
					}
					gsets = kegg.pathways$human$kg.sets
				} else
				if  (input$MSigDB == "Custom Geneset") {

					CustomGeneset <- input$CustomGeneset

					if(grepl("\n", CustomGeneset)) {
						CustomGeneset <-  stringr::str_split(CustomGeneset, "\n")[[1]]
					} else if (grepl(",", CustomGeneset)) {
						CustomGeneset <-  stringr::str_split(CustomGeneset, ",")[[1]]
					}

					CustomGeneset <- gsub(" ", "", CustomGeneset, fixed = TRUE)
					CustomGeneset <- toupper(unique(CustomGeneset[CustomGeneset != ""]))

					validate(need(length(CustomGeneset)>2, message = "input custom geneset >2"))

					CustomGeneseList <- hgnc %>%
					dplyr::filter(!is.na(entrez_id)) %>%
					dplyr::filter(symbol %in% CustomGeneset) %>%
					dplyr::pull(entrez_id) %>% unique()  %>% as.character()

					gsets = list(CustomGeneset = CustomGeneseList)

				}	else {
					if(!("gmtlist" %in% names(LoadedData))){
						load("db/gmtlist.RData")
						LoadedData[["gmtlist"]]  <-  gmtlist
					} else {
						gmtlist <- LoadedData[["gmtlist"]]
					}

					gsets <- gmtlist[[input$MSigDB]]
				}

				gsa <- ORAEnrichment (deGenes=names(sig_genes), universe=all_genes, detected_genes=detected_genes, gsets, logFC =sig_genes, Dir=input$direction)
				res <- 	gsa %>%
				tibble::rownames_to_column(var="ID") %>%
				dplyr::filter(p.value < orapvalcut)

				res[,sapply(res,is.numeric)] <- signif(res[,sapply(res,is.numeric)],3)

				if (nrow(res) == 0) {
					res = data.frame("1"="", Notice = "No Gene Set passed P value cutoff, you can try increase P value cutoff")
				}

				#validate(need(nrow(res) > 0, message = "No Gene Set passed P value cutoff, you can try increase P value cutoff"))

				output$MSigDB <-  DT::renderDT(server=FALSE, {
					DT::datatable(res,  extensions = 'Buttons', selection = 'none', class = 'cell-border strip hover',
						options = list( dom = 'lBfrtip', pageLength = 15,
							buttons = list(
								list(extend = "csv", text = "Download Page", filename = "Page_results", exportOptions = list(modifier = list(page = "current"))
								),
								list(extend = "csv", text = "Download All", filename = "All_Results", exportOptions = list(modifier = list(page = "all"))
								)
							)
						)
					)  %>% formatStyle(1, cursor = 'pointer',color='blue')
					})
				})
			})

			observeEvent(input$MSigDB_cell_clicked, {
				info = input$MSigDB_cell_clicked
				if (is.null(info$value) || info$col != 1) return()
				updateTabsetPanel(session, 'tabset', selected = 'Gene Expression')
				updateTextInput(session, 'x1', value = info$value)
				updateTextInput(session, 'x2', value = info$value)
				updateTextInput(session, 'x3', value = info$value)
			})

			keggView_out <- reactive({withProgress(message = 'Making KEGG Pathway View...', value = 0, {
				req(DataInSets[[working_project()]]$tests)

				validate(need(input$MSigDB == "KEGG", message = "Only works on KEGG."))
				ID = input$x2
				validate(need(ID!="", message = "Select one geneset by clicking geneset name from 'Gene Set Enrichment' tab."))

				data_results = DataInSets[[working_project()]]$data_results
				Prj_species = DataInSets[[working_project()]]$Species

				getresults <- DataGenesetReactive()
				sig_genes <- 	getresults$sig_genes
				pid <- strsplit(ID," ")[[1]][1]

				tests=input$test
				if (input$kegg_more_tests=="Yes") {
					if (input$test2!="None") {tests=c(tests, input$test2)}
					if (input$test3!="None") {tests=c(tests, input$test3)}
					if (input$test4!="None") {tests=c(tests, input$test4)}
					if (input$test5!="None") {tests=c(tests, input$test5)}
				}

				selCol=rep(NA, length(tests))
				all_names=names(data_results)
				for (i in 1:length(tests)) {
					sel_i=which(str_detect(all_names, regex(str_c("^", tests[i]), ignore_case=T)) & str_detect(all_names, regex("logFC$", ignore_case=T)) )
					if (length(sel_i)==1) {
						selCol[i]=sel_i
					}
				}


				if (length(tests) > 1) {
					#			if (sum(is.na(selCol))==0) {
					img.file <- paste(pid,"pathview.multi.png",sep=".")
					if (file.exists(img.file)) {
						file.remove(img.file)
					}

					if (Prj_species=="human") {
						data_results<-data_results%>%left_join(hgnc, by=c("Gene.Name"="symbol"))
					} else if  (Prj_species=="mouse") {
						data_results<-data_results%>%mutate(mouse_Gene=Gene.Name)%>%left_join(M_match%>%transmute(mouse_Gene=mouse_symbol, symbol=human_symbol))%>%
						mutate(symbol=ifelse(is.na(symbol), toupper(mouse_Gene), symbol))%>%left_join(hgnc)
					} else if  (Prj_species=="rat") {
						data_results<-data_results%>%mutate(rat_Gene=Gene.Name)%>%left_join(R_match%>%transmute(rat_Gene=rat_symbol, symbol=human_symbol))%>%
						mutate(symbol=ifelse(is.na(symbol), toupper(rat_Gene), symbol))%>%left_join(hgnc)
					} else { #other species, use human
						data_results<-data_results%>%mutate(Gene.Name=toupper(Gene.Name))%>%left_join(hgnc, by=c("Gene.Name"="symbol"))
					}

					sel_gene=which(!is.na(data_results$entrez_id))
					FCdata=data.matrix(data_results[sel_gene, selCol])
					rownames(FCdata)=data_results$entrez_id[sel_gene]

					results_long = DataInSets[[working_project()]]$results_long

					FCdata <- results_long %>%
					dplyr::filter(test %in% tests) %>%
					dplyr::select(UniqueID, test, logFC) %>%
					tidyr::pivot_wider(names_from = test, values_from = logFC)


					tmp <- pathview(gene.data=FCdata, pathway.id=pid, kegg.dir="./kegg", kegg.native = T, species="hsa",low = "green", mid = "yellow", high = "red",
					same.layer = F, map.symbol=as.logical(input$kegg_mapsample), limit=list(gene=as.numeric(input$kegg_logFC), cpd=1) )
						#	if (ncol(FCdata)>1) {
						#		img.file <- paste(pid,"pathview.multi.png",sep=".")
						#	}

					} else {
						img.file <- paste(pid,"pathview","png",sep=".")
						if (file.exists(img.file)) {
							file.remove(img.file)
						}

						tmp <- pathview(gene.data=sig_genes, pathway.id=pid, kegg.dir="./kegg", kegg.native = T, species="hsa",low = "green", mid = "yellow", high = "red",
							same.layer = F, map.symbol=as.logical(input$kegg_mapsample), limit=list(gene=as.numeric(input$kegg_logFC), cpd=1)
						)
						#img.file <- paste(pid,"pathview","png",sep=".")
					}
					return(img.file)
				})
			})

			observeEvent(input$keggSave, {
				ID = input$x3
				img.file <- keggView_out()
				if (file.exists(img.file)) {
					img <- readPNG(img.file)
					saved_plots$keggSave[[ID]] <- img
				}
			})

			output$keggView = renderImage({
				img.file <- keggView_out()
				if (file.exists(img.file)) {
					list(src = img.file, contentType = 'image/png',	alt = "This is alternate text")
				}
			}, deleteFile = FALSE)

			genesetheatmap_out <- reactive({ withProgress(message = 'Making heatmap...', value = 0, {
				ID = input$x3
				validate(need(ID!="", message = "Select one geneset by clicking geneset name from 'Gene Set Enrichment' tab."))
				sel_group <- DataInSets[[working_project()]]$group_order
				getresults <- DataGenesetReactive()
				sig_genes <- 	getresults$sig_genes
				all_genes <- 	getresults$all_genes
				terminals.df <- getresults$terminals.df

				data_long = DataInSets[[working_project()]]$data_long
				ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
				Prj_species=DataInSets[[working_project()]]$Species

				if (input$MSigDB == "KEGG") {
					kegg.pathways <- LoadedData[["kegg.pathways"]]
					GenesetSig = kegg.pathways$human$kg.sets[[ID]]
				} else if  (input$MSigDB == "Custom Geneset") {
					CustomGeneset <- input$CustomGeneset

					if(grepl("\n", CustomGeneset)) {
						CustomGeneset <-  stringr::str_split(CustomGeneset, "\n")[[1]]
					} else if (grepl(",", CustomGeneset)) {
						CustomGeneset <-  stringr::str_split(CustomGeneset, ",")[[1]]
					}

					CustomGeneset <- gsub(" ", "", CustomGeneset, fixed = TRUE)
					CustomGeneset <- toupper(unique(CustomGeneset[CustomGeneset != ""]))

					validate(need(length(CustomGeneset)>2, message = "input custom geneset >2"))

					GenesetSig <- hgnc %>%
					dplyr::filter(!is.na(entrez_id)) %>%
					dplyr::filter(symbol %in% CustomGeneset) %>%
					dplyr::pull(entrez_id) %>% unique()

				}	else {
					GenesetSig <- gmtlist[[input$MSigDB]][[ID]]
				}

				terminalsdf.set <- dplyr::filter(terminals.df, entrez_id %in% GenesetSig)

				if (Prj_species=="human") {
					terminals_id <- dplyr::filter(ProteinGeneName, toupper(Gene.Name) %in% terminalsdf.set$symbol)  %>%
					dplyr::pull(UniqueID)
				} else if (Prj_species=="mouse"){
					terminals_id <- dplyr::filter(ProteinGeneName, Gene.Name %in% terminalsdf.set$mouse_Gene)  %>%
					dplyr::pull(UniqueID)
				}	else if (Prj_species=="rat"){
					terminals_id <- dplyr::filter(ProteinGeneName, Gene.Name %in% terminalsdf.set$rat_Gene)  %>%
					dplyr::pull(UniqueID)
				} else {  #no match, just use human
					terminals_id <- dplyr::filter(ProteinGeneName, toupper(Gene.Name) %in% terminalsdf.set$symbol)  %>%
					dplyr::pull(UniqueID)
				}

				subdatlong <- dplyr::filter(data_long, UniqueID %in% terminals_id ) %>%
				group_by(., group, UniqueID) %>%
				dplyr::summarise(mean=mean(expr, na.rm = TRUE), .groups = 'drop')
				subdatlong<-subdatlong%>%left_join(ProteinGeneName, by = "UniqueID")
				subdatwide <- subdatlong  %>%
				tidyr::spread(.,group, mean, fill = 0) %>%
				as.data.frame() %>%
				remove_rownames(.) %>%
				column_to_rownames(.,var="UniqueID") %>%
				dplyr::select(one_of(as.character(sel_group)))
				subdatwide=data.matrix(subdatwide)
				if (input$gs_heatmap_label=="Gene.Name") {
					sel_col=match(rownames(subdatwide), ProteinGeneName$UniqueID)
					rownames(subdatwide)=ProteinGeneName$Gene.Name[sel_col]
				}

				#remove rows with same values across all samples, which can cause hcluster error
				row_SD=apply(subdatwide, 1, function(x) sd(x,na.rm=T))
				subdatwide=subdatwide[row_SD!=0, ]

				scaled_data=t(scale(t(subdatwide))); scaled_data=pmin(scaled_data, 3); scaled_data=pmax(scaled_data, -3)
				p <- ComplexHeatmap::Heatmap(scaled_data,
					column_names_gp = gpar(fontsize = as.numeric(as.character(input$hxfontsize_gsh))),
					row_names_gp = gpar(fontsize = as.numeric(as.character(input$hyfontsize_gsh))),
					column_title_gp = gpar(fontsize = as.numeric(as.character(input$htfontsize_gsh))),
					heatmap_legend_param = list(
						title = "Z Score",
						color_bar = "continuous",
						title_gp = gpar(fontsize = as.numeric(as.character(input$hlfontsize_gsh))),
						labels_gp = gpar(fontsize = as.numeric(as.character(input$hlfontsize_gsh))-1)
					),column_title = ID, cluster_columns =F)
					return(p)
				})
			})

			output$SetHeatMap = renderPlot({
				ht <- genesetheatmap_out()
				ComplexHeatmap::draw(ht,  merge_legend=T,  auto_adjust = FALSE)
			})

			observeEvent(input$genesetheatmap, {
				#ID = input$x3
				# saved_plots$genesetheatmap[[ID]] <- genesetheatmap_out() #this only works on R4.0
				saved_plots$genesetheatmap<- genesetheatmap_out() #this works on R3.5 - 3.6
				#saved_plots$genesetheatmap <- genesetheatmap_out()$gtable
			})

			output$Expression <-  DT::renderDT(server=FALSE,{
				ID = input$x1
				validate(need(ID!="", message = "Select one geneset by clicking geneset name from 'Gene Set Enrichment' tab."))

				getresults <- DataGenesetReactive()
				terminals.df <- getresults$terminals.df

				if (input$MSigDB == "KEGG") {
					kegg.pathways <- LoadedData[["kegg.pathways"]]
					GenesetSig = kegg.pathways$human$kg.sets[[ID]]
				} else if  (input$MSigDB == "Custom Geneset") {
					CustomGeneset <- input$CustomGeneset

					if(grepl("\n", CustomGeneset)) {
						CustomGeneset <-  stringr::str_split(CustomGeneset, "\n")[[1]]
					} else if (grepl(",", CustomGeneset)) {
						CustomGeneset <-  stringr::str_split(CustomGeneset, ",")[[1]]
					}

					CustomGeneset <- gsub(" ", "", CustomGeneset, fixed = TRUE)
					CustomGeneset <- toupper(unique(CustomGeneset[CustomGeneset != ""]))
					validate(need(length(CustomGeneset)>2, message = "input custom geneset >2"))

					GenesetSig <- hgnc %>%
					dplyr::filter(!is.na(entrez_id)) %>%
					dplyr::filter(symbol %in% CustomGeneset) %>%
					dplyr::pull(entrez_id) %>% unique()

				}	else {
					GenesetSig <- gmtlist[[input$MSigDB]][[ID]]
				}

				terminalsdf.set <- dplyr::filter(terminals.df, entrez_id %in% GenesetSig)
				terminalsdf.set[,sapply(terminalsdf.set,is.numeric)] <- signif(terminalsdf.set[,sapply(terminalsdf.set,is.numeric)],3)

				DT::datatable(terminalsdf.set,  extensions = 'Buttons', options = list(
					dom = 'lBfrtip', pageLength = 15,
					buttons = list(
						list(extend = "csv", text = "Download Page", filename = "Page_results",
						exportOptions = list(modifier = list(page = "current"))),
							list(extend = "csv", text = "Download All", filename = "All_Results",
								exportOptions = list(modifier = list(page = "all"))
							)
						)
					)
				)

			})

			observeEvent(input$RunEnrichR, {
				getresults <- DataGenesetReactive()
				sig_genes <- 	getresults$sig_genes
				terminals.df <- getresults$terminals.df
				sig_genes_Dir=sig_genes
				if (input$direction=="Up") {
					sig_genes_Dir=sig_genes[terminals.df$logFC>0]
				} else if (input$direction=="Down") {
					sig_genes_Dir=sig_genes[terminals.df$logFC<0]
				}

				sig_symbols <- hgnc %>%
				dplyr::filter(entrez_id %in% names(sig_genes_Dir)) %>%
				dplyr::pull(symbol) %>% unique()

				dbs <- listEnrichrDbs()
				enrichRLive <- TRUE

				if (is.null(dbs)) {
					enrichRLive <- FALSE
				}

				dbs <- input$SelectedPathway
				enriched <- enrichr(sig_symbols , dbs)
				output$RegPathway <- renderPlot({
					plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value") + theme_bw(base_size = 20)
				})
				output$RegPathwayTable <- DT::renderDataTable({
					DT::datatable(enriched[[1]], options = list(lengthMenu = c(10,20,50,100,1000), pageLength = 10, scrollX = TRUE),
						selection=list(mode = "multiple")
					)
				})
			})
		}
	)
}
