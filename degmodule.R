###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: degmodule.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 05/23/2023
##@version 3.0
###########################################################################################################
#pkgs:"DT", "shinyjqui","dplyr","stringr", "plotly", "rlang","tidyr","coop","ggrastr"
#data req:
library("coop")

deg_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
		column(3,
			wellPanel(
				conditionalPanel(ns = ns, "input.tabset =='DEG Counts' || input.tabset =='Volcano Plot' || input.tabset == 'Data Table'",
					uiOutput(ns('loadedprojects')),
					selectInput(ns("test"), label="Select Comparison Groups for Volcano Plot", choices=NULL)
				),
				conditionalPanel(ns = ns, "input.tabset == 'Multi Volcano Plots' || input.tabset == 'DEGs in Two Comparisons'",
					uiOutput(ns('loaddatasets')),
					conditionalPanel(ns = ns, "input.tabset=='Multi Volcano Plots'",
						radioButtons(ns("SeparateOnePlot"), label="Separate or One Plot", inline = TRUE, choices = c("Separate" = "Separate", "OnePlot" = "OnePlot")),
						conditionalPanel(ns = ns, "input.SeparateOnePlot=='OnePlot'",
							radioButtons(ns("OnePlotColorBy"), label="Color By", inline = TRUE, choices = c("Project and Test" = "project_test", "Test" = "test", "Project" = "ProjectID")),
						)
					),
					conditionalPanel(ns = ns, "input.tabset=='DEGs in Two Comparisons'",
						radioButtons(ns("merged_by"),label="Select merged by", inline = TRUE, choices="")
					)
				),
				conditionalPanel(ns = ns, "input.tabset != 'Help'",
					fluidRow(
						column(width=6, numericInput(ns("FCcut"), label= "Fold Change Cutoff", value = 1.2, min=1, step=0.1)),
						column(width=6, numericInput(ns("pvalcut"), label= "P Value Cutoff", value=0.01, min=0, step=0.001))
					),
					radioButtons(ns("psel"), label= "P value or P.adj Value?", choices= c("Pval"="Pval", "Padj"="Padj"), inline = TRUE)
				),
				conditionalPanel(ns = ns, "input.tabset =='Volcano Plot'",
					span(textOutput(ns("filteredgene")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
					span(textOutput(ns("filteredgene2")), style = "color:red; font-size:15px; font-family:arial; font-style:italic")
				),
				conditionalPanel(ns = ns, "input.tabset =='Volcano Plot' || input.tabset == 'Multi Volcano Plots' || input.tabset == 'DEGs in Two Comparisons'",
					radioButtons(ns("interactive"), label="Static or Interactive Plot", inline = TRUE, choices = c("Static" = "static", "Interactive" = "interactive")),
					radioButtons(ns("genelabel"),label="Select Gene Label", inline = TRUE, choices=""),
					radioButtons(ns("label"), label="Label Genes:", inline = TRUE, choices = c("DEGs", "Upload", "Geneset"), selected = "DEGs"),
					sliderInput(ns("Ngenes"), "# of Genes to Label", min = 0, max = 200, step = 5, value = 50, width = "80%"),
					conditionalPanel(ns = ns, "input.label=='Upload'", textAreaInput(ns("gene_list"), "List of genes to label", "", cols = 5, rows=6)),
					conditionalPanel(ns = ns, "input.label=='Geneset'",
						selectizeInput(ns("sel_geneset"), label="Available GeneSet", choices = NULL, multiple = FALSE),
						textAreaInput(ns("geneset_genes"), "Genes in Geneset", "", cols = 5, rows=6)
					),
					fluidRow(
						column(width=5, numericInput(ns("Max_logFC"), label= "Max log2FC", value=0, min=0)),
						column(width=7, numericInput(ns("Max_Pvalue"), label= "Max -log10(Pvalue)", value=0, min=0))
					),
					fluidRow(
						column(width=6, sliderInput(ns("lfontsize"), "Label Font Size:", min = 1, max = 10, step = 1, value = 4, width = "80%")),
						column(width=6, sliderInput(ns("yfontsize"), "Legend Font Size:", min = 8, max = 24, step = 1, value = 14,  width = "80%"))
					),
					fluidRow(
						column(width=6, radioButtons(ns("vlegendpos"), label="Legend position", inline = TRUE, choices = c("bottom","right"), selected = "bottom")),
						column(width=6, radioButtons(ns("rasterize"), label="Rasterize plot", inline = TRUE, choices = c("Yes","No"), selected = "No"))
					),
					conditionalPanel(ns = ns, "input.tabset =='Volcano Plot'",
						fluidRow(
							column(width=6, radioButtons(ns("showpep"), label="Show Unique Peptide (if unique.peptides column is in result table) ", inline = TRUE, choices = c("Yes","No"), selected = "No"))
						)
					),
					conditionalPanel(ns = ns, "input.tabset=='DEGs in Two Comparisons'",
						fluidRow(
							column(width=6, radioButtons(ns("DEG_comp_XY"), label="Same X and Y scale?", inline = TRUE, choices = c("Yes","No"), selected = "Yes")),
							column(width=6, radioButtons(ns("DEG_comp_color"), label="Same label and dot color?", inline = TRUE, choices = c("Yes","No"), selected = "Yes"))
						),
						radioButtons(ns("LM"), label="Add linear regression?", inline = TRUE, choices = c("Yes","No"), selected = "Yes"),
						fluidRow(
							column(width=6, textInput(ns("Xlab"), "X label", width = "100%")),
							column(width=6, textInput(ns("Ylab"), "Y label", width = "100%"))
						),
						sliderInput(ns("FC_Cutoff"), "Remove Data by FC Cutoff:", min = 1, max = 2, step = 0.1, value = 1)
					)
				)
			)
		),
		column(9,
			tabsetPanel(id=ns("tabset"),
				tabPanel(title="DEG Counts", value ="DEG Counts", tags$p("Click a comparison name to view volcano plot."), tags$hr(style="border-color: black;") , DT::dataTableOutput(ns("deg_counts"))
				),
				tabPanel(title="Volcano Plot", value ="Volcano Plot",
					tabsetPanel(id=ns("Volcano_tabset"),
						tabPanel(title="Volcano Plot", value ="Volcano Plot", actionButton(ns("plot_volcano"), "Plot/Refresh"), actionButton(ns("volcano"), "Save to output"),
							uiOutput(ns("plot.volcano"))
						),
						tabPanel(title="Data Table", value = "Data Table",
							DT::dataTableOutput(ns("volcanoData")),
							actionButton(ns("DEG_data"), "Save data for output as one Excel file")
						)
					)
				),
				tabPanel(title="DEGs in Two Comparisons", value = "DEGs in Two Comparisons",
					tabsetPanel(id=ns("Comparisons_tabset"),
						tabPanel(title="DEGs in Two Comparisons", value = "DEGs in Two Comparisons", actionButton(ns("plot_DEG_comp"), "Plot/Refresh"),  actionButton(ns("DEG_comp"), "Save to output"),
							fluidRow(
								uiOutput(ns("plot.DEG_Compare"))
							),
							fluidRow(
								hr(),
								DT::dataTableOutput(ns("LMResult"))
							)
						),
						tabPanel(title="Comparisons Data Table", value = "Comparisons Data Table",
							DT::dataTableOutput(ns("DEGCompareData")),
							actionButton(ns("DEGCompareData"), "Save data for output as one Excel file"),
							shiny::downloadButton(outputId = ns("download_data_button"),  label = "Download All as CSV file")
						)
					)
				),
				tabPanel(title="Multi Volcano Plots", value = "Multi Volcano Plots", actionButton(ns("plot_Multivolcano"), "Plot/Refresh"), actionButton(ns("Multivolcano"), "Save to output"),
					uiOutput(ns("plot.Multivolcano"))
				),
				tabPanel(title="Help", value = "Help", htmlOutput("help_volcano")
				)
			)
		)
	)
}

deg_server <- function(id) {
	shiny::moduleServer(id,
		function(input, output, session) {
			ns <- shiny::NS(id)

			observeEvent(input$current_dataset, {
				working_project(input$current_dataset)
			})

			output$loadedprojects <- renderUI({
				req(length(working_project()) > 0)
				radioButtons(ns("current_dataset"), label = "Change Working Dataset", choices=DS_names(), inline = F, selected=working_project())
			})

			output$loaddatasets <- renderUI({
				req(length(working_project()) > 0)
				projectlist <- list()
				for (project in DS_names()) {
					projectlist <-	append(projectlist, paste(project, DataInSets[[project]]$tests_order, sep="->"))
				}
				tagList(
					shinyjqui::sortableCheckboxGroupInput(ns("add_dataset"), label = "Add Dataset to Compare", choices = projectlist)
				)
			})

			maxcomparison <- reactiveVal(NULL)
			observeEvent(input$tabset, {
				if (input$tabset == "DEGs in Two Comparisons") {
					maxcomparison(2)
				} else {
					maxcomparison(8)
				}
			})

			observe({
				req(length(working_project()) > 0)
				req(input$tabset)
				req(input$add_dataset)
				projectlist <- list()
				for (project in DS_names()) {
					projectlist <-	append(projectlist, paste(project, DataInSets[[project]]$tests_order, sep="->"))
				}
				if(length(input$add_dataset) > maxcomparison()){
					updateCheckboxGroupInput(session, "add_dataset", selected= head(input$add_dataset, maxcomparison()))
				}
			})

			observeEvent(input$label , {
				req(length(working_project()) > 0)
				req(input$label == "Geneset")
				genesetnames <- GetGeneSetNames()
				updateSelectizeInput(session, "sel_geneset", choices =  c('Type to Search' = '', genesetnames), server = TRUE)
			})

			observeEvent(input$sel_geneset, {
				req(length(working_project()) > 0)
				req(input$label == "Geneset")
				req(input$sel_geneset!="")
				sel_geneset <- input$sel_geneset
				geneset_genenames <- GetGenesFromGeneSet(sel_geneset)
				updateTextAreaInput(session, "geneset_genes", value=paste(geneset_genenames, collapse=","))
			})

			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$tests_order)
				tests = DataInSets[[working_project()]]$tests_order
				ProteinGeneNameHeader = DataInSets[[working_project()]]$ProteinGeneNameHeader
				updateRadioButtons(session,'genelabel', inline = TRUE, choices=ProteinGeneNameHeader, selected="Gene.Name")
				updateRadioButtons(session,'merged_by', inline = TRUE, choices=ProteinGeneNameHeader, selected="UniqueID")
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

				DEGs=tmpdat$Gene.Name

				if (nrow(tmpdat) > input$Ngenes) {
					DEGs=sample(DEGs, input$Ngenes)
				}
				updateTextAreaInput(session, "gene_list", value=paste(DEGs, collapse="\n"))
			})

			DatavolcanoReactive <- reactive({
				req(length(working_project()) > 0)
				results_long = DataInSets[[working_project()]]$results_long
				data_results = DataInSets[[working_project()]]$data_results
				sel_test = input$test
				FCcut = log2(as.numeric(input$FCcut))
				FCcut_rd = round(FCcut*1000)/1000
				pvalcut = as.numeric(input$pvalcut)
				genelabel = input$genelabel
				Max_Pvalue <- input$Max_Pvalue
				Max_logFC <- input$Max_logFC

				psel <- input$psel
				if (psel == "Pval")
				psel = "P.Value"
				if (psel == "Padj")
				psel = "Adj.P.Value"

				res = results_long %>%
				dplyr::filter(test==sel_test) %>%
				dplyr::mutate_if(is.factor, as.character) %>%
				dplyr::mutate(labelid = !!sym(genelabel)) %>%
				dplyr::filter(!is.na(P.Value)) %>%
				dplyr::mutate(
					Significance = case_when(
						abs(logFC) >= FCcut & !!sym(psel) <= pvalcut ~ paste0(psel,"<", pvalcut," & abs(log2FC)>=", FCcut_rd),
						abs(logFC) <  FCcut & !!sym(psel) <= pvalcut ~ paste0(psel,"<", pvalcut," & abs(log2FC)<", FCcut_rd),
					TRUE ~ "Not Significant")
				) %>%
				dplyr::mutate(Significance = factor(Significance,	levels = c("Not Significant",	paste0(psel,"<", pvalcut," & abs(log2FC)>=", FCcut_rd),	paste0(psel,"<", pvalcut," & abs(log2FC)<", FCcut_rd)))
				) %>%
				dplyr::mutate(logFC_ori = logFC)


				if (Max_Pvalue > 0) {
					if (psel == "Padj")
					res <- res %>% dplyr::mutate(Adj.P.Value=pmax(Adj.P.Value, 10^(0-Max_Pvalue) ))
					if (psel == "P.Value")
					res <- res %>% dplyr::mutate(P.Value=pmax(P.Value, 10^(0-Max_Pvalue) ))
				}

				if (Max_logFC > 0) {
					res <- res %>% dplyr::mutate(logFC=ifelse(logFC>=0, pmin(Max_logFC, logFC), pmax(0-Max_logFC, logFC)))
				}

				if (input$showpep == "Yes" & ("unique.peptides" %in% colnames(data_results))) {
					res <- res %>%
					dplyr::left_join(.,(data_results %>% dplyr::select(id, unique.peptides)), by="id")
				}
				return(res)
			})

			volcanoplot_out <- eventReactive(input$plot_volcano, {
				req(length(working_project()) > 0)
				res = DatavolcanoReactive()
				ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
				sel_test = input$test
				FCcut = log2(as.numeric(input$FCcut))
				FCcut_rd=round(FCcut*1000)/1000
				pvalcut = as.numeric(input$pvalcut)

				psel <- input$psel
				if (psel == "Pval")
				psel = "P.Value"
				if (psel == "Padj")
				psel = "Adj.P.Value"

				if (psel == "Adj.P.Value") {
					if (input$showpep == "Yes" & ("unique.peptides" %in% colnames(res))) {  #02122024 by bgao
						res <- res %>%
						dplyr::select(UniqueID, labelid, logFC, Adj.P.Value, Significance,unique.peptides,logFC_ori)
					} else {
						res <- res %>%
						dplyr::select(UniqueID, labelid, logFC, Adj.P.Value, Significance, logFC_ori)
					}

					p <- ggplot(res, aes(x = logFC, y = -log10(Adj.P.Value), text = labelid))
					ylab <- "-log10(Padj.Value)"

					filterSig <- paste0(psel,"<", pvalcut," & abs(log2FC)>=", FCcut_rd)
					data.label <- filter(res, Significance == filterSig)

					if (nrow(data.label) > input$Ngenes) {
						data.label <- top_n(data.label, input$Ngenes, abs(logFC_ori))
					}

				} else {
					if (input$showpep == "Yes" & ("unique.peptides" %in% colnames(res))) { #02122024 by bgao
						res <- res %>%
						dplyr::select(UniqueID, labelid, logFC, P.Value, Significance, unique.peptides, logFC_ori)
					} else {
						res <- res %>%
						dplyr::select(UniqueID, labelid, logFC, P.Value, Significance, logFC_ori)
					}

					filterSig <- paste0(psel,"<", pvalcut," & abs(log2FC)>=", FCcut_rd)
					data.label <- filter(res, Significance == filterSig)

					if (nrow(data.label) > input$Ngenes) {
						data.label <- top_n(data.label, input$Ngenes, abs(logFC_ori))
					}
					p <- ggplot(res, aes(x = logFC, y = -log10(P.Value), text = labelid))
					ylab <- "-log10(P.Value)"
				}

				if (input$label=="Upload" | input$label=="Geneset") {
					if (input$label=="Upload") {
						req(input$gene_list)
						gene_list <- input$gene_list
					} else {
						req(input$geneset_genes)
						gene_list <- input$geneset_genes
					}

					gene_list <- ProcessUploadGeneList(gene_list)
					validate(need(length(gene_list)>0, message = "input gene list"))
					uploadlist <- dplyr::filter(ProteinGeneName, (UniqueID %in% gene_list) | (Protein.ID %in% gene_list) | (toupper(Gene.Name) %in% toupper(gene_list )))  %>%
					dplyr::pull(UniqueID)

					data.label <- res %>%
					dplyr::filter(UniqueID %in% uploadlist)
					validate(need(nrow(data.label)>0, message = "no gene found in result"))
				}

				p <- p +
				scale_color_manual(values = c("grey", "red2","green2"))

				if (input$rasterize=="Yes") {
					p <- p + geom_point_rast(aes(color = Significance ), size=0.7, alpha=0.6, na.rm=TRUE, dev="ragg")
				} else {
					p <- p + geom_point(aes(color = Significance ), size=0.7)
				}

				if (input$showpep == "Yes" & ("unique.peptides" %in% colnames(res))) {  #02122024 by bgao  rasterize doesn't work here
					res$labelid <- paste(
						"Name:", res$labelid,
						"<br>",
						"Unique Peptides:", res$unique.peptides
					)
					res <- res %>%
					dplyr::mutate(unique.peptides = ifelse(unique.peptides >= 3, 3, unique.peptides))

					p <- ggplot(res, aes(x = logFC, y = -log10(P.Value), text = labelid)) +
					scale_color_manual(values = c("grey", "red2","green2")) +
					geom_point(aes(color = Significance , size = unique.peptides)) +
					scale_size(range = c(1,3), breaks = c(1, 2, 3), labels=c("1","2",">=3"), guide = "legend") +
					theme_bw(base_size = 20) +
					geom_hline(yintercept = -log10(pvalcut), colour="grey") +
					geom_vline(xintercept = c(-FCcut,0,FCcut), colour="grey") +
					ylab(ylab) + xlab("log2 Fold Change") +
					ggtitle(sel_test) +
					theme(legend.position = input$vlegendpos, legend.text=element_text(size=input$yfontsize)) +
					guides(color = guide_legend(override.aes = list(alpha = 1, size = 4)), size = "legend")
					pl <- p

					if (input$label!="None") {
						p = p + geom_text_repel(data = data.label, aes(label=labelid),	size = input$lfontsize,	box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines") )
					}
				} else {
					p <- p +
					theme_bw(base_size = 20) +
					geom_hline(yintercept = -log10(pvalcut), colour="grey") +
					geom_vline(xintercept = c(-FCcut,0,FCcut), colour="grey") +
					ylab(ylab) + xlab("log2 Fold Change") +
					ggtitle(sel_test) +
					theme(legend.position = input$vlegendpos, legend.text=element_text(size=input$yfontsize))
					p <- p + guides(color = guide_legend(override.aes = list(alpha = 1, size = 4)))
					pl <- p

					if (input$label!="None") {
						p = p + geom_text_repel(data = data.label, aes(label=labelid),	size = input$lfontsize,	box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines") )
					}

				}
				return(list(pl = pl, p = p))
			})

			output$volcanoplotstatic <- renderPlot({
				volcanoplot_out()$p
			})

			output$volcanoplotinteractive <- renderPlotly({
				pl <-	volcanoplot_out()$pl
				ggplotly((pl + theme_bw(base_size = 16)),  tooltip = c("text")) %>% layout(legend = list(orientation = "h", y = -0.2))
			})

			output$plot.volcano <- renderUI({
				if (input$interactive == "static"){
					plotOutput(ns("volcanoplotstatic"), height=800)
				}	else {
					plotly::plotlyOutput(ns("volcanoplotinteractive"), height=800)
				}
			})

			observeEvent(input$volcano, {
				sel_test = input$test
				saved_plots$volcano[[sel_test]] <- volcanoplot_out()$p
			})

			DEG_data <- reactive ({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$results_long)
				results_long = DataInSets[[working_project()]]$results_long

				req(input$psel)
				p_sel   <- input$psel
				test_sel <- input$test
				FCcut <- log2(as.numeric(input$FCcut))
				pvalcut <- as.numeric(input$pvalcut)
				sel_label <- "UniqueID"
				direction <- "UpDown"
				tmpdat <- GeneFilter(results_long, test_sel, p_sel, direction, pvalcut, FCcut, sel_label)
				tmpdat[,sapply(tmpdat,is.numeric)] <- signif(tmpdat[,sapply(tmpdat,is.numeric)],3)
				return(tmpdat)
			})

			output$volcanoData <- DT::renderDataTable(server=FALSE,{
				DT::datatable(DEG_data() %>% dplyr::select(-labelid), extensions = 'Buttons',  selection = 'none', class = 'cell-border strip hover',
					options = list(
						dom = 'lBfrtip', pageLength = 20,
						buttons = list(
							list(extend = "csv", text = "Download Page", filename = "Page_results",
								exportOptions = list(modifier = list(page = "current"))
							),
							list(extend = "csv", text = "Download All", filename = "All_Results",
								exportOptions = list(modifier = list(page = "all"))
							)
						)
					),
				rownames= FALSE)
			})

			observeEvent(input$DEG_data, {
				sel_test <- input$test
				DEG_data_tab <- paste("DEG_data_",sel_test,sep="")
				saved_table[[DEG_data_tab]] <- DEG_data()
			})

			deg_counts_data <- reactive ({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$results_long)
				results_long <-  DataInSets[[working_project()]]$results_long

				p_sel   <- input$psel
				FCcut <- log2(as.numeric(input$FCcut))
				pvalcut <- as.numeric(input$pvalcut)

				sel_label <- "UniqueID"
				direction <- "UpDown"
				test_sel <- NA
				tmpdat <- GeneFilter(results_long, test_sel, p_sel, direction, pvalcut, FCcut, sel_label)

				deg_stat <- tmpdat %>%
				dplyr::filter(test %in% DataInSets[[working_project()]]$tests_order) %>%
				dplyr::group_by(test) %>%
				dplyr::summarize(DEG=n(), Up=sum(logFC>0), Down=sum(logFC<0)) %>%
				dplyr::ungroup()

				return(deg_stat)
			})

			output$deg_counts <- DT::renderDT(server=FALSE, {
				DT::datatable(deg_counts_data(), extensions = 'Buttons',  selection = 'none', class = 'cell-border strip hover',
					options = list(
						dom = 'lBfrtip', pageLength = 20,
						buttons = list(
							list(extend = "csv", text = "Download Page", filename = "Page_results",
								exportOptions = list(modifier = list(page = "current"))
							),
							list(extend = "csv", text = "Download All", filename = "All_Results",
								exportOptions = list(modifier = list(page = "all"))
							)
						)
					),
				rownames= FALSE) %>%
				formatStyle(1, cursor = 'pointer', color='blue')
			})

			observeEvent(input$deg_counts_cell_clicked, {
				info = input$deg_counts_cell_clicked
				if (is.null(info$value) || info$col != 0) return()
				updateTabsetPanel(session, 'tabset', selected = 'Volcano Plot')
				updateSelectizeInput(session, 'test', selected = info$value)
			})

			###Multivolcano
			DataMultivolcanoReactive <- reactive({
				req(length(working_project()) > 0)
				req(input$add_dataset)
				selectdatasets <- input$add_dataset
				results_long <- DataInSets[[working_project()]]$results_long

				FCcut <- log2(as.numeric(input$FCcut))
				FCcut_rd <- round(FCcut*1000)/1000
				psel <- input$psel
				pvalcut <- as.numeric(input$pvalcut)
				genelabel <- input$genelabel
				Max_Pvalue <- input$Max_Pvalue
				Max_logFC <- input$Max_logFC

				if (psel == "Pval")
				psel = "P.Value"
				if (psel == "Padj")
				psel = "Adj.P.Value"

				selectdatasetsdf <- stringr::str_split(selectdatasets, "->", simplify = TRUE) %>%
				as.data.frame() %>%
				rlang::set_names(c("project", "tests")) %>%
				dplyr::group_by(project) %>%
				dplyr::summarise(tests = list(tests))

				results_long_sel <- list()
				for (row in 1:nrow(selectdatasetsdf)) {
					ProjectID <- selectdatasetsdf[row, "project"] %>% as.character()
					ShortName <- DataInSets[[ProjectID]]$ShortName
					tests <- selectdatasetsdf[row, "tests"][[1]][[1]]
					res = DataInSets[[ProjectID]]$results_long %>%
					dplyr::filter(test %in% tests) %>%
					dplyr::select(one_of(c("UniqueID","Gene.Name", "Protein.ID", "test","logFC", "P.Value", "Adj.P.Value")) ) %>%
					dplyr::mutate (ProjectID = ShortName)
					results_long_sel[[row]] <- res
				}

				res <- dplyr::bind_rows(results_long_sel) %>%
				dplyr::mutate_if(is.factor, as.character) %>%
				dplyr::mutate(labelid = !!sym(genelabel)) %>%
				dplyr::filter(!is.na(P.Value)) %>%
				dplyr::mutate(
					Significance = case_when(
						abs(logFC) >= FCcut & !!sym(psel) <= pvalcut ~ paste0(psel,"<", pvalcut," & abs(log2FC)>=", FCcut_rd),
						abs(logFC) <  FCcut & !!sym(psel) <= pvalcut ~ paste0(psel,"<", pvalcut," & abs(log2FC)<", FCcut_rd),
					TRUE ~ "Not Significant")
				) %>%
				dplyr::mutate(Significance = factor(Significance,	levels = c("Not Significant",	paste0(psel,"<", pvalcut," & abs(log2FC)>=", FCcut_rd),	paste0(psel,"<", pvalcut," & abs(log2FC)<", FCcut_rd)))
				) %>%
				dplyr::mutate(logFC_ori = logFC)

				if (Max_Pvalue > 0) {
					if (psel == "Adj.P.Value")
					res <- res %>% dplyr::mutate(Adj.P.Value=pmax(Adj.P.Value, 10^(0-Max_Pvalue) ))
					if (psel == "P.Value")
					res <- res %>% dplyr::mutate(P.Value=pmax(P.Value, 10^(0-Max_Pvalue) ))
				}

				if (Max_logFC > 0) {
					res <- res %>% dplyr::mutate(logFC=ifelse(logFC>=0, pmin(Max_logFC, logFC), pmax(0-Max_logFC, logFC)))
				}
				return(res)
			})

			Multivolcanoplot_out <- eventReactive(input$plot_Multivolcano, {
				req(DataInSets[[working_project()]]$ProteinGeneName)
				ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName

				res = DataMultivolcanoReactive()
				res <- res %>%
				tidyr::unite(project_test, ProjectID, 'test', sep=".", remove=FALSE)

				FCcut <- log2(as.numeric(input$FCcut))
				FCcut_rd <- round(FCcut*1000)/1000
				pvalcut <- as.numeric(input$pvalcut)
				psel <- input$psel

				Ngenes <- as.numeric(input$Ngenes)
				vlegendpos <- input$vlegendpos
				yfontsize <- as.numeric(input$yfontsize)

				res$logFC_ori=res$logFC
				if (psel == "Pval")
				psel = "P.Value"
				if (psel == "Padj")
				psel = "Adj.P.Value"

				if (psel == "Adj.P.Value") {
					p <- ggplot(res, aes(x = logFC, y = -log10(Adj.P.Value), text = labelid))
					ylab <- "-log10(Padj.Value)"

					filterSig <- paste0(psel,"<", pvalcut," & abs(log2FC)>=", FCcut_rd)
					data.label <- filter(res, Significance == filterSig)

					if (nrow(data.label) > Ngenes) {
						data.label <- top_n(data.label, Ngenes, abs(logFC_ori))
					}

				} else {
					filterSig <- paste0(psel,"<", pvalcut," & abs(log2FC)>=", FCcut_rd)
					data.label <- filter(res, Significance == filterSig)
					if (nrow(data.label) > Ngenes) {
						data.label <- top_n(data.label, Ngenes, abs(logFC_ori))
					}
					p <- ggplot(res, aes(x = logFC, y = -log10(P.Value), text = labelid))
					ylab <- "-log10(P.Value)"
				}

				if (input$label=="Upload" | input$label=="Geneset") {
					if (input$label=="Upload") {
						req(input$gene_list)
						gene_list <- input$gene_list
					} else {
						req(input$geneset_genes)
						gene_list <- input$geneset_genes
					}

					if(grepl("\n",gene_list)) {
						gene_list <-  stringr::str_split(gene_list, "\n")[[1]]
					} else if(grepl(",",gene_list)) {
						gene_list <-  stringr::str_split(gene_list, ",")[[1]]
					}
					gene_list <- gsub(" ", "", gene_list, fixed = TRUE)
					gene_list  <- unique(gene_list[gene_list != ""])

					validate(need(length(gene_list)>0, message = "input gene list"))

					uploadlist <- dplyr::filter(ProteinGeneName, (UniqueID %in% gene_list) | (Protein.ID %in% gene_list) | (Gene.Name %in% gene_list))  %>%
					dplyr::pull(UniqueID)

					data.label <- res %>%
					dplyr::filter(UniqueID %in% uploadlist)
					validate(need(nrow(data.label)>0, message = "no gene found in result"))
				}

				OnePlotColorBy = input$OnePlotColorBy

				if (input$rasterize=="Yes") {
					p <- p + geom_point_rast(aes(color = !!sym(OnePlotColorBy)), size=0.7, alpha=0.6, na.rm=TRUE, dev="ragg")
				} else {
					p <- p + geom_point(aes(color = !!sym(OnePlotColorBy)), size=0.7, alpha=0.6, na.rm=TRUE)
				}

				if (input$SeparateOnePlot == "Separate")
				p <- p + 	facet_grid(~ ProjectID + test)

				p <- p + theme_bw(base_size = 20) +
				geom_hline(yintercept = -log10(pvalcut), colour="grey") +
				geom_vline(xintercept = c(-FCcut,0,FCcut), colour="grey") +
				ylab(ylab) + xlab("log2 Fold Change") +
				theme(legend.position = vlegendpos, legend.text=element_text(size=yfontsize)) +
				guides(color = guide_legend(title = "Comparison", override.aes = list(alpha = 1, size = 4)))
				pl <- p

				if (input$label!="None") {
					p <- p + geom_text_repel(data = data.label,	aes(x = logFC, y = -log10(P.Value),label=labelid),	size = 5,	box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines") )
					if (input$SeparateOnePlot == "OnePlotOnePlot")
					p <- p +  geom_text_repel(data = data.label, aes(x = logFC, y = -log10(P.Value), label = labelid),	size = 5,	box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines") )
				}

				return(list(pl = pl, p = p))
			})

			output$Multivolcano_Static <- renderPlot({
				Multivolcanoplot_out()$p %>% print
			})

			output$Multivolcano_Interactive <- renderPlotly({
				pl = Multivolcanoplot_out()$pl
				ggplotly((pl + theme_bw(base_size = 16)),  tooltip = c("text")) %>% layout(legend = list(orientation = "h", y = -0.2))
			})

			output$plot.Multivolcano <- renderUI({
				if (input$interactive == "static"){
					plotOutput(ns("Multivolcano_Static"), height=800)
				}	else {
					plotly::plotlyOutput(ns("Multivolcano_Interactive"),	height=800)
				}
			})

			observeEvent(input$Multivolcano, {
				sel_test = "Multivolcano"
				saved_plots$volcano[[sel_test]] <- Multivolcanoplot_out()$p
			})

			####DEG_Compare
			DataDEGCompareReactive <- reactive({
				req(length(working_project()) > 0)
				selectdatasets <- input$add_dataset
				req(length(selectdatasets)==2)
				FCcut = log2(as.numeric(input$FCcut))
				FCcut_rd=round(FCcut*1000)/1000
				pvalcut = as.numeric(input$pvalcut)
				genelabel = input$genelabel

				selectdatasetsdf <- stringr::str_split(selectdatasets, "->", simplify = TRUE) %>%
				as.data.frame() %>%
				rlang::set_names(c("project", "tests"))

				for(i in 1:nrow(selectdatasetsdf)) {
					if (i==1) {
						notsigstr <- "X_notsig"
						sigstr <-  "X_sig"
					}  else {
						notsigstr <- "Y_notsig"
						sigstr <-  "Y_sig"
					}

					ProjectID <- selectdatasetsdf[i, "project"] %>% as.character()
					tests <- selectdatasetsdf[i, "tests"]
					sel_test = paste(DataInSets[[ProjectID]]$ShortName,"-",tests,sep="")
					assign(paste0("sel_test", i), sel_test)

					res = DataInSets[[ProjectID]]$results_long %>%
					dplyr::filter(test == tests) %>%
					dplyr::mutate (ProjectID = ProjectID) %>%
					dplyr::filter(!is.na(P.Value)) %>%
					dplyr::mutate (color = "Not Significant") %>%
					dplyr::mutate (Sig = notsigstr) %>%
					as.data.frame()

					res$labelid = res[,match(genelabel, colnames(res))]

					if (input$psel == "Padj") {
						res$Sig[which((abs(res$logFC)>FCcut)*(res$Adj.P.Value<pvalcut)==1)] = sigstr

						if (input$Max_Pvalue>0) {
							res <- res %>% mutate(Adj.P.Value=pmax(Adj.P.Value, 10^(0-input$Max_Pvalue) ))
						}
					} else {
						res$Sig[which((abs(res$logFC)>FCcut)*(res$P.Value<pvalcut)==1)] = sigstr
						if (input$Max_Pvalue>0) {
							res <- res %>% mutate(P.Value=pmax(P.Value, 10^(0-input$Max_Pvalue) ))
						}
					}

					if (input$Max_logFC>0) {
						res <- res %>% mutate(logFC=ifelse(logFC>=0, pmin(input$Max_logFC, logFC), pmax(0-input$Max_logFC, logFC) ) )
					}
					assign(paste0("res", i), res)
				}

				merged_by <- input$merged_by
				res = dplyr::inner_join(res1, res2, by = merged_by, relationship = "many-to-many")
				validate(need(nrow(res) > 1, message = "Length of merged data is 0."))
				FC_Cutoff <- log2(as.numeric(input$FC_Cutoff))

				res <- res %>% dplyr::mutate(color = paste(Sig.x, Sig.y)) %>%
				dplyr::select(-c(Sig.x, Sig.y)) %>%
				dplyr::filter(abs(logFC.x) > FC_Cutoff | abs(logFC.y) > FC_Cutoff)
				c.res <- coop::pcor(res$logFC.x, res$logFC.y, use = "complete.obs")
				cor_string <- paste("Pearson corr:", format(c.res, digits=3), sep="")

				if (input$LM == "Yes") {
					model_equation <- function(model, ...) {
						format_args <- list(...)

						model_coeff <- model$coefficients
						format_args$x <- abs(model$coefficients)
						model_coeff_sign <- sign(model_coeff)
						model_coeff_prefix <- case_when(model_coeff_sign == -1 ~ " - ",	model_coeff_sign == 1 ~ " + ", model_coeff_sign == 0 ~ " + ")

						model_eqn <- paste(strsplit(as.character(model$call$formula), "~")[[2]], "=",
							paste(if_else(model_coeff[1]<0, "- ", ""),
								do.call(format, format_args)[1],
								paste(model_coeff_prefix[-1],
									do.call(format, format_args)[-1],
									" * ",
									names(model_coeff[-1]),
								sep = "", collapse = ""),
							sep = "")
						)
						return(model_eqn)
					}

					model <- lm(logFC.y ~ logFC.x, data = res)
					formula.str <- model_equation(model, digits = 2)
					cor_string <- paste("Pearson corr:", format(c.res, digits=3), " Formula:", formula.str, sep="")
				} else {
					model = NA
				}
				return(list("res" = res, "model" = model, "cor_string" = cor_string, "sel_test1" = sel_test1, "sel_test2" = sel_test2))
			})

			DEG_Compare <- eventReactive(input$plot_DEG_comp, {
				withProgress(message = 'Making Plot. It may take a while...', value = 0, {

					req(length(working_project()) > 0)
					validate(need(length(input$add_dataset) == 2, message = "Please Select 2 sets"))

					rescompare = DataDEGCompareReactive()
					ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
					res <- rescompare$res

					cor_string <- rescompare$cor_string
					sel_test1 = rescompare$sel_test1
					sel_test2 = rescompare$sel_test2

					if (input$Xlab != "")
					sel_test1 = input$Xlab
					if (input$Ylab != "")
					sel_test2 = input$Ylab

					FCcut = log2(as.numeric(input$FCcut))
					pvalcut = as.numeric(input$pvalcut)

					data.label <- res %>%
					dplyr::filter(str_detect(color, "_sig")) %>%
					dplyr::mutate(H_logFC=pmax(abs(logFC.x), abs(logFC.y)))  #at least one is sig

					if (input$psel == "Padj") {
						res <- res %>%
						dplyr::mutate(size=-pmin(log10(Adj.P.Value.x),log10(Adj.P.Value.y)))

						p <- ggplot(res, aes(x=logFC.x, y=logFC.y, color=color, size=size, text=labelid.x))

						fitresult <- data.frame()
						if (input$LM == "Yes") {
							#get intercept and slope value
							model <- rescompare$model
							coeff <- coefficients(model)
							intercept <- coeff[1]
							slope <- coeff[2]
							p <- p + geom_abline(intercept = intercept, slope = slope, color="red",  size=1)

							fitresult <- as.data.frame(summary(model)$coefficients) %>%
							mutate_if(is.numeric, round, digits = 4)
						}

						if (input$rasterize=="Yes") {
							p <- p + geom_point_rast(na.rm=TRUE, dev="ragg")
						} else {
							p <- p + geom_point()
						}

						p <- p +
						theme_bw(base_size = 16) + ylab(str_c("log2FC in ", sel_test2)) + xlab(str_c("log2FC in ", sel_test1))
						p <- p +
						labs(color='Significance', size='-log10 min Adj.P.Value', title=cor_string)
					} else {
						res <- res %>%
						dplyr::mutate(size=-pmin(log10(P.Value.x),log10(P.Value.y))) 

						p <- ggplot(res, aes(x=logFC.x, y=logFC.y, color=color,	size=size, text=labelid.x))
						fitresult <- data.frame()
						if (input$LM == "Yes") {
							#get intercept and slope value
							model <- rescompare$model
							coeff <- coefficients(model)
							intercept <- coeff[1]
							slope <- coeff[2]
							p <- p + geom_abline(intercept = intercept, slope = slope, color="red",  size=1)

							fitresult <- as.data.frame(summary(model)$coefficients) %>%
							mutate_if(is.numeric, round, digits = 4)
						}

						if (input$rasterize=="Yes") {
							p <- p + geom_point_rast(na.rm=TRUE, dev="ragg")
						} else {
							p <- p + geom_point()
						}
						p <- p +
						theme_bw(base_size = 16) + ylab(str_c("log2FC in ", sel_test2)) + xlab(str_c("log2FC in ", sel_test1))
						p <- p +
						labs(color='Significance', size='-log10 min P.Value', title=cor_string)
					}

					if (input$label=="Upload" | input$label=="Geneset") {
						if (input$label=="Upload") {
							req(input$gene_list)
							gene_list <- input$gene_list
						} else {
							req(input$geneset_genes)
							gene_list <- input$geneset_genes
						}

						gene_list <- ProcessUploadGeneList(gene_list)

						validate(need(length(gene_list)>0, message = "input gene list"))

						uploadlist <- dplyr::filter(ProteinGeneName, (UniqueID %in% gene_list) | (Protein.ID %in% gene_list) | (Gene.Name %in% gene_list))  %>%
						dplyr::pull(UniqueID)

						data.label <- res %>%
						dplyr::filter(if_any(starts_with("UniqueID"), ~ . %in% uploadlist))
						validate(need(nrow(data.label)>0, message = "no gene found in result"))
					}

					p <- p + scale_color_manual(values=c('X_sig Y_sig'='blue3','X_sig Y_notsig'='green3', 'X_notsig Y_sig'='orange','X_notsig Y_notsig'='#00000022')) +
					theme(legend.position=input$vlegendpos, legend.text=element_text(size=input$yfontsize), legend.title=element_text(size=input$yfontsize+1))

					if (input$DEG_comp_XY=="Yes"){
						XY_min <- min(min(res$logFC.x), min(res$logFC.y))
						XY_max <- max(max(res$logFC.x), max(res$logFC.y))
						p <- p + xlim(XY_min, XY_max) + ylim(XY_min, XY_max)
					}

					p <- p + guides(color = guide_legend(override.aes = list(alpha = 1, size = 4))) # Remove strange "a" from legend
					pl <- p

					if (nrow(data.label) > input$Ngenes) {
						data.label <- data.label %>%
						dplyr::mutate(H_logFC=pmax(abs(logFC.x), abs(logFC.y)))  %>%
						top_n(input$Ngenes, H_logFC)
					}

					if (input$label=="Upload" || input$label=="Geneset" || input$DEG_comp_color=="No") {
						p <- p +
						geom_text_repel(data = data.label,  aes(label=labelid.x),	size = input$lfontsize,	box.padding = unit(0.35, "lines"), color="coral3",  point.padding = unit(0.3, "lines"))
					} else {
						p <- p +
						geom_text_repel(data = data.label, aes(label=labelid.x),	size = input$lfontsize,	box.padding = unit(0.35, "lines"),	point.padding = unit(0.3, "lines"))
					}

					return(list(pl = pl, p = p, fitresult = fitresult))
				})
			})

			output$DEG_Compare_Static <- renderPlot({
				DEG_Compare()$p %>% print
			})

			output$DEG_Compare_Interactive <- renderPlotly({
				pl = DEG_Compare()$pl
				pl$labels$size <- NULL
				ggplotly((pl + theme_bw(base_size = 16)),  tooltip = c("text")) %>% layout(legend = list(orientation = "h", y = -0.2))
			})

			output$plot.DEG_Compare <- renderUI({
				if (input$interactive == "static"){
					plotOutput(ns("DEG_Compare_Static"), height=800)
				}	else {
					plotly::plotlyOutput(ns("DEG_Compare_Interactive"),	height=800)
				}
			})

			output$LMResult <- renderDataTable({
				fitresult = DEG_Compare()$fitresult
				DT::datatable(fitresult)
			})

			output$DEGCompareData <- DT::renderDataTable(server=TRUE,{
				req(length(working_project()) > 0)
				validate(need(length(input$add_dataset) == 2, message = "Please Select 2 sets"))
				rescompare = DataDEGCompareReactive()

				DEGCompareRes <- rescompare$res %>%
				dplyr::select(-any_of(c("id.x", "color.x","labelid.x","id.y", "color.y",	"labelid.y"))) 	%>%
				dplyr::mutate_if(is.numeric, round, digits = 4)

				DT::datatable(DEGCompareRes, extensions = 'Buttons',  selection = 'none', class = 'cell-border strip hover',
					options = list(
						dom = 'lBfrtip', pageLength = 20,
						buttons = list(
							list(extend = "csv", text = "Download Page", filename = "Page_results",
								exportOptions = list(modifier = list(page = "current"))
							)
						)
					),
				rownames= FALSE)
			})

			output$download_data_button <- shiny::downloadHandler(
				filename <- function() {
					paste("Data-", Sys.Date(), ".csv", sep="")
				},
				content = function(file) {
					req(length(working_project()) > 0)
					validate(need(length(input$add_dataset) == 2, message = "Please Select 2 sets"))
					rescompare = DataDEGCompareReactive()
					DEGCompareRes <- rescompare$res %>%
					dplyr::select(-any_of(c("id.x", "color.x","labelid.x","id.y", "color.y",	"labelid.y"))) 	%>%
					dplyr::mutate_if(is.numeric, round, digits = 4)
					write.csv(DEGCompareRes, file)
				}
			)

			observeEvent(input$DEGCompareData, {
				req(length(working_project()) > 0)
				validate(need(length(input$add_dataset) == 2, message = "Please Select 2 sets"))
				rescompare = DataDEGCompareReactive()
				DEGCompareRes <- rescompare$res %>%
				dplyr::select(-any_of(c("id.x", "color.x","labelid.x","id.y", "color.y",	"labelid.y"))) 	%>%
				dplyr::mutate_if(is.numeric, round, digits = 4)

				saved_table$DEGCompareData <- DEGCompareRes
			})

			observeEvent(input$DEG_comp, {
				selectdatasets <- input$add_dataset
				selectdatasetsdf <- stringr::str_split(selectdatasets, "->", simplify = TRUE) %>%
				as.data.frame() %>%	rlang::set_names(c("project", "tests"))
				row <- 1
				ProjectID1 <- selectdatasetsdf[row, "project"] %>% as.character()
				tests1 <- selectdatasetsdf[row, "tests"]
				row <- 2
				ProjectID2 <- selectdatasetsdf[row, "project"] %>% as.character()
				tests2 <- selectdatasetsdf[row, "tests"]

				sel_test = paste(ProjectID1,"_",tests1, "vs", ProjectID2,"_",tests2)
				saved_plots$volcano[[sel_test]] <- DEG_Compare()$p
			})
		}
	)
}
