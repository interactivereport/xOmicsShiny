###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: expressionmodule.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 08/28/2024
##@version 3.0

##########################################################################################################
## Gene/Protein Expression Plot
##########################################################################################################
#pkgs: "ggprism", "DT","dplyr", "stringr", "tidyr", "scales", "ggExtra", "png", "grid"

library(ggprism)

expression_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
		column(4,
			wellPanel(
				uiOutput(ns('loadedprojects')),
				uiOutput(ns("selectGroupSample")),
				conditionalPanel(ns = ns, "input.tabset !='Help'",
					radioButtons(ns("subset"), label="Genes Used in Plot", choices=c("Select", "Browsing", "Upload Genes", "Geneset"), inline = TRUE, selected="Select"),
					conditionalPanel(ns = ns, "input.subset=='Upload Genes'",
						textAreaInput(ns("uploadlist"), "Enter Gene List", "", cols = 5, rows=6)
					)
				),
				conditionalPanel(ns = ns, "input.subset=='Select'",
					radioButtons(ns("label"),label="Search Gene or UniqueID", inline = TRUE, choices=c("UniqueID", "Gene.Name"), selected="UniqueID"),
					selectizeInput(ns("sel_gene"),	label="Gene Name (Select 1 or more)",	choices = NULL,	multiple=TRUE, options = list(placeholder =	'Type to search'))
				),
				conditionalPanel(ns = ns, "input.subset=='Geneset'",
					selectizeInput(ns("sel_geneset"), label="Available GeneSet", choices = NULL, multiple = FALSE),
					textAreaInput(ns("geneset_genes"), "Genes in Geneset", "", cols = 5, rows=6)
				),
				conditionalPanel(ns = ns, "input.subset=='Browsing'",
					selectInput(ns("sel_test"), label="Select Test", choices=NULL),
					fluidRow(
						column(width=6, radioButtons(ns("psel"), label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"), inline = TRUE)),
						column(width=6, radioButtons(ns("updown"), label= "All, Up or Down?", choices= c("All"="All","Up"="Up","Down"="Down"), inline = TRUE))
					),
					fluidRow(
						column(width=6, numericInput(ns("fccut"), label= "Fold Change Threshold", value = 1.2, min=1, step=0.1)),
						column(width=6, numericInput(ns("pvalcut"), label= "P-value Threshold", value=0.01, min=0, step=0.01))
					),
					radioButtons(ns("browsing_gene_order"), label="Order Genes by", inline = TRUE, choices = c("abs(Fold Change)","P value"))
				),
				radioButtons(ns("sel_geneid"), label="Select Gene Label", inline = TRUE, choices=""),
				span(textOutput(ns("filteredgene")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
				span(textOutput(ns("filteredUniqueID")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
				tags$hr(style="border-color: black;"),
				conditionalPanel(ns = ns, "input.tabset=='Expression Plot'",
					radioButtons(ns("SeparateOnePlot"), label="Separate or One Plot", inline = TRUE, choices = c("Separate Plot" = "Separate", "One Plot" = "OnePlot")),
					radioButtons(ns("plotformat"), label="Select Plot Format", inline = TRUE, choices = c("Box Plot" = "boxplot","Bar Plot" = "barplot", "Violin" = "violin","Line" = "line")),
					radioButtons(ns("prismstyle"), label="Graphpad Prism Style", inline = TRUE, choices = c("YES" = "YES","NO" = "NO"), selected = "NO"),
					conditionalPanel(ns = ns, "input.prismstyle=='NO'",
						selectInput(ns("colpalette"), label= "Select palette", choices=""),
						conditionalPanel(ns = ns, "input.colpalette=='Single'",	colourpicker::colourInput(ns("barcol"), "Select colour", "#1E90FF", palette = "limited")),
					),
					conditionalPanel(ns = ns, "input.prismstyle=='YES'",
						fluidRow(
							column(width=6,selectInput(ns("prismcolpalette"), label= "Select palette", choices="")),
							column(width=6,selectInput(ns("prismfillpalette"), label= "Select palette", choices=""))
						)
					)
				),
				conditionalPanel(ns = ns, "input.tabset=='Expression Plot' | input.tabset=='Laser Capture Microdissection'",
					fluidRow(
						column(width=4,sliderInput(ns("plot_ncol"), label= "Column Number", min = 1, max = 6, step = 1, value = 3)),
						column(width=4,sliderInput(ns("plot_nrow"), label= "Row Number", min = 1, max = 6, step = 1, value = 3)),
						column(width=4,selectInput(ns("sel_page"), label="Select Page",	choices = NULL,	selected=1))
					)
				),
				conditionalPanel(ns = ns, "input.tabset=='Expression Plot'",
					radioButtons(ns("IndividualPoint"), label="Individual Point?", inline = TRUE, choices = c("YES" = "YES","NO" = "NO"), selected = "YES"),
					radioButtons(ns("ShowErrorBar"), label="Error Bar?(except box plot)", inline = TRUE, choices = c("SD" = "SD","SEM" = "SEM","NO" = "NO"), selected = "SD"),
					radioButtons(ns("PvalueBar"), label="Show P-Values?", inline = TRUE, choices = c("YES" = "YES","NO" = "NO"), selected = "NO"),
					conditionalPanel(ns = ns, "input.PvalueBar=='YES'",
						fluidRow(
							column(width=6,radioButtons(inputId = ns("pval_sel"), label = "P-value or P.adj Value?",inline = TRUE, choices = c("Pval" = "P.Value", "Padj" = "Adj.P.Value"), selected = "P.Value")),
							column(width=6,numericInput(ns("pv_threshold"), label= "P Value Threshold:",  value= 0.01, min=0, step=0.01))
						),
						radioButtons(inputId = ns("plabelpos"), label = "P Value Position?", inline = TRUE, choices = c("Top" = "y2", "inside" = "y1","Top 3 SD" = "y3"), selected = "y2")
					),
					fluidRow(
						column(width=6, selectizeInput(ns("plotx"), label="X Axis", choices=NULL, multiple = FALSE)),
						column(width=6, selectizeInput(ns("colorby"), label="Color By:", choices=NULL, multiple = FALSE))
					),
					fluidRow(
						column(width=6, radioButtons(ns("xlogbase"), label="X Log tansform", inline = TRUE,  choices=c("Non"= 1, "10"= 10, "2"= 2), selected=1)),
						column(width=6,	radioButtons(ns("ylogbase"), label="Y Log tansform", inline = TRUE,  choices=c("Non"= 1, "10"= 10, "2"= 2), selected=1))
					),
					fluidRow(
						column(width=6,sliderInput(ns("titlefontsize"), "Title Font Size:", min = 12, max = 36, step = 2, value = 20)),
						column(width=6,sliderInput(ns("axisfontsize"), "Axis Font Size:", min = 12, max = 36, step = 2, value = 20))
					),
					fluidRow(
						column(width=6,textInput(ns("Ylab"), "Y label", width = "100%")),
						column(width=6,textInput(ns("Xlab"), "X label", width = "100%"))
					),
					sliderInput(ns("Xangle"), label= "X Angle", min = 0, max = 90, step = 5, value = 90),
					radioButtons(ns("plot_Y_scale"), label="Y Axis Scale", inline = TRUE, choices = c("Auto","Manual"), selected = "Auto"),
					radioButtons(ns("interactive"), label="Static or Interactive Plot", inline = TRUE, choices = c("Static" = "static", "Interactive" = "interactive")),
					conditionalPanel(ns = ns,"input.plot_Y_scale=='Manual'",
						fluidRow(
							column(width=6,numericInput(ns("plot_Ymin"), label= "Y Min",  value = 0, step=0.1)),
							column(width=6,numericInput(ns("plot_Ymax"), label= "Y Max",  value=5, step=0.1))
						)
					),
					conditionalPanel(ns = ns, "input.tabset=='Expression Plot'",
						h5("After changing parameters, please click Plot/Refresh button in the plot panel to generate expression plot.")
					)
				),
				conditionalPanel(ns = ns, "input.tabset=='Rank Abundance Curve'",
					fluidRow(
						column(width=6, sliderInput(ns("scurve_axisfontsize"), "Axis Font Size:", min = 12, max = 28, step = 4, value = 16)),
						column(width=6, sliderInput(ns("scurve_labelfontsize"), "Label Font Size:", min = 2, max = 12, step = 1, value = 6))
					),
					fluidRow(
						column(width=6, textInput(ns("scurveYlab"), "Y label", value="Abundance", width = "100%")),
						column(width=6, textInput(ns("scurveXlab"), "X label", value="Rank", width = "100%"))
					),
					sliderInput(ns("scurveXangle"), label= "X Angle", min = 0, max = 90, step = 15, value = 45),
					radioButtons(ns("scurveright"), label="Density or histogram on Right:", inline = TRUE, choices = c("densigram" = "densigram", "density" = "density","histogram" = "histogram","boxplot" = "boxplot", "violin"= "violin"), selected = "densigram")
				),
				conditionalPanel(ns = ns, "input.tabset=='Expression Correlation'",
					conditionalPanel(ns = ns, condition = "input.sel_gene.length != 2",
						tags$h4("Visualization of a correlation matrix options:"),
						fluidRow(
							column(width=6, sliderInput(ns("tlcex"), "Label Font Size:", min = 0.4, max = 2, step = 0.1, value = 0.6)),
							column(width=6, numericInput(ns("siglevel"), label= "sig.level:", value= 0.05, min=0, max = 1, step=0.01))
						),
						radioButtons(ns("method"), label="Visualization Method:", inline = TRUE, choices = c("circle" = "circle", "square" = "square", "ellipse" = "ellipse", "number" = "number", "shade" = "shade", "pie" = "pie"), selected = "circle"),
						radioButtons(ns("corrcol"), label="Color:", inline = TRUE, choices = c("RdBu" = "RdBu", "BrBG" = "BrBG", "PiYG" = "PiYG", "PRGn" = "PRGn", "PuOr" = "PuOr", "RdYlBu" = "RdYlBu"), selected = "RdBu"),
						radioButtons(ns("type"),  label="Plot Type", inline = TRUE, choices = c("full"="full", "lower"="lower", "upper"="upper"), selected = "upper"),
						radioButtons(ns("order"),  label="Order by", inline = TRUE, choices = c("original"="original", "Angular order of eigenvectors" = "AOE", "First principal component"="FPC", "Hierarchical clustering"="hclust", "alphabet"="alphabet"),selected = "AOE"),
						radioButtons(ns("hclustmethod"), label="Cluster Method", inline = TRUE, choices = c("complete"="complete", "ward"="ward", "ward.D"="ward.D", "ward.D2"="ward.D2", "single"="single", "average"="average","mcquitty"="mcquitty", "median"="median", "centroid"="centroid"),selected = "complete"),
						radioButtons(ns("diag"), label="Show Principal Diagonal", inline = TRUE, choices = c("YES" = "YES","NO" = "NO"), selected = "YES")
					)
				),
				conditionalPanel(ns = ns, "input.tabset=='Laser Capture Microdissection'",
					selectizeInput(ns("sel_image"), label="Select Images",	choices = NULL,	multiple=FALSE, options = list(placeholder =	'Type to search')),
					radioButtons(ns("convert2ratio"), label="Convert to Ratio to Average", inline = TRUE, choices=c("expr", "ratio"), selected="expr"),
					conditionalPanel(ns = ns,"input.convert2ratio=='ratio'",
						fluidRow(
							column(width=6,numericInput(ns("LowLimit"), label= "Low Limit",  value = -5, step=1)),
							column(width=6,numericInput(ns("HighLimit"), label= "High Limit",  value=5, step=1))
						)
					),
					conditionalPanel(ns = ns,"input.convert2ratio=='expr'",
						selectInput(ns("lcmcolpalette"), label= "Select Color Palette", choices="")
					),
					conditionalPanel(ns = ns, condition = "input.sel_gene.length == 1",
						radioButtons(ns("lcm_format"), label="Plot Format", inline = TRUE, choices = c( "Expressoin (side by side)" = "exp_side","Group (side by side)" = "group_side", "Expression (overlay)"="exp_overlay",  "Group (overlay)" = "group_overlay"), selected = "exp_side"),
						fluidRow(
							column(width=6, sliderInput(ns("LCM_alpha"), "Transparancy", min = 0, max = 1, step = 0.1, value = 0.5)),
							column(width=6, sliderInput(ns("LCM_axisfontsize"), "Axis Font Size:", min = 12, max = 28, step = 4, value = 16))
						)
					),
					fluidRow(
						column(width=6, textInput(ns("LCMYlab"), "Y label", value="Row", width = "100%")),
						column(width=6, textInput(ns("LCMXlab"), "X label", value="Column", width = "100%"))
					)
				)
			)
		),
		column(8,
			tabsetPanel(id=ns("tabset"),
				tabPanel(title="Expression Plot", value ="Expression Plot",
					actionButton(ns("plot_exp"), "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
					actionButton(ns("saveboxplot"), "Save to output"),
					uiOutput(ns("plot.exp"))
				),
				tabPanel(title="Rank Abundance Curve", value ="Rank Abundance Curve",
					actionButton(ns("plot_SCurve"), "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
					actionButton(ns("AbundanceCurve"), "Save to output"),
					plotOutput(ns("SCurve"), height=800)
				),
				tabPanel(title="Expression Correlation", value ="Expression Correlation",
					actionButton(ns("plot_ExprCorr"), "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
					#actionButton(ns("ExprCorr"), "Save to output"),#to do
					plotOutput(ns("ExprCorr"), height=800)
				),
				tabPanel(title="Laser Capture Microdissection", value ="Laser Capture Microdissection",
					actionButton(ns("plot_LCM"), "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
					#actionButton(ns("LCM"), "Save to output"), #to do
					plotOutput(ns("LCM"), height=800)
				),
				tabPanel(title="Data Table", value ="Data Table",	DT::dataTableOutput(ns("dat_dotplot"))
				),
				tabPanel(title="Result Table", value ="Result Table",	DT::dataTableOutput(ns("res_dotplot"))
				),
				tabPanel(title="Help", value ="Help", htmlOutput("help_expression")
				)
			)
		)
	)
}

expression_server <- function(id) {
	shiny::moduleServer(id,
		function(input, output, session) {
			ns <- shiny::NS(id)
			output$loadedprojects <- renderUI({
				req(length(working_project()) > 0)
				radioButtons(ns("current_dataset"), label = "Change Working Dataset", choices=DS_names(), inline = F, selected=working_project())
			})

			observeEvent(input$current_dataset, {
				working_project(input$current_dataset)
			})

			output$selectGroupSample <-
			renderUI({
				req(length(working_project()) > 0)
				sample_info <- paste("Selected ", length(DataInSets[[working_project()]]$group_order), " out of ", length(DataInSets[[working_project()]]$groups), " Groups, ",
					length(DataInSets[[working_project()]]$sample_order), " out of ", length(DataInSets[[working_project()]]$samples),
				" Samples. (Update Selection at: Top Menu -> Groups and Samples.)", sep="")
				tagList(
					tags$p(sample_info),
					tags$hr(style="border-color: black;")
				)
			})

			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$ProteinGeneNameHeader)
				ProteinGeneNameHeader = DataInSets[[working_project()]]$ProteinGeneNameHeader
				req(DataInSets[[working_project()]]$tests_order)
				tests = DataInSets[[working_project()]]$tests_order
				updateRadioButtons(session,'sel_geneid', inline = TRUE, choices=ProteinGeneNameHeader, selected="UniqueID")
				updateSelectizeInput(session,'sel_test', choices=tests, selected=tests[1])
			})

			observe({
				req(input$prismstyle=='YES')
				updateSelectizeInput(session,'prismcolpalette', choices=names(ggprism_data$colour_palettes), selected="floral")
				updateSelectizeInput(session,'prismfillpalette', choices=names(ggprism_data$fill_palettes), selected="floral")
			})

			observe({
				req(input$prismstyle=='NO')
				updateSelectizeInput(session,'colpalette', choices=c("Single",rownames(brewer.pal.info)), selected="Dark2")
			})

			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$MetaData)
				req(DataInSets[[working_project()]]$ProteinGeneName)
				MetaData = DataInSets[[working_project()]]$MetaData
				ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
				exp_unit = DataInSets[[working_project()]]$exp_unit

				req(input$label)
				label = sym(input$label)
				DataIngenes <- ProteinGeneName %>% dplyr::pull(!!label)
				updateSelectizeInput(session,'sel_gene', choices= DataIngenes, server=TRUE)

				attributes=sort(setdiff(colnames(MetaData), c("sampleid", "Order", "ComparePairs") ))
				updateSelectInput(session, "colorby", choices=attributes, selected="group")
				updateSelectInput(session, "plotx", choices=attributes, selected="group")
				updateTextInput(session, "Ylab", value=exp_unit)
				updateTextInput(session, "Xlab", value="group")
			})

			toListen <- reactive({
				list(input$subset,input$sel_geneset)
			})

			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$results_long)
				results_long = DataInSets[[working_project()]]$results_long
				ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
				req(input$subset)
				if (input$subset == "Select") {
					req(input$sel_gene)
					sel_gene = input$sel_gene
					validate(need(length(input$sel_gene)>0,"Please select a gene."))
					ProteinGeneName_sel <- dplyr::filter(ProteinGeneName, (UniqueID %in% sel_gene) | (Protein.ID %in% sel_gene) | (toupper(Gene.Name) %in% toupper(sel_gene)))
				}

				if (input$subset == "Upload Genes"){
					req(input$uploadlist)
					gene_list <- input$uploadlist
					gene_list <- ProcessUploadGeneList(gene_list)
					validate(need(length(gene_list)>0, message = "Please input at least 1 valid gene."))
					ProteinGeneName_sel <- dplyr::filter(ProteinGeneName, (UniqueID %in% gene_list) | (Protein.ID %in% gene_list) | (toupper(Gene.Name) %in% toupper(gene_list)))
				}

				if (input$subset == "Geneset" ) {
					req(input$sel_geneset!="")
					sel_geneset <- input$sel_geneset
					gene_list <- GetGenesFromGeneSet(sel_geneset)
					validate(need(length(gene_list)>0, message = "Please input at least 1 valid gene."))
					ProteinGeneName_sel <- dplyr::filter(ProteinGeneName, (UniqueID %in% gene_list) | (Protein.ID %in% gene_list) | (toupper(Gene.Name) %in% toupper(gene_list)))
				}

				if (input$subset == "Browsing") {
					req(input$sel_test)
					req(input$psel)
					p_sel   <- input$psel
					test_sel <- input$sel_test
					FCcut <- log2(as.numeric(input$fccut))
					pvalcut <- as.numeric(input$pvalcut)
					sel_label <- "UniqueID"
					direction <- "UpDown"
					tmpdat <- GeneFilter(results_long, test_sel, p_sel, direction, pvalcut, FCcut, sel_label)
					tmpids <- tmpdat %>% as.data.frame() %>% dplyr::pull(UniqueID)
					ProteinGeneName_sel <- dplyr::filter(ProteinGeneName, UniqueID %in% tmpids)
				}

				validate(need(nrow(ProteinGeneName_sel) > 0, message = "Please input at least 1 matched gene."))
				tmpids <- ProteinGeneName_sel %>% dplyr::pull(UniqueID)

				msg_filter1 <- paste("Selected Genes:",length(ProteinGeneName_sel %>% dplyr::pull(Gene.Name) %>% unique()),sep="")
				if (length(ProteinGeneName_sel %>% dplyr::pull(UniqueID) %>% unique())  > 100)
				msg_filter1 <- paste(msg_filter1, " Too many genes, try < 100 genes", sep="")
				output$filteredgene <- renderText({msg_filter1})

				msg_filter <- paste("Selected IDs:", length(tmpids),sep="")
				if (length(ProteinGeneName_sel %>% dplyr::pull(UniqueID) %>% unique()) > length(ProteinGeneName_sel %>% dplyr::pull(Gene.Name) %>% unique()))
				msg_filter <- paste(msg_filter, " (Gene(s) match multiple IDs, use UniqueID for Gene Label)", sep="")
				output$filteredUniqueID <- renderText({msg_filter})

				numberpage = as.numeric(input$plot_ncol) * as.numeric(input$plot_nrow)
				updateSelectInput(session,'sel_page', choices= seq_len(ceiling(length(tmpids)/numberpage)))
			})

			###############
			observeEvent(input$subset , {
				req(length(working_project()) > 0)
				req(input$subset == "Geneset")
				genesetnames <- GetGeneSetNames()
				updateSelectizeInput(session, "sel_geneset", choices =  c('Type to Search' = '', genesetnames), server = TRUE)
			})

			observeEvent(input$sel_geneset, {
				req(length(working_project()) > 0)
				req(input$subset == "Geneset")
				req(input$sel_geneset!="")
				sel_geneset <- input$sel_geneset
				geneset_genenames <- GetGenesFromGeneSet(sel_geneset)
				updateTextAreaInput(session, "geneset_genes", value=paste(geneset_genenames, collapse=","))
			})

			observe({
				req(length(working_project()) > 0)
				MetaData = DataInSets[[working_project()]]$MetaData
				req(all(sapply(MetaData %>% pull('sampleid'), grepl, pattern = "^.+(_)[A-Za-z]+[0-9]+$")))

				image_ids <- MetaData %>%
				tidyr::separate(sampleid, into = c("Sample", "RowCol"), sep = "_") %>%
				pull(Sample) %>%
				unique()

				updateSelectizeInput(session,'sel_image',  choices=image_ids, selected=image_ids[1])
			})

			observeEvent(input$lcm_format, {
				req(length(working_project()) > 0)
				req(input$lcm_format)
				lcm_format = input$lcm_format
				if (lcm_format == "exp_side")
				updateSelectizeInput(session,'lcmcolpalette', choices=rownames(brewer.pal.info), selected="OrRd")
				if (lcm_format == "group_side")
				updateSelectizeInput(session,'lcmcolpalette', choices=rownames(brewer.pal.info), selected="Dark2")
				if (lcm_format == "exp_overlay")
				updateSelectizeInput(session,'lcmcolpalette', choices=rownames(brewer.pal.info), selected="OrRd")
				if (lcm_format == "group_overlay")
				updateSelectizeInput(session,'lcmcolpalette', choices=rownames(brewer.pal.info), selected="Dark2")
			})

			###############
			DataExpReactive <- reactive({
				req(length(working_project()) > 0)
				validate(need(length(DataInSets[[working_project()]]$group_order)>0,"Please select group(s)."))
				req(input$subset)
				data_long = DataInSets[[working_project()]]$data_long
				req(DataInSets[[working_project()]]$results_long)
				results_long = DataInSets[[working_project()]]$results_long
				ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
				sel_group = DataInSets[[working_project()]]$group_order
				sel_samples = DataInSets[[working_project()]]$sample_order

				genelabel=input$sel_geneid

				if (input$subset == "Select") {
					req(input$sel_gene)
					sel_gene = input$sel_gene
					validate(need(length(input$sel_gene)>0,"Please select a gene."))
					ProteinGeneName_sel <- dplyr::filter(ProteinGeneName, (UniqueID %in% sel_gene) | (Protein.ID %in% sel_gene) | (toupper(Gene.Name) %in% toupper(sel_gene)))
				}

				if (input$subset == "Upload Genes"){
					req(input$uploadlist)
					gene_list <- input$uploadlist
					gene_list <- ProcessUploadGeneList(gene_list)
					validate(need(length(gene_list)>0, message = "Please input at least 1 valid gene."))
					ProteinGeneName_sel <- dplyr::filter(ProteinGeneName, (UniqueID %in% gene_list) | (Protein.ID %in% gene_list) | (toupper(Gene.Name) %in% toupper(gene_list)))
				}

				if (input$subset == "Geneset" ) {
					req(input$sel_geneset!="")
					sel_geneset <- input$sel_geneset
					gene_list <- GetGenesFromGeneSet(sel_geneset)
					validate(need(length(gene_list)>0, message = "Please input at least 1 valid gene."))
					ProteinGeneName_sel <- dplyr::filter(ProteinGeneName, (UniqueID %in% gene_list) | (Protein.ID %in% gene_list) | (toupper(Gene.Name) %in% toupper(gene_list)))
				}

				if (input$subset == "Browsing") {
					req(input$sel_test)
					req(input$psel)
					p_sel   <- input$psel
					test_sel <- input$sel_test
					FCcut <- log2(as.numeric(input$fccut))
					pvalcut <- as.numeric(input$pvalcut)
					sel_label <- "UniqueID"
					direction <- "UpDown"
					tmpdat <- GeneFilter(results_long, test_sel, p_sel, direction, pvalcut, FCcut, sel_label)
					tmpids <- tmpdat %>% as.data.frame() %>% dplyr::pull(UniqueID)
					ProteinGeneName_sel <- dplyr::filter(ProteinGeneName, UniqueID %in% tmpids)
				}

				validate(need(nrow(ProteinGeneName_sel) > 0, message = "Please input at least 1 matched gene."))
				tmpids <- ProteinGeneName_sel %>% dplyr::pull(UniqueID)

				if (!input$tabset %in% c("Expression Correlation","Rank Abundance Curve")) {
					numberpage = as.numeric(input$plot_ncol) * as.numeric(input$plot_nrow)
					updateSelectInput(session,'sel_page', choices= seq_len(ceiling(length(tmpids)/numberpage)))

					req(input$sel_page)
					sel_page = as.numeric(input$sel_page) - 1
					startslice = sel_page * numberpage  + 1
					endslice = startslice + numberpage - 1
					if (endslice > length(tmpids))
					endslice <- length(tmpids)
					tmpids <- tmpids[startslice: endslice]
				}

				data_long_tmp = dplyr::filter(data_long, UniqueID %in% tmpids, group %in% sel_group, sampleid %in% sel_samples) %>%
				dplyr::filter(!is.na(expr)) %>%
				dplyr::filter(expr != 0) %>% ## considering change to expr + 1, if the values are not distribute around 1 (ratio)
				as.data.frame()

				data_long_tmp$UniqueID <- factor(data_long_tmp$UniqueID, levels = tmpids)
				data_long_tmp <- data_long_tmp %>% arrange(UniqueID)

				ylogbase = as.numeric(input$ylogbase)
				if (ylogbase == 2){
					data_long_tmp <- data_long_tmp %>% dplyr::mutate_at(vars(expr), log2)
				}
				if (ylogbase == 10){
					data_long_tmp <- data_long_tmp %>% dplyr::mutate_at(vars(expr), log10)
				}

				data_long_tmp$labelgeneid = data_long_tmp[,match(genelabel,colnames(data_long_tmp))]
				data_long_tmp$labelgeneid <- factor(data_long_tmp$labelgeneid, levels = unique(data_long_tmp$labelgeneid))
				data_long_tmp$group = factor(data_long_tmp$group, levels = sel_group)

				result_long_tmp = dplyr::filter(results_long, UniqueID %in% tmpids) %>% as.data.frame()
				result_long_tmp$labelgeneid = result_long_tmp[,match(genelabel,colnames(result_long_tmp))]
				result_long_tmp$labelgeneid <- factor(result_long_tmp$labelgeneid, levels = unique(result_long_tmp$labelgeneid))

				return(list("data_long_tmp" = data_long_tmp,"result_long_tmp"= result_long_tmp, "tmpids"=tmpids))
			})

			############### boxplot
			boxplot_out <-	eventReactive(input$plot_exp, {
				withProgress(message = 'Making Plot. It may take a while...', value = 0, {
					req(length(working_project()) > 0)
					sel_group = DataInSets[[working_project()]]$group_order
					tests_order = DataInSets[[working_project()]]$tests_order
					MetaData = DataInSets[[working_project()]]$MetaData

					results_long_tmp <- DataExpReactive()$result_long_tmp %>%
					dplyr::filter(test %in% tests_order)

					data_long_tmp <- DataExpReactive()$data_long_tmp %>%
					dplyr::select(-any_of(c("id","UniqueID", "Gene.Name", "Protein.ID")))

					colorby = sym(input$colorby)
					plotx = sym(input$plotx)
					Val_colorby = input$colorby
					Val_plotx = input$plotx

					xlogbase = as.numeric(input$xlogbase)
					ylogbase = as.numeric(input$ylogbase)

					ncol = as.numeric(input$plot_ncol)
					nrow = as.numeric(input$plot_nrow)
					numberpage = ncol * nrow

					genelabel = input$sel_geneid
					barcol = input$barcol

					if (Val_colorby != Val_plotx) {
						if (Val_colorby!="group" ) {
							data_long_tmp <- data_long_tmp %>% dplyr::left_join((MetaData %>% dplyr::select(sampleid, !!colorby)), by = "sampleid")
						}

						if (Val_plotx!="group" ) {
							data_long_tmp <- data_long_tmp %>% dplyr::left_join((MetaData %>% dplyr::select(sampleid, !!plotx)), by = "sampleid")
						}
					} else {
						if (Val_colorby!="group" ) {
							data_long_tmp <- data_long_tmp %>% dplyr::left_join((MetaData %>% dplyr::select(sampleid, !!colorby)), by = "sampleid")
						}
					}

					showPvalues = input$PvalueBar
					value = sym("expr")
					pval_sel = input$pval_sel
					psel = sym(pval_sel)

					results_long_tmp[,sapply(results_long_tmp,is.numeric)] <- signif(results_long_tmp[,sapply(results_long_tmp,is.numeric)],3)
					results_long_tmp <- results_long_tmp %>%
					dplyr::select(-any_of(c("id","UniqueID", "Gene.Name", "Protein.ID")))

					if (Val_plotx == Val_colorby) {
						df.summary <- data_long_tmp %>%
						group_by(!!plotx, labelgeneid) %>%
						dplyr::summarise(sd = sd(!!value), expr = mean(!!value),  n = n(), se = sd / sqrt(n), .groups = "drop") %>%
						dplyr::mutate(ypos = sum(expr, sd, na.rm=TRUE)) %>%
						dplyr::mutate(ypos3sd = sum(expr,3*sd, na.rm=TRUE)) %>%
						dplyr::select(-n)
					} else {
						df.summary <- data_long_tmp %>%
						group_by(!!plotx, !!colorby, labelgeneid) %>%
						dplyr::summarise(sd = sd(!!value), expr = mean(!!value),  n = n(), se = sd / sqrt(n), .groups = "drop") %>%
						dplyr::mutate(ypos = sum(expr, sd, na.rm=TRUE))  %>%
						dplyr::mutate(ypos3sd = sum(expr,3*sd, na.rm=TRUE)) %>%
						dplyr::select(-n)
					}

					annotation_df <- data.frame()
					if((showPvalues == "YES") & (Val_plotx == "group") & (Val_colorby == "group")){
						pv_thresh = as.numeric(input$pv_threshold)

						y_position <- df.summary %>%
						rowwise() %>%
						#	dplyr::mutate(ypos = expr + sd) %>% #sd is NA for n =1
						dplyr::mutate(ypos = sum(expr, sd, na.rm=TRUE))  %>%
						dplyr::select(-c(sd,se,expr)) %>%
						dplyr::group_by(labelgeneid) %>%
						dplyr::mutate(max_ypos = max(ypos, na.rm=TRUE)) %>%
						dplyr::mutate(max_ypos_3sd = max(ypos3sd, na.rm=TRUE)) %>%
						dplyr::ungroup() %>%
						as.data.frame()

						annotation_df = results_long_tmp %>%
						dplyr::filter(!!psel < pv_thresh) %>%
						as.data.frame()

						if (nrow(annotation_df) > 0)  {
							annotation_df <- annotation_df %>%
							tidyr::separate(test, c('start', 'end'),  sep = "vs") %>%
							dplyr::mutate_at(c("start","end"), trimws, which = "both") %>%
							dplyr::mutate(Label = case_when(!!psel <= 0.001 ~ "***",
								!!psel > 0.001 & !!psel <= 0.01 ~ "**",
								!!psel > 0.01 & !!psel <= 0.05 ~ "*",
								!!psel > 0.05 & !!psel <= 0.10 ~ ".",
							TRUE ~"")) %>%
							dplyr::mutate(group = start) %>%
							dplyr::inner_join(y_position, by=c("start"="group", "labelgeneid"="labelgeneid")) %>%
							dplyr::inner_join(y_position %>% dplyr::select(-c(max_ypos,max_ypos_3sd)), by=c("end"="group", "labelgeneid"="labelgeneid")) %>%
							dplyr::group_by(labelgeneid) %>%
							dplyr::mutate(count = n()) %>%
							dplyr::mutate(y= 1.02 * pmax(ypos.x, ypos.y)) %>%
							dplyr::arrange(y, .by_group = TRUE) %>%
							#dplyr::mutate_if(is.numeric, ~replace(., is.na(.), 0))  %>%
							dplyr::mutate(y1 = seq(min(y), max(y), length.out = count[1])) %>%
							dplyr::mutate(y2 = seq(max_ypos[1], 1.15 * max_ypos[1], length.out = count[1])) %>%
							dplyr::mutate(y3 = seq(max_ypos_3sd[1], 1.15 * max_ypos_3sd[1], length.out = count[1])) %>%
							dplyr::select(-c(ypos.x, ypos.y,ypos3sd.x, ypos3sd.y, count, logFC)) %>%
							dplyr::ungroup() %>%
							dplyr::filter(start %in% sel_group & end %in% sel_group)
						}
					}

					data_long_tmp <- data_long_tmp %>%
					dplyr::mutate_at(vars(!!plotx), as.factor)

					pd <- position_dodge(0.8)

					if (input$SeparateOnePlot == "Separate") {
						p <- ggplot(data_long_tmp, aes(x=!!plotx, y=!!value, alpha = !!colorby,  color=!!colorby))
						p <- p + facet_wrap(~labelgeneid, scales = "free", nrow = nrow, ncol = ncol)
					} else {
						p <- ggplot(data_long_tmp, aes(x=!!plotx, y=!!value,  alpha = labelgeneid, color = labelgeneid))
					}

					if (input$plotformat == "boxplot") {
						p <- p + geom_boxplot(outlier.colour=NA)
					}

					if (input$plotformat == "barplot") {
						if (input$SeparateOnePlot == "Separate") {
							p <- p + geom_col(data = df.summary, aes(x = !!plotx, y=!!value, fill = !!colorby), position = pd, width = 0.7)
						} else {
							p <- p + geom_col(data = df.summary, aes(x = !!plotx, y=!!value, fill = labelgeneid), position = pd, width = 0.7)
						}
					}

					if (input$plotformat == "violin") {
						p <- p + geom_violin()
					}

					if (input$plotformat == "line") {
						if (input$SeparateOnePlot == "Separate") {
							if (input$plotx == input$colorby) {
								p <- ggplot(data_long_tmp, aes(x=!!plotx, y=!!value))
							} else {
								p <- ggplot(data_long_tmp, aes(x=!!plotx, y=!!value, fill=!!colorby, color=!!colorby))
							}
							p <- p + geom_line(data = df.summary, aes(x = !!plotx, y=!!value), group=1) +
							facet_wrap(~labelgeneid, scales = "free", nrow = nrow, ncol = ncol)

							if (input$ShowErrorBar == "SD") {
								p <- p + geom_errorbar(aes(ymin = !!value-sd, ymax = !!value+sd), data = df.summary, width = 0.2, position = pd)
							}
							if (input$ShowErrorBar == "SEM") {
								p <- p + geom_errorbar(aes(ymin = !!value-se, ymax = !!value+se), data = df.summary, width = 0.2, position = pd)
							}
						} else {
							p <- ggplot(data_long_tmp, aes(x = !!plotx, y= !!value)) +
							geom_line(data = df.summary, aes(x = !!plotx, y = !!value, color = labelgeneid, group = labelgeneid))

							if (input$ShowErrorBar == "SD") {
								p <- p + geom_errorbar(aes(ymin = !!value-sd, ymax = !!value+sd,  color = labelgeneid,), data = df.summary, width = 0.2, position = pd )
							}
							if (input$ShowErrorBar == "SEM") {
								p <- p + geom_errorbar(aes(ymin = !!value-se, ymax = !!value+se,  color = labelgeneid), data = df.summary, width = 0.2, position = pd)
							}
						}
					}

					#######################
					if (input$plotformat != "boxplot" & input$plotformat != "line") {
						if (input$ShowErrorBar == "SD") {
							p <- p + geom_errorbar(aes(ymin = !!value-sd, ymax = !!value+sd), data = df.summary, width = 0.2, position = pd)
						}
						if (input$ShowErrorBar == "SEM") {
							p <- p + geom_errorbar(aes(ymin = !!value-se, ymax = !!value+se), data = df.summary, width = 0.2, position = pd)
						}
					}

					if((nrow(annotation_df) > 0) & (showPvalues == "YES") & (Val_plotx == "group") &  (Val_colorby == "group")){
						plabelpos = sym(input$plabelpos)
						if (input$SeparateOnePlot == "Separate") {
							p <- p + geom_signif(
								data = annotation_df,
								aes(xmin = start, xmax = end, annotations = Label, color = !!colorby,  y_position = !!plabelpos),
								textsize = 6,
								tip_length = 0,
								vjust = 0.8,
								margin_top = 0.05,
								step_increase = 0.5,
								manual = TRUE,
								inherit.aes = TRUE,
								show.legend = FALSE
							)
						} else {
							p <- p + geom_signif(
								data = annotation_df,
								aes(xmin = start, xmax = end, annotations = Label, color =labelgeneid,  y_position = !!plabelpos),
								textsize = 6,
								tip_length = 0,
								vjust = 0.8,
								margin_top = 0.05,
								step_increase = 0.5,
								manual = TRUE,
								inherit.aes = TRUE,
								show.legend = FALSE
							)
						}
					}

					if (input$IndividualPoint == "YES") {
						if (input$SeparateOnePlot == "Separate") {
							p <- p + geom_jitter(aes(color=!!colorby), position = position_jitterdodge())
						} else {
							p <- p + geom_jitter(aes(color=labelgeneid), position = position_jitterdodge())
						}
					}

					if (xlogbase != 1){
						p <- p +  scale_x_continuous(trans=scales::pseudo_log_trans(base = xlogbase))
					}

					if (input$SeparateOnePlot == "Separate") {
						if (input$colpalette == "Single") {
							colorpal = rep(barcol, length(sel_group))
						} else {
							colorpal <- UserColorPlalette(colpalette = input$colpalette, items = sel_group)
						}
						user_alpha = rep(1, length(sel_group))
					} else {
						if (input$colpalette == "Single") {
							colorpal = rep(barcol, length(unique(data_long_tmp$labelgeneid)))
						} else {
							colorpal <- UserColorPlalette(colpalette = input$colpalette, items = unique(data_long_tmp$labelgeneid))
						}
						user_alpha = rep(1, length(unique(data_long_tmp$labelgeneid)))
					}

					p <- p +
					scale_alpha_manual(values = user_alpha) +
					scale_fill_manual(values = colorpal) +
					scale_colour_manual(values = colorpal) +
					theme_bw(base_size = 16) +
					ylab(input$Ylab) +
					xlab(input$Xlab)

					if (input$colpalette == "Single") {
						p <- p +  theme (legend.position="none")
					}

					if (input$plot_Y_scale=="Manual") {
						p <- p + ylim(input$plot_Ymin, input$plot_Ymax)
					}

					if (input$prismstyle=="YES") {
						p <- p + scale_color_prism(input$prismcolpalette) +
						scale_fill_prism(input$prismfillpalette) +
						guides(y = "prism_offset_minor") +
						theme_prism(base_size = 16)
					}

					p <- p +
					theme (plot.margin = unit(c(1,1,1,1), "cm"),
						text = element_text(size=input$axisfontsize),
						axis.text.x = element_text(angle = input$Xangle, hjust=0.5, vjust=0.5),
						strip.text.x = element_text(size=input$titlefontsize)
					)
					return(p)
				})
			})

			output$plot.exp <- renderUI({
				if (input$SeparateOnePlot == "Separate")  {
					nrow = isolate(as.numeric(input$plot_nrow))
					height = 300 * nrow
				} else {
					height=800
				}

				if (input$interactive == "static"){
					plotOutput(ns("boxplot"), height = height)
				}	else {
					plotly::plotlyOutput(ns("boxplot_Interactive"),	height=height)
				}
			})

			output$boxplot_Interactive <- renderPlotly({
				p <- boxplot_out()
				pl <- ggplotly(p) %>% layout(dragmode = "pan")
			})

			output$boxplot <- renderPlot({
				boxplot_out() %>% print
			})

			observeEvent(input$saveboxplot, {
				saved.num <- length(saved_plots$boxplot) + 1
				saved_plots$boxplot[[saved.num]] <- boxplot_out()
			})

			output$dat_dotplot <- DT::renderDT(server=FALSE, {
				req(length(working_project()) > 0)
				data_long_tmp <- DataExpReactive()$data_long_tmp
				data_long_tmp <- data_long_tmp %>%dplyr::select(-labelgeneid)
				data_long_tmp[,sapply(data_long_tmp,is.numeric)] <- signif(data_long_tmp[,sapply(data_long_tmp,is.numeric)],3)
				DT::datatable(data_long_tmp,  extensions = 'Buttons',  options = list(
					dom = 'lBfrtip', pageLength = 15,
					buttons = list(
						list(extend = "csv", text = "Download Page", filename = "Page_results",
							exportOptions = list(modifier = list(page = "current"))
						), list(extend = "csv", text = "Download All", filename = "All_Results",
						exportOptions = list(modifier = list(page = "all"))))
						)
					)
				}
			)

			output$res_dotplot <- DT::renderDT(server=FALSE,{
				req(length(working_project()) > 0)
				result_long_tmp <- DataExpReactive()$result_long_tmp
				result_long_tmp[,sapply(result_long_tmp,is.numeric)] <- signif(result_long_tmp[,sapply(result_long_tmp,is.numeric)],3)
				DT::datatable(result_long_tmp,  extensions = 'Buttons',  options = list(
					dom = 'lBfrtip', pageLength = 15,
					buttons = list(
						list(extend = "csv", text = "Download Page", filename = "Page_results",
							exportOptions = list(modifier = list(page = "current"))
						),
						list(extend = "csv", text = "Download All", filename = "All_Results",
							exportOptions = list(modifier = list(page = "all"))
						)
					)
				))
			})

			############### S curve
			Scurve_out <- eventReactive(input$plot_SCurve, {
				req(length(working_project()) > 0)
				data_results <- DataInSets[[working_project()]]$data_results
				gene_list <- DataExpReactive()$tmpids
				validate(need(length(gene_list)>0,"Please select a gene or input gene."))
				validate(need(("Intensity" %in% colnames(data_results)),"Need Intensity in data_result table."))

				scurve.data <- data_results %>%
				dplyr::select(any_of(c("UniqueID","Gene.Name","Protein.ID","Intensity"))) %>%
				dplyr::filter((!is.na(Intensity)) & Intensity > 0) %>%
				dplyr::arrange(desc(Intensity)) %>%
				dplyr::mutate(RANK = row_number())

				genelabel=input$sel_geneid
				scurve.data$labelgeneid = scurve.data[,match(genelabel,colnames(scurve.data))]
				scurve.data$labelgeneid <- factor(scurve.data$labelgeneid, levels = unique(scurve.data$labelgeneid))

				scurve.label <- scurve.data %>%
				dplyr::filter((UniqueID %in% gene_list) | (Gene.Name %in% gene_list) | (Protein.ID %in% gene_list))

				if (nrow(scurve.label) > 0) {
					p <- ggplot(scurve.data, aes(x = RANK, y = Intensity)) +
					geom_point(colour = "blue")  +
					geom_point(data=scurve.label, aes(x = RANK, y = Intensity, colour='red')) +
					scale_y_continuous(trans='log10')  +
					geom_label_repel(
						data = scurve.label,aes(x = RANK, y = Intensity,  label = labelgeneid),
						fontface = 'bold', color = 'red', size = input$scurve_labelfontsize,
						box.padding = unit(0.35, "lines"),
						point.padding = unit(0.5, "lines"),
						segment.color = 'red'
					) +
					theme_bw(base_size = 14) +
					ylab(input$scurveYlab) +
					xlab(input$scurveXlab) +
					theme (plot.margin = unit(c(1,1,1,1), "cm"),
						text = element_text(size=input$scurve_axisfontsize),
						axis.text.x = element_text(angle = input$scurveXangle),
					legend.position="none")

					p <- ggExtra::ggMarginal(p, type = input$scurveright, margins = "y", colour="blue", fill = "gray")
				}
				return(p)
			})

			output$SCurve <- renderPlot({
				Scurve_out()
			})

			observeEvent(input$AbundanceCurve, {
				saved_plots$AbundanceCurve <- Scurve_out()
			})

			############### Laser capture microdissection
			LCM_out <- eventReactive(input$plot_LCM, {
				req(length(working_project()) > 0)
				sel_group = DataInSets[[working_project()]]$group_order
				tests_order = DataInSets[[working_project()]]$tests_order
				MetaData = DataInSets[[working_project()]]$MetaData
				validate(need(all(sapply(MetaData %>% pull('sampleid'), grepl, pattern = "^.+(_)[A-Za-z]+[0-9]+$")),"Only work for LCM format"))

				data_long_tmp <- DataExpReactive()$data_long_tmp %>%
				dplyr::select(-any_of(c("id","UniqueID", "Gene.Name", "Protein.ID")))
				tmpids <- DataExpReactive()$tmpid

				ncol = as.numeric(input$plot_ncol)
				nrow = as.numeric(input$plot_nrow)
				numberpage = ncol * nrow

				library(png)
				library(grid)
				library(gridExtra)

				LCM_alpha = as.numeric(input$LCM_alpha)
				sel_image = input$sel_image

				pngfile = paste("data/", sel_image, ".png", sep="")
				if (file.exists(pngfile)) {
					im <- readPNG(pngfile)
				}

				data_long_tmp <-  data_long_tmp %>%
				tidyr::separate(sampleid, into = c("Sample", "RowCol"), sep = "_")%>%
				tidyr::separate(RowCol, into = c("Row", "Col"), "(?<=[A-Za-z])(?=[0-9])") %>%
				dplyr::mutate(Row = forcats::fct_rev(Row)) %>%
				dplyr::group_by(labelgeneid) %>%
				dplyr::mutate(ratio = log2(expr/mean(expr))) %>%
				dplyr::ungroup()

				convert2ratio <- input$convert2ratio

				if (length(tmpids) > 1) {
					p1 <- ggplot(data_long_tmp, aes(Col, y= Row, fill=!!sym(convert2ratio))) +
					geom_tile()
					if (convert2ratio == "ratio")
					p1 <- p1 + scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", limits=c(input$LowLimit,input$HighLimit))
					else
					p1 <- p1 + scale_fill_distiller(palette = input$lcmcolpalette, direction = 1)

					p1 <- p1 +
					theme_bw(base_size = 14) +
					ylab(input$LCMYlab) +
					xlab(input$LCMXlab) +
					theme(plot.margin = unit(c(1,1,1,1), "cm"), text = element_text(size=input$LCM_axisfontsize), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

					p <- p1 + facet_wrap(~labelgeneid, scales = "free", nrow = nrow, ncol = ncol)
				}	else {
					p1 <- ggplot(data_long_tmp, aes(Col, y= Row, fill=!!sym(convert2ratio))) +
					geom_tile()

					if (convert2ratio == "ratio")
					p1 <- p1 + scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", limits=c(input$LowLimit,input$HighLimit))
					else
					p1 <- p1 + scale_fill_distiller(palette = input$lcmcolpalette, direction = 1)

					p1 <- p1 +
					coord_fixed(ratio = nrow(im)/ncol(im)) +
					theme_bw(base_size = 14) +
					ylab(input$LCMYlab) +
					xlab(input$LCMXlab) +
					theme(plot.margin = unit(c(1,1,1,1), "cm"), text = element_text(size=input$LCM_axisfontsize), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

					colorpal <- UserColorPlalette(colpalette = input$lcmcolpalette, items = unique(data_long_tmp$group))
					p2 <- ggplot(data_long_tmp, aes(Col, y= Row, fill=group)) +
					geom_tile() +
					scale_fill_manual(values =  colorpal) +
					coord_fixed(ratio = nrow(im)/ncol(im)) +
					theme_bw(base_size = 14) +
					ylab(input$LCMYlab) +
					xlab(input$LCMXlab) +
					theme(plot.margin = unit(c(1,1,1,1), "cm"), text = element_text(size=input$LCM_axisfontsize), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

					im2 <- matrix(rgb(im[,,1],im[,,2],im[,,3],im[,,4] * 1), nrow=dim(im)[1])

					lcm_format = input$lcm_format
					if (lcm_format == "exp_side")
					p =  grid.arrange(grid::rasterGrob(im2), p1, ncol=2)
					if (lcm_format == "group_side")
					p =  grid.arrange(grid::rasterGrob(im2), p2, ncol=2)
					if (lcm_format == "exp_overlay"){
						im2 <- matrix(rgb(im[,,1],im[,,2],im[,,3],im[,,4] * LCM_alpha), nrow=dim(im)[1])
						p <- p1 +  annotation_custom(rasterGrob(im2,	width = unit(1,"npc"), height = unit(1,"npc")),	-Inf, Inf, -Inf, Inf)
					}
					if (lcm_format == "group_overlay") {
						im2 <- matrix(rgb(im[,,1],im[,,2],im[,,3],im[,,4] * LCM_alpha), nrow=dim(im)[1])
						p <- p2 + annotation_custom(rasterGrob(im2,	width = unit(1,"npc"), height = unit(1,"npc")),	-Inf, Inf, -Inf, Inf)
					}

				}
				return(p)
			})

			output$LCM <- renderPlot({
				LCM_out()
			})

			observeEvent(input$LCM, {
				saved_plots$LCM <- LCM_out()
			})

			############### Expr Corr
			ExprCorr_out <- eventReactive(input$plot_ExprCorr, {
				req(length(working_project()) > 0)
				sel_group = DataInSets[[working_project()]]$group_order
				tests_order = DataInSets[[working_project()]]$tests_order
				MetaData = DataInSets[[working_project()]]$MetaData

				data_long_tmp <- DataExpReactive()$data_long_tmp %>%
				dplyr::select(-any_of(c("id","UniqueID", "Gene.Name", "Protein.ID")))
				tmpids <- DataExpReactive()$tmpid
				validate(need(length(tmpids) >= 2,"select > 2 Ids"))

				library(corrplot)

				if (length(tmpids) > 2) {
					cor_matrix <-  data_long_tmp  %>% as.data.frame() %>%
					dplyr::select(-c(group)) %>%
					pivot_wider(names_from = labelgeneid, values_from = expr) %>%
					dplyr::select(-c(sampleid)) %>% as.matrix() %>%
					cor(.,use = "p")
					testRes = cor.mtest(cor_matrix, conf.level = 0.95)

					if(input$diag=="YES")
					diag = TRUE
					else
					diag = FALSE
					
					p <- corrplot(cor_matrix,
						method=input$method,
						order=input$order,
						tl.pos="lt",
						type=input$type ,
						tl.col="black",
						tl.cex=input$tlcex,
						tl.srt=45, addCoef.col="black", 
						addCoefasPercent = FALSE,
						hclust.method=input$hclustmethod,
						plotCI="n",
						p.mat = testRes$p,
						sig.level=input$siglevel,
						col = COL2(input$corrcol),
						diag = diag,
					insig = "blank")
				}	else if (length(tmpids)  == 2){
					wide_temp <-  data_long_tmp  %>% as.data.frame() %>%
					dplyr::select(-c(group)) %>%
					pivot_wider(names_from = labelgeneid, values_from = expr) %>%
					dplyr::select(-c(sampleid)) #%>% as.matrix()

					c.res <- cor(wide_temp[,1], wide_temp[,2], use = "complete.obs")
					cor_string <- paste("Pearson corr:", format(c.res[1,1], digits=3), sep="")
					p <- ggplot(data = wide_temp,	aes(x = eval(parse(text = colnames(wide_temp)[1])),	y =eval(parse(text = colnames(wide_temp)[2])))) +
					geom_smooth(method=lm) +
					geom_point(size = 2) +
					labs(title=cor_string) +
					xlab(colnames(wide_temp)[1]) + ylab(colnames(wide_temp)[2])
				}
				return(p)
			})

			output$ExprCorr<- renderPlot({
				ExprCorr_out()
			})

			observeEvent(input$ExprCorr, {
				saved_plots$ExprCorr <- ExprCorr_out()
			})
		}
	)
}

