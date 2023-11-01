###########################################################################################################
## Multiomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: app.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 05/23/2023
##@version 3.0
###########################################################################################################
source("global.R",local = TRUE)$value

ui <- fluidPage(
	theme = shinytheme("cerulean"),
	titlePanel(
		fluidRow(
			column(4, img(height =75 , src = "")),
			column(8,  h2(strong("Multiomics Visualization"), align = 'left'))
		),
	windowTitle = "Multiomics Visualization"),
	includeCSS("menuhexagonal1.css"),
	navbarPage(title = "", id="menu", selected = "Dataset",
		#tabPanel(textOutput('project')),
		##########################################################################################################
		## Select Dataset
		##########################################################################################################
		tabPanel("Dataset",
			fluidRow(
				column(12,
					tabsetPanel(id="Tables",
						tabPanel(title="Select Dataset",
							fluidRow(
								column(3,
									#radioButtons("select_dataset", label="Select data set", choices=c("Saved Projects in CSV file", "Saved Projects in Database", "Upload RData File", "Upload Data Files (csv)", "Public Data(DiseaseLand)", "Long Format", "Wide Format"), inline = F, selected="Saved Projects in CSV file"),
									radioButtons("select_dataset", label="Select data set", choices=c("Saved Projects in CSV file", "Saved Projects in Database", "Upload RData File", "Upload Data Files (csv)", "Public Data(DiseaseLand)"), inline = F, selected="Saved Projects in CSV file"),
									uiOutput('loadedprojects')
								),
								column(5,
									conditionalPanel("input.select_dataset=='Saved Projects in Database' || input.select_dataset=='Saved Projects in CSV file' || input.select_dataset=='Public Data(DiseaseLand)'",
										selectInput("sel_project", label="Available Dataset",	choices=NULL, width = "80%"),
										uiOutput("loaddata")
									),
									conditionalPanel("input.select_dataset=='Upload RData File'",
										fileInput("file1", "Choose data file"),
										fileInput("file2", "(Optional) Choose network file"),
										uiOutput('ui.action')#,
										#uiOutput("loaddata")
									),
									conditionalPanel("input.select_dataset=='Upload Data Files (csv)'",
										uiOutput('upload.files.ui')#,
										#uiOutput("loaddata")
									)#,
									#	conditionalPanel(condition = "input.select_dataset == 'Wide Format'",
									#	fileInput('datafile_meta', '1st Step: Input meta file', accept = c('.csv', '.txt')),
									#	fileInput('datafile_wide', '2nd Step: Input data file', accept = c('.csv', '.txt')
									#		)
									#),
									#conditionalPanel(condition = "input.select_dataset == 'Long Format'",
									#	fileInput('datafile_long', 'Input data file', accept = c('.csv', '.txt')),
									#	textOutput('selcolumns'),
									#	uiOutput('uniqueID'),
									#	uiOutput('group'),
									#	uiOutput('conc'),
									#	uiOutput('response'),
									#	uiOutput('actionbutton')
									#)
								),
								column(4,
									htmlOutput("summary"),
									tableOutput('group_table'), tags$br(),
									uiOutput('exp_unit'),
									uiOutput('ShortName'),
									uiOutput('comp_info'),
									uiOutput("removedata")
								)
							)
						),
						tabPanel(title="Sample Table", value="Sample Table",
							DT::dataTableOutput('sample'),
							actionButton("sample", "Save data for output as one Excel file")
						),
						tabPanel(title="Result Table", value="Result Table",
							DT::dataTableOutput('results'),
							shiny::downloadButton(outputId = "download_result_button",  label = "Download All as CSV file"),
							actionButton("results", "Save data for output as one Excel file")
						),
						tabPanel(title="Data Table", value="Data Table",
							DT::dataTableOutput('data_wide'),
							shiny::downloadButton(outputId = "download_data_button",  label = "Download All as CSV file"),
							actionButton("data_wide", "Save data for output as one Excel file")
						),
						tabPanel(title="Protein Gene Names", value="Protein Gene Names",
							DT::dataTableOutput('ProteinGeneName'),
							shiny::downloadButton(outputId = "download_ProteinGeneName_button",  label = "Download All as CSV file"),
							actionButton("ProteinGeneName", "Save data for output as one Excel file")
						),
						tabPanel(title="Help", value="Help", htmlOutput('help_input')
						)
					)
				)
			)
		),
		##########################################################################################################
		##########################################################################################################
		## footer
		##########################################################################################################
		footer= HTML(footer_text)
	)
)

server <- function(input, output, session) {
	saved_setting <- reactiveValues(value = c())
	source("inputdata.R",local = TRUE)
	source("help.R",local = TRUE)

	modulelist <- c("QC Plots", "DEGs", "Heatmap", "Expression Plot", "Geneset Enrichment", "Pattern",
	"Correlation Network", "Venn diagram", "PCSF", "WGCNA", "Dose Response", "Dromics", "Monotonic Trend", "Merge Data", "Groups and Samples")
	moduleFilelist <- c("qcplotmodule.R", "degmodule.R", "heatmapmodule.R", "expressionmodule.R",
		"genesetmodule.R", "patternmodule.R", "networkmodule.R", "vennmodule.R", "pcsfmodule.R", "wgcnamodule.R", "drcmodule.R",
	"dromicsmodule.R", "monotonicmodule.R", "mergedatamodule.R","groupsamplemodule.R")

	observe({
		req(input$select_dataset=='Upload Data Files (csv)')
		library(stringi)
		library(biomaRt)

		#Upload CVS Data Files
		output$upload.files.ui <- renderUI({
			tagList(tags$div(
			tags$p("Prepare your own data files in Excel, save them as csv files and upload here. The system will automatically process the files and create the R data files. You need sample metadata file, expression data file, and comparison data file. The system can create gene/protein annotation based on the IDs from data files, or you can upload your own Gene/Protein Name file.")),
				tags$a(href="RNA_Seq_Demo.zip", "Download RNA-Seq example csv files (200 genes from mouse microglia dataset)"),
				tags$br(),
				tags$a(href="Proteomics_Demo.zip", "Download Proteomics example csv files (200 proteins from AD PD dataset)"),
				tags$hr(),
				textInput("F_project_name", label="Project Name", value=""),
				radioButtons("Fspecies",label="Select Species", choices=c("human","mouse", "rat"), inline = T, selected="human"),
				tags$hr(),
				tags$p("Sample MetaData must have sampleid and group columns, with additional columns optional. The sample names in sampleid column must match the expression data file."),
				fileInput("F_sample", "Sample MetaData File"),
				tags$hr(),
				tags$p("Expression data should be matrix of expression values with genes/proteins as rows, and samples as columns.  The unique IDs for genes/proteins are in the first column. We recommend using log of normalized expression values (e.g. log2(TPM+1). Upload csv file, can be compressed as .gz or .zip file."),
				fileInput("F_exp", "Expression Data File"),
				tags$hr(),
				tags$p("Comparison data should have five columns, UniqueID, test, Adj.P.Value, P.Value and logFC. The comparison names are listed in test column. Upload csv file, can be compressed as .gz or .zip file."),
				fileInput("F_comp", "Comparison Data File"),
				tags$hr(),
				checkboxInput("F_annot_auto", "Create Gene/Protein Name File automatically (or uncheck to upload your own file)", TRUE, width="90%"),
				#radioButtons("F_annot_auto", label="Create Gene/Protein Names automatically:", inline = TRUE, choices = c("Yes","No"), selected = "Yes"),
				conditionalPanel(condition="input.F_annot_auto==1",
					radioButtons("F_ID_type",label="Unique ID Type in the Data Files", choices=c("Ensembl Gene ID", "Gene Symbol", "NCBI GeneID","UniProtKB Protein ID", "UniProt Protein Name"), inline = T, selected="Ensembl Gene ID"),
					checkboxInput("F_ID_info", "Show ID type examples", FALSE, width="90%"),
					conditionalPanel(condition="input.F_ID_info==1",
						tags$p("The system will extract gene/protein names based on the unique IDs from your data using BioMart or UniProt database. The unique IDs from the expression and comparison data can be one of the following formats:"),
						tags$ol(
							tags$li("Ensembl Gene ID (e.g. ENSG00000118260, or with version number ENSG00000118260.14)"),
							tags$li("Gene Symbol (e.g. CREB1)"),
							tags$li("NCBI GeneID (e.g. 1385)"),
							tags$li("UniProtKB Protein ID (e.g. P16220, C9J4L5, P16220-2)"),
							tags$li("UniProt Protein Name (e.g. CREB1_HUMAN, C9J276_HUMAN)")
						)
					),
					checkboxInput("F_fillName", "Fill in uniqueID when Gene.Name not found", TRUE),
					checkboxInput("F_description", "Add gene/protein description", FALSE)
				),
				# radioButtons("F_description",label="Add gene/protein description?", choices=c("Yes", "No"), inline = T, selected="No"),
				conditionalPanel(condition="input.F_annot_auto==0",
					tags$p("The Gene/Protein Name csv file must have four columns: id (sequential numbers), UniqueID (match with the IDs in the expression and comparison data file), Gene.Name (official gene symbols), Protein.ID (UniProt protein IDs, or enter empty values for RNA-Seq data). Additional columns (e.g. gene biotype) can be added."),
					fileInput("F_annot", "Gene/Protein Name File")
				),
				radioButtons("savetoserver", label="Save to Server", choices=c("YES","NO"), inline = T, selected="NO"),
				radioButtons("folder_name", label="Folder Name", choices=c("unlisted", "data"), inline = T, selected="unlisted"),
				tags$hr(),

				actionButton("uploadData", "Submit Data")
			)
		})
	})


	### Custom Data
	output$ui.action <- renderUI({
		req(input$file1)
		ProjectID = str_replace(input$file1$name, regex(".RData", ignore_case = TRUE), "")
		tagList(
			textInput("project_name", label="Rename Project", value=ProjectID ),
			radioButtons("species",label="Select species", choices=c("human","mouse", "rat"), inline = T, selected="human"),
			radioButtons("savetoserver", label="Save to Server", choices=c("YES","NO"), inline = T, selected="NO"),
			radioButtons("folder_name", label="Folder Name", choices=c("unlisted", "data"), inline = T, selected="unlisted"),
			actionButton("customData", "Submit Data")
			#uiOutput("loaddata")
		)
	})


	observe({
		req(length(reactiveValuesToList(saved_plots)) !=0 | length(reactiveValuesToList(saved_table)) !=0 )
		if (!("output" %in% saved_setting$value)) {
			source("outputmodule.R",  local = TRUE)
			appendTab(session=session, inputId = "menu",  tabPanel("output", output_ui(id = "1")))
			output_server(id = "1")
			saved_setting$value <- c(saved_setting$value, "output")
		}
	})

	observe({
		#req(length(working_project()) > 0)
		if (!("Setting" %in% saved_setting$value)) {
			saved_setting$value <- c(saved_setting$value, "Setting")
			insertTab(session=session, inputId = "menu", target = "Dataset",  position = "after", tabPanel("Setting",	uiOutput("settingControls")))
			output$settingControls <- renderUI({
				tags$div('id'="hexcontainer",
					tags$div('class'="hex1",
						tags$div('class'="hex2",
							tags$div('class'="hexlink", 'id'="hl1",
								if (!(modulelist[1] %in% saved_setting$value)) {
									tags$a('id'="add1", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Add QC", tags$br(), "Module"),
											tags$div('class'="plus")
										)
									)
								} else {
									tags$a('id'="remove1", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Remove QC", tags$br(), "Module"),
											tags$div('class'="minus")
										)
									)
								}
							)
						)
					),
					tags$div('class'="hex1",
						tags$div('class'="hex2",
							tags$div('class'="hexlink", 'id'="hl2",
								if (!(modulelist[2] %in% saved_setting$value)) {
									tags$a('id'="add2", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Add DEGs", tags$br(), "Module"),
											tags$div('class'="plus")
										)
									)
								} else {
									tags$a('id'="remove2", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Remove DEGs", tags$br(), "Module"),
											tags$div('class'="minus")
										)
									)
								}
							)
						)
					),
					tags$div('class'="hex1",
						tags$div('class'="hex2",
							tags$div('class'="hexlink", 'id'="hl3",
								if (!(modulelist[3] %in% saved_setting$value)) {
									tags$a('id'="add3", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Add Heatmap", tags$br(), "Module"),
											tags$div('class'="plus")
										)
									)
								} else {
									tags$a('id'="remove3", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Remove Heatmap", tags$br(), "Module"),
											tags$div('class'="minus")
										)
									)
								}
							)
						)
					),
					tags$div('class'="hex1",
						tags$div('class'="hex2",
							tags$div('class'="hexlink", 'id'="hl4",
								if (!(modulelist[4] %in% saved_setting$value)) {
									tags$a('id'="add4", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Add Expression", tags$br(), "Module"),
											tags$div('class'="plus")
										)
									)
								} else {
									tags$a('id'="remove4", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Remove Expression", tags$br(), "Module"),
											tags$div('class'="minus")
										)
									)
								}
							)
						)
					),
					tags$div('class'="hex1",
						tags$div('class'="hex2",
							tags$div('class'="hexlink", 'id'="hl5",
								if (!(modulelist[5] %in% saved_setting$value)) {
									tags$a('id'="add5", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Add Geneset", tags$br(), "Enrichment Module"),
											tags$div('class'="plus")
										)
									)
								} else {
									tags$a('id'="remove5", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Remove Geneset", tags$br(), "Enrichment Module"),
											tags$div('class'="minus")
										)
									)
								}
							)
						)
					),
tags$br(),
					tags$div('class'="hex1",
						tags$div('class'="hex2",
							tags$div('class'="hexlink", 'id'="hl6",
								if (!(modulelist[6] %in% saved_setting$value)) {
									tags$a('id'="add6", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Add Pattern", tags$br(), "Module"),
											tags$div('class'="plus")
										)
									)
								} else {
									tags$a('id'="remove6", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Remove Pattern", tags$br(), "Module"),
											tags$div('class'="minus")
										)
									)
								}
							)
						)
					),
					tags$div('class'="hex1",
						tags$div('class'="hex2",
							tags$div('class'="hexlink", 'id'="hl10",
								if (!(modulelist[10] %in% saved_setting$value)) {
									tags$a('id'="add10", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Add WGCNA", tags$br(), "Module"),
											tags$div('class'="plus")
										)
									)
								} else {
									tags$a('id'="remove10", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Remove WGCNA", tags$br(), "Module"),
											tags$div('class'="minus")
										)
									)
								}
							)
						)
					),
					tags$div('class'="hex1",
						tags$div('class'="hex2",
							tags$div('class'="hexlink", 'id'="hl7",
								if (!(modulelist[7] %in% saved_setting$value)) {
									tags$a('id'="add7", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Add Correlation", tags$br(), "Network Module"),
											tags$div('class'="plus")
										)
									)
								} else {
									tags$a('id'="remove7", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Remove Correlation", tags$br(), "Network Module"),
											tags$div('class'="minus")
										)
									)
								}
							)
						)
					),
					tags$div('class'="hex1",
						tags$div('class'="hex2",
							tags$div('class'="hexlink", 'id'="hl8",
								if (!(modulelist[8] %in% saved_setting$value)) {
									tags$a('id'="add8", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Add Venn", tags$br(), "Diagram Module"),
											tags$div('class'="plus")
										)
									)
								} else {
									tags$a('id'="remove8", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Remove Venn", tags$br(), "Diagram Module"),
											tags$div('class'="minus")
										)
									)
								}
							)
						)
					),
tags$br(),
					tags$div('class'="hex1",
						tags$div('class'="hex2",
							tags$div('class'="hexlink", 'id'="hl9",
								if (!(modulelist[9] %in% saved_setting$value)) {
									tags$a('id'="add9", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Add PCSF", tags$br(), "Module"),
											tags$div('class'="plus")
										)
									)
								} else {
									tags$a('id'="remove9", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Remove PCSF", tags$br(), "Module"),
											tags$div('class'="minus")
										)
									)
								}
							)
						)
					),
					tags$div('class'="hex1",
						tags$div('class'="hex2",
							tags$div('class'="hexlink", 'id'="hl11",
								if (!(modulelist[11] %in% saved_setting$value)) {
									tags$a('id'="add11", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Add Curve Fitting ", tags$br(), "Dose Response"),
											tags$div('class'="plus")
										)
									)
								} else {
									tags$a('id'="remove11", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Remove Curve Fitting", tags$br(), "Dose Response"),
											tags$div('class'="minus")
										)
									)
								}
							)
						)
					),
					tags$div('class'="hex1",
						tags$div('class'="hex2",
							tags$div('class'="hexlink", 'id'="hl12",
								if (!(modulelist[12] %in% saved_setting$value)) {
									tags$a('id'="add12", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Add Curve Fitting", tags$br(), "Time Course"),
											tags$div('class'="plus")
										)
									)
								} else {
									tags$a('id'="remove12", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Remove Curve Fitting", tags$br(), "Time Course"),
											tags$div('class'="minus")
										)
									)
								}
							)
						)
					),
					tags$div('class'="hex1",
						tags$div('class'="hex2",
							tags$div('class'="hexlink", 'id'="hl13",
								if (!(modulelist[13] %in% saved_setting$value)) {
									tags$a('id'="add13", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Add Monotonic", tags$br(), "Trend Module"),
											tags$div('class'="plus")
										)
									)
								} else {
									tags$a('id'="remove13", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Remove Monotonic", tags$br(), "Trend Module"),
											tags$div('class'="minus")
										)
									)
								}
							)
						)
					),
					tags$div('class'="hex1",
						tags$div('class'="hex2",
							tags$div('class'="hexlink", 'id'="hl14",
								if (!(modulelist[14] %in% saved_setting$value)) {
									tags$a('id'="add14", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Add Merge Data", tags$br(), "Module"),
											tags$div('class'="plus")
										)
									)
								} else {
									tags$a('id'="remove14", 'class'="btn btn-secondary action-button",
										tags$div('class'="hexcover",
											tags$h3("Remove Merge Data", tags$br(), "Module"),
											tags$div('class'="minus")
										)
									)
								}
							)
						)
					)#,
					#tags$div('class'="hex1",
					#	tags$div('class'="hex2",
					#		tags$div('class'="hexlink", 'id'="hl15",
					#			if (!(modulelist[15] %in% saved_setting$value)) {
					#				tags$a('id'="add15", 'class'="btn btn-secondary action-button",
					#					tags$div('class'="hexcover",
					#						tags$h3("Add Future", tags$br(), "Module"),
					#						tags$div('class'="plus")
					#					)
					#				)
					#			} else {
					#				tags$a('id'="remove15", 'class'="btn btn-secondary action-button",
					#					tags$div('class'="hexcover",
					#						tags$h3("Remove Future", tags$br(), "Module"),
					#						tags$div('class'="minus")
					#					)
					#				)
					#			}
					#		)
					#	)
					#)
				)
			})
		}
	})

	####
	observeEvent(input$add1, {
		if (modulelist[1] %in% saved_setting$value) {
			updateTabsetPanel(session=session, "menu", selected = modulelist[1])
		} else {
			saved_setting$value <- c(saved_setting$value, modulelist[1])
			source(moduleFilelist[1],  local = TRUE)
			insertTab(session=session, inputId = "menu", target = "Groups and Samples",  position = "after", tabPanel(modulelist[1],	qcplot_ui("1")))
			qcplot_server(id = "1")
		}
	})

	observeEvent(input$remove1, {
		removeTab(session=session, inputId = "menu", target = modulelist[1])
		saved_setting$value = saved_setting$value[saved_setting$value!= modulelist[1]]
	})

	###
	observeEvent(input$add2, {
		if (modulelist[2] %in% saved_setting$value) {
			updateTabsetPanel(session=session, "menu", selected = modulelist[2])
		} else {
			saved_setting$value <- c(saved_setting$value, modulelist[2])
			source(moduleFilelist[2],  local = TRUE)
			insertTab(session=session, inputId = "menu", target = "Groups and Samples",  position = "after", tabPanel(modulelist[2],	deg_ui("deg")))
			deg_server(id = "deg")
		}
	})

	observeEvent(input$remove2, {
		removeTab(session=session, inputId = "menu", target = modulelist[2])
		saved_setting$value = saved_setting$value[saved_setting$value!= modulelist[2]]
	})

	###
	observeEvent(input$add3, {
		if (modulelist[3] %in% saved_setting$value) {
			updateTabsetPanel(session=session,  "menu", selected = modulelist[3])
		} else {
			saved_setting$value <- c(saved_setting$value, modulelist[3])
			source(moduleFilelist[3],  local = TRUE)
			insertTab(session=session, inputId = "menu", target = "Groups and Samples",  position = "after", tabPanel(modulelist[3],	heatmap_ui("3")))
			heatmap_server(id = "3")
		}
	})

	observeEvent(input$remove3, {
		removeTab(session=session, inputId = "menu", target = modulelist[3])
		saved_setting$value = saved_setting$value[saved_setting$value!= modulelist[3]]
	})

	###
	observeEvent(input$add4, {
		if (modulelist[4] %in% saved_setting$value) {
			updateTabsetPanel(session=session,  "menu", selected = modulelist[4])
		} else {
			saved_setting$value <- c(saved_setting$value, modulelist[4])
			source(moduleFilelist[4], local = TRUE)
			insertTab(session=session,  inputId = "menu", target = "Groups and Samples",  position = "after", tabPanel(modulelist[4],	expression_ui("4")))
			expression_server(id = "4")
		}
	})

	observeEvent(input$remove4, {
		removeTab(session=session, inputId = "menu", target = modulelist[4])
		saved_setting$value = saved_setting$value[saved_setting$value!= modulelist[4]]
	})

	###
	observeEvent(input$add5, {
		if (modulelist[5] %in% saved_setting$value) {
			updateTabsetPanel(session=session,  "menu", selected = modulelist[5])
		} else {
			saved_setting$value <- c(saved_setting$value, modulelist[5])
			source(moduleFilelist[5],  local = TRUE)
			insertTab(session=session, inputId = "menu", target = "Groups and Samples",  position = "after", tabPanel(modulelist[5],	geneset_ui("5")))
			geneset_server(id = "5")
		}
	})

	observeEvent(input$remove5, {
		removeTab(session=session, inputId = "menu", target = modulelist[5])
		saved_setting$value = saved_setting$value[saved_setting$value!= modulelist[5]]
	})

	###
	observeEvent(input$add6, {
		if (modulelist[6] %in% saved_setting$value) {
			updateTabsetPanel(session=session,  "menu", selected = modulelist[6])
		} else {
			saved_setting$value <- c(saved_setting$value, modulelist[6])
			source(moduleFilelist[6],  local = TRUE)
			insertTab(session=session, inputId = "menu", target = "Groups and Samples",  position = "after", tabPanel(modulelist[6],	pattern_ui("6")))
			pattern_server(id = "6")
		}
	})

	observeEvent(input$remove6, {
		removeTab(session=session, inputId = "menu", target = modulelist[6])
		saved_setting$value = saved_setting$value[saved_setting$value!= modulelist[6]]
	})

	###
	observeEvent(input$add7, {
		if (modulelist[7] %in% saved_setting$value) {
			updateTabsetPanel(session=session,  "menu", selected = modulelist[7])
		} else {
			saved_setting$value <- c(saved_setting$value, modulelist[7])
			source(moduleFilelist[7],  local = TRUE)
			insertTab(session=session,  inputId = "menu", target = "Groups and Samples",  position = "after", tabPanel(modulelist[7],	network_ui("7")))
			network_server(id = "7")
		}
	})

	observeEvent(input$remove7, {
		removeTab(session=session, inputId = "menu", target = modulelist[7])
		saved_setting$value = saved_setting$value[saved_setting$value!= modulelist[7]]
	})

	###
	observeEvent(input$add8, {
		if (modulelist[8] %in% saved_setting$value) {
			updateTabsetPanel(session=session,  "menu", selected = modulelist[8])
		} else {
			saved_setting$value <- c(saved_setting$value, modulelist[8])
			source(moduleFilelist[8],  local = TRUE)
			insertTab(session=session,  inputId = "menu", target = "Groups and Samples",  position = "after", tabPanel(modulelist[8],	venn_ui("8")))
			venn_server(id = "8")
		}
	})

	observeEvent(input$remove8, {
		removeTab(session=session, inputId = "menu", target = modulelist[8])
		saved_setting$value = saved_setting$value[saved_setting$value!= modulelist[8]]
	})

	###
	observeEvent(input$add9, {
		if (modulelist[9] %in% saved_setting$value) {
			updateTabsetPanel(session=session,  "menu", selected = modulelist[9])
		} else {
			saved_setting$value <- c(saved_setting$value, modulelist[9])
			source(moduleFilelist[9],  local = TRUE)
			insertTab(session=session,  inputId = "menu", target = "Groups and Samples",  position = "after", tabPanel(modulelist[9],	pcsf_ui("9")))
			pcsf_server(id = "9")
		}
	})

	observeEvent(input$remove9, {
		removeTab(session=session, inputId = "menu", target = modulelist[9])
		saved_setting$value = saved_setting$value[saved_setting$value!= modulelist[9]]
	})

	###
	observeEvent(input$add10, {
		if (modulelist[10] %in% saved_setting$value) {
			updateTabsetPanel(session=session,  "menu", selected = modulelist[10])
		} else {
			saved_setting$value <- c(saved_setting$value, modulelist[10])
			source(moduleFilelist[10],  local = TRUE)
			insertTab(session=session,  inputId = "menu", target = "Groups and Samples",  position = "after", tabPanel(modulelist[10],	wgcna_ui("10")))
			wgcna_server(id = "10")
		}
	})

	observeEvent(input$remove10, {
		removeTab(session=session, inputId = "menu", target = modulelist[10])
		saved_setting$value = saved_setting$value[saved_setting$value!= modulelist[10]]
	})

	###
	observeEvent(input$add11, {
		if (modulelist[11] %in% saved_setting$value) {
			updateTabsetPanel(session=session,  "menu", selected = modulelist[11])
		} else {
			saved_setting$value <- c(saved_setting$value, modulelist[11])
			source(moduleFilelist[11],  local = TRUE)
			insertTab(session=session,  inputId = "menu", target = "Groups and Samples",  position = "after", tabPanel(modulelist[11],	drc_ui("11")))
			drc_server(id = "11")
		}
	})

	observeEvent(input$remove11, {
		removeTab(session=session, inputId = "menu", target = modulelist[11])
		saved_setting$value = saved_setting$value[saved_setting$value!= modulelist[11]]
	})

	###
	observeEvent(input$add12, {
		if (modulelist[12] %in% saved_setting$value) {
			updateTabsetPanel(session=session,  "menu", selected = modulelist[12])
		} else {
			saved_setting$value <- c(saved_setting$value, modulelist[12])
			source("util-basicandfitfunc.R",local=TRUE)$value
			source(moduleFilelist[12],  local = TRUE)
			insertTab(session=session, inputId = "menu", target = "Groups and Samples",  position = "after", tabPanel(modulelist[12],	dromics_ui("12")))
			dromics_server(id = "12")
		}
	})

	observeEvent(input$remove12, {
		removeTab(session=session, inputId = "menu", target = modulelist[12])
		saved_setting$value = saved_setting$value[saved_setting$value!= modulelist[12]]
	})

	###
	observeEvent(input$add13, {
		if (modulelist[13] %in% saved_setting$value) {
			updateTabsetPanel(session=session,  "menu", selected = modulelist[13])
		} else {
			saved_setting$value <- c(saved_setting$value, modulelist[13])
			source(moduleFilelist[13],  local = TRUE)
			insertTab(session=session,  inputId = "menu", target = "Groups and Samples",  position = "after", tabPanel(modulelist[13],	monotonic_ui("13")))
			monotonic_server(id = "13")
		}
	})

	observeEvent(input$remove13, {
		removeTab(session=session, inputId = "menu", target = modulelist[13])
		saved_setting$value = saved_setting$value[saved_setting$value!= modulelist[13]]
	})

	###
	observeEvent(input$add14, {
		if (modulelist[14] %in% saved_setting$value) {
			updateTabsetPanel(session=session,  "menu", selected = modulelist[14])
		} else {
			saved_setting$value <- c(saved_setting$value, modulelist[14])
			source(moduleFilelist[14],  local = TRUE)
			insertTab(session=session,  inputId = "menu", target = "Groups and Samples",  position = "after", tabPanel(modulelist[14],	mergedata_ui("14")))
			mergedata_server(id = "14")
		}
	})

	observeEvent(input$remove14, {
		removeTab(session=session, inputId = "menu", target = modulelist[14])
		saved_setting$value = saved_setting$value[saved_setting$value!= modulelist[14]]
	})
}

shinyApp(ui, server)





