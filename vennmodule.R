###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: venn.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 05/23/2023
##@version 3.0

##########################################################################################################
## Venn Diagram
##########################################################################################################
#pkgs: "VennDiagram", "ComplexHeatmap", "gplots", "colourpicker",  "DT", "shinyjqui", "rlang", "dplyr", "stringr",  "futile.logger",  "tibble"

library(VennDiagram)
library(ComplexHeatmap)
library(gplots)

venn_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
	  rclipboard::rclipboardSetup(),
		column(3,
			wellPanel(
				uiOutput(ns('loaddatasets')),
				fluidRow(
					column(width=6, numericInput(ns("venn_fccut"), label= "Fold Change Cutoff", value = 1.2, min=1, step=0.1)),
					column(width=6, numericInput(ns("venn_pvalcut"), label= "P Value Cutoff", value=0.01, min=0, step=0.001))
				),
				radioButtons(ns("venn_psel"), label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"), inline = TRUE),
				radioButtons(ns("venn_updown"), label= "All, Up or Down?", choices= c("All"="All","Up"="Up","Down"="Down"), inline = TRUE),
				radioButtons(ns("sel_label"), label= "Label name", choices= c("Gene.Name"="Gene.Name", "UniqueID"="UniqueID"), inline = TRUE, selected = "Gene.Name"),
				conditionalPanel(ns = ns, "input.venn_tabset=='UpSets' || input.venn_tabset=='UpSets Table'",
					radioButtons(ns("combinationmode"), label="Combination Mode", inline = TRUE,	choices = c("intersect" = "intersect", "distinct" = "distinct", "union" = "union"), selected = "intersect"),
					radioButtons(ns("combinationtranspose"), label="Transposing Combination Matrix", inline = TRUE,	choices = c("no" = "no", "yes" = "yes"), selected = "no")
				)
			)
		),
		column(9,
			tabsetPanel(id=ns("venn_tabset"),
				tabPanel(title="vennDiagram", value = "vennDiagram",
					fluidRow(
						column(9, align = "center", actionButton(ns("vennDiagram"), "Save to output"),
							plotOutput(ns("vennDiagram"), height = 600, width = 700)
						),
						column(3,
							wellPanel(
								textInput(ns("venn_title"), "Title", width = "100%"),
								fluidRow(
									column(width=6, sliderInput(ns("venn_mainx"), "Title X Pos", min = 0.1, max = 1, step=0.1,value = 0.5)),
									column(width=6, sliderInput(ns("venn_mainy"), "Title Y Pos", min = 0.8, max = 1.2, step=0.05, value = 1.05))
								),
								fluidRow(
									column(width=6, sliderInput(ns("venn_maincex"), "Title Size", min = 0, max = 6, step=1, value = 3)),
									column(width=6, sliderInput(ns("venn_margin"), "Margin", min = 0.05, max = 0.3, step=0.05, value = 0.1))
								),
								fluidRow(
									column(width=6, sliderInput(ns("venn_alpha"), "Opacity", min = 0, max = 1, value = 0.4)),
									column(width=6, sliderInput(ns("venn_rotation"), "Rotation Degree", min = -180, max = 180, step = 30, value = 0))
								),
								fluidRow(
									column(width=6, sliderInput(ns("venn_lwd"), "Line Thick", min = 1, max = 4, value = 1)),
									column(width=6, sliderInput(ns("venn_lty"), "Line Type", min = 1, max = 6, value = 1))
								),
								fluidRow(
									column(width=6, radioButtons(ns("venn_fontface"),"Number Font Face",list("plain", "bold"),selected = "plain",	inline = TRUE)),
									column(width=6, sliderInput(ns("venn_cex"), "Font size", min = 1, max = 4, value = 2))
								),
								fluidRow(
									column(width=6, radioButtons(ns("venn_catfontface"),"Label Font Face",list("plain", "bold"),	selected = "plain", inline = TRUE)),
									column(width=6, sliderInput(ns("venn_catcex"), "Font Size", min = 0.5, max = 3, step=0.5, value = 2))
								)
							)
						)
					),
					hr(),
					fluidRow(
						column(2,
							wellPanel(
								fluidRow(
									column(width=8,textInput(ns("venn_cat1"), "Rename:", width = "100%")),
									column(width=4,colourpicker::colourInput(ns("col1"), "Color", "#0000FF",palette = "limited"))
								),
								sliderInput(ns("venn_pos1"), "Position:", min = -360, max = 360, step = 30, value = 0),
								sliderInput(ns("venn_dis1"), "Distance:", min = -0.5, max =0.5, step = 0.01, value = 0.2)
							)
						),
						column(2,
							wellPanel(
								fluidRow(
									column(width=8,textInput(ns("venn_cat2"), "Rename:", width = "100%")),
									column(width=4,colourpicker::colourInput(ns("col2"), "Color", "#FF7F00",palette = "limited"))
								),
								sliderInput(ns("venn_pos2"), "Position:", min = -360, max = 360, step = 30, value = 287.5),
								sliderInput(ns("venn_dis2"), "Distance:", min = -0.5, max =0.5, step = 0.01, value =0.2)
							)
						),
						column(2,
							wellPanel(
								fluidRow(
									column(width=8,textInput(ns("venn_cat3"), "Rename:", width = "100%")),
									column(width=4,colourpicker::colourInput(ns("col3"), "Color", "#00FF00",palette = "limited"))
								),
								sliderInput(ns("venn_pos3"), "Position:", min = -360, max = 360, step = 30, value = 215),
								sliderInput(ns("venn_dis3"), "Distance:", min = -0.5, max =0.5, step = 0.01, value = 0.2)
							)
						),
						column(2,
							wellPanel(
								fluidRow(
									column(width=8,textInput(ns("venn_cat4"), "Rename:", width = "100%")),
									column(width=4,colourpicker::colourInput(ns("col4"), "Color", "#FF00FF",palette = "limited"))
								),
								sliderInput(ns("venn_pos4"), "Position:", min = -360, max = 360, step = 30, value = 145),
								sliderInput(ns("venn_dis4"), "Distance:", min = -0.5, max =0.5, step = 0.01, value = 0.2)
							)
						),
						column(2,
							wellPanel(
								fluidRow(
									column(width=8,textInput(ns("venn_cat5"), "Rename:", width = "100%")),
									column(width=4,colourpicker::colourInput(ns("col5"), "Color", "#FFFF00",palette = "limited"))
								),
								sliderInput(ns("venn_pos5"), "Position:", min = -360, max = 360, step = 30, value = 70),
								sliderInput(ns("venn_dis5"), "Distance:", min = -0.5, max =0.5, step = 0.01, value = 0.2)
							)
						)
					)
				),
				#	tabPanel(title="VennDiagram(black & white)", plotOutput(ns("SvennDiagram"), height = 800, width = 800)),
				tabPanel(title="Intersection Output", htmlOutput(ns("vennHTML"))),
				tabPanel(title="UpSets", value = "UpSets",  actionButton(ns("upsetsplot"), "Save to output"), jqui_resizable(plotOutput(outputId = ns("upsetsplot"), height = 700, width = 1200))),
				tabPanel(title="UpSets Table", value = "UpSets Table", DT::dataTableOutput(ns("UpSetsList"))),
				tabPanel(title="Help", htmlOutput("help_venn"))
			)
		)
	)
}

venn_server <- function(id) {
	shiny::moduleServer(id,
		function(input, output, session) {
			ns <- session$ns

			output$loaddatasets <- renderUI({
				req(length(working_project()) > 0)
				projectlist <- list()
				for (project in DS_names()) {
					projectlist <-	append(projectlist, paste(project, DataInSets[[project]]$tests_order, sep="->"))
				}
				tagList(
					shinyjqui::sortableCheckboxGroupInput(ns("add_dataset"), label = "Drag to order, then select Datasets", choices = projectlist)#,
				)
			})

			maxcomparison <- reactiveVal(NULL)
			observeEvent(input$venn_tabset, {
				if (input$venn_tabset == "vennDiagram") {
					maxcomparison(5)
				} else {
					maxcomparison(12)
				}
			})

			observe({
				req(length(working_project()) > 0)
				req(input$add_dataset)
				projectlist <- list()
				for (project in DS_names()) {
					projectlist <-	append(projectlist, paste(project, DataInSets[[project]]$tests_order, sep="->"))
				}
				if(length(input$add_dataset) > maxcomparison()){
					updateCheckboxGroupInput(session, "add_dataset", selected= head(input$add_dataset, maxcomparison()))
				}
			})

			observe({
				req(input$add_dataset)
				selectdatasets <- input$add_dataset
				selectdatasetsdf <- str_split(selectdatasets, "->", simplify = TRUE) %>%
				as.data.frame() %>%
				rlang::set_names(c("project", "tests")) %>%
				dplyr::group_by(project) %>%
				dplyr::summarise(tests = list(tests))

				if (nrow(selectdatasetsdf)  == 1) {
					catnames  <- selectdatasetsdf$tests[[1]]
				} else {
					catnames <- c()
					selectdatasetsdf <- stringr::str_split(selectdatasets, "->", simplify = TRUE) %>%
					as.data.frame() %>%
					rlang::set_names(c("project", "tests"))
					for (row in 1:nrow(selectdatasetsdf)) {
						ProjectID <- selectdatasetsdf[row, "project"] %>% as.character()
						ShortName <- DataInSets[[ProjectID]]$ShortName
						testname <- selectdatasetsdf[row, "tests"]
						catname <- paste(ShortName, "->", testname, sep="")
						catnames  <- c(catnames, catname)
					}
				}

				updateTextInput(session, "venn_cat1", value=catnames[1])
				updateTextInput(session, "venn_cat2", value=catnames[2])
				updateTextInput(session, "venn_cat3", value=catnames[3])
				updateTextInput(session, "venn_cat4", value=catnames[4])
				updateTextInput(session, "venn_cat5", value=catnames[5])

				if (length(selectdatasets)==1) {
					updateSliderInput(session, "venn_pos1", value=0)
					updateSliderInput(session, "venn_dis1", value=0.025)
				}
				if (length(selectdatasets)==2) {
					updateSliderInput(session, "venn_pos1", value=-20)
					updateSliderInput(session, "venn_pos2", value=20)
					updateSliderInput(session, "venn_dis1", value=0.05)
					updateSliderInput(session, "venn_dis2", value=0.05)
				}
				if (length(selectdatasets)==3) {
					updateSliderInput(session, "venn_pos1", value=-40)
					updateSliderInput(session, "venn_pos2", value=40)
					updateSliderInput(session, "venn_pos3", value=180)
					updateSliderInput(session, "venn_dis1", value=0.05)
					updateSliderInput(session, "venn_dis2", value=0.05)
					updateSliderInput(session, "venn_dis3", value=0.025)
				}
				if (length(selectdatasets)==4) {
					updateSliderInput(session, "venn_pos1", value=-15)
					updateSliderInput(session, "venn_pos2", value=15)
					updateSliderInput(session, "venn_pos3", value=0)
					updateSliderInput(session, "venn_pos4", value=0)
					updateSliderInput(session, "venn_dis1", value=0.22)
					updateSliderInput(session, "venn_dis2", value=0.22)
					updateSliderInput(session, "venn_dis3", value=0.11)
					updateSliderInput(session, "venn_dis4", value=0.11)
				}
				if (length(selectdatasets)==5) {
					updateSliderInput(session, "venn_pos1", value=0)
					updateSliderInput(session, "venn_pos2", value=287)
					updateSliderInput(session, "venn_pos3", value=215)
					updateSliderInput(session, "venn_pos4", value=145)
					updateSliderInput(session, "venn_pos5", value=70)
					updateSliderInput(session, "venn_dis1", value=0.2)
					updateSliderInput(session, "venn_dis2", value=0.2)
					updateSliderInput(session, "venn_dis3", value=0.2)
					updateSliderInput(session, "venn_dis4", value=0.2)
					updateSliderInput(session, "venn_dis5", value=0.2)
				}
			})

			DataVennReactive <- reactive({

				req(input$add_dataset)
				selectdatasets <- input$add_dataset
				selectdatasetsdf <- str_split(selectdatasets, "->", simplify = TRUE) %>%
				as.data.frame() %>%
				rlang::set_names(c("project", "tests")) %>%
				dplyr::group_by(project) %>%
				dplyr::summarise(tests = list(tests))

				results_long_sel <- list()
				for (row in 1:nrow(selectdatasetsdf)) {
					ProjectID <- selectdatasetsdf[row, "project"] %>% as.character()
					testname <- selectdatasetsdf[row, "tests"][[1]][[1]]

					res = DataInSets[[ProjectID]]$results_long %>%
					dplyr::filter(test %in%testname) %>%
					dplyr::select(one_of(c("UniqueID","Gene.Name","test","logFC", "P.Value","Adj.P.Value")) ) %>%
					dplyr::mutate (ProjectID = ProjectID)
					results_long_sel[[row]] <- res
				}

				results_long  <- bind_rows(results_long_sel)  %>% dplyr::mutate_if(is.factor, as.character)

				test_sel <- NA
				p_sel <- input$venn_psel
				pcut = as.numeric(input$venn_pvalcut)
				FCcut = log2(as.numeric(input$venn_fccut))
				direction <- input$venn_updown
				sel_label <- input$sel_label
				resultfiltered <- resultfilter(results_long, test_sel, p_sel, direction, pcut, FCcut, sel_label)
				results_long <- resultfiltered[["results_long_filtered"]]

				vennlist <- fill <- venn_cat <- venn_pos <- venn_dis <- list()
				selectdatasets <- input$add_dataset
				for (i in 1:length(selectdatasets)) {
					selectdataset <- selectdatasets[i]
					ProjectIDname <- str_split(selectdataset, "->", simplify = TRUE)[,1]
					testname <-  str_split(selectdataset, "->", simplify = TRUE)[,2]
					labelidlist = results_long %>% as.data.frame() %>% dplyr::filter((test == testname) & (ProjectID == ProjectIDname)) %>% dplyr::pull(labelid) %>% unique()
					vennlist[[selectdataset]]  <- labelidlist
					fill[[selectdataset]] <- input[[paste0("col", i)]]
					venn_cat[[selectdataset]] <- input[[paste0("venn_cat", i)]]
					venn_pos[[selectdataset]] <- input[[paste0("venn_pos", i)]]
					venn_dis[[selectdataset]] <- input[[paste0("venn_dis", i)]]
				}
				return(venndata = list("vennlist"=vennlist, "fillcor"=fill, "venn_cat" = venn_cat, "venn_pos" = venn_pos, "venn_dis" = venn_dis ))
			})

			vennDiagram_out <- reactive({
				venndata <- DataVennReactive()
				vennlist <- venndata$vennlist
				fillcor <- unlist(venndata$fillcor)
				venn_cat <- unlist(venndata$venn_cat)
				venn_pos <- unlist(venndata$venn_pos)
				venn_dis <- unlist(venndata$venn_dis)

				SetNum = length(vennlist)
				futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
				venn.plot <- venn.diagram(
					#height = 3000, width = 3000,
					x = vennlist,
					fill = fillcor,
					lty = input$venn_lty,
					lwd = input$venn_lwd,
					alpha = input$venn_alpha,
					category.names = venn_cat,
					cex = input$venn_cex,
					fontface = input$venn_fontface,
					cat.cex = input$venn_catcex,
					cat.fontface = input$venn_catfontface,
					cat.dist = venn_dis,
					cat.pos = venn_pos,
					rotation.degree = input$venn_rotation,
					#reverse = FALSE,
					#euler.d = TRUE,
					margin = input$venn_margin,
					main = input$venn_title,
					main.cex = input$venn_maincex,
					main.pos = c(input$venn_mainx, input$venn_mainy),
					main.fontface = "bold",
					filename = NULL
				)

				return(venn.plot)
			})

			output$vennDiagram <- renderPlot({
				grid.draw(vennDiagram_out())
			})

			observeEvent(input$vennDiagram, {
				saved.num <- length(saved_plots$vennDiagram) + 1
				saved_plots$vennDiagram[[saved.num]] <- vennDiagram_out()
			})

			############################
			DataUpsetsReactive <- reactive({
				req(input$add_dataset)
				selectdatasets <- input$add_dataset
				selectdatasetsdf <- str_split(selectdatasets, "->", simplify = TRUE) %>%
				as.data.frame() %>%
				rlang::set_names(c("project", "tests")) %>%
				dplyr::group_by(project) %>%
				dplyr::summarise(tests = list(tests))

				results_long_sel <- list()
				for (row in 1:nrow(selectdatasetsdf)) {
					ProjectID <- selectdatasetsdf[row, "project"] %>% as.character()
					tests <- selectdatasetsdf[row, "tests"][[1]][[1]]

					res = DataInSets[[ProjectID]]$results_long %>%
					dplyr::filter(test %in% tests) %>%
					dplyr::select(one_of(c("UniqueID","Gene.Name","test","logFC", "P.Value","Adj.P.Value")) ) %>%
					dplyr::mutate (ProjectID = ProjectID)
					results_long_sel[[row]] <- res
				}

				results_long  <- bind_rows(results_long_sel)  %>% dplyr::mutate_if(is.factor, as.character)

				test_sel <- NA
				p_sel <- input$venn_psel
				pcut = as.numeric(input$venn_pvalcut)
				FCcut = log2(as.numeric(input$venn_fccut))
				direction <- input$venn_updown
				sel_label <- input$sel_label
				resultfiltered <- resultfilter(results_long, test_sel, p_sel, direction, pcut, FCcut, sel_label)
				results_long_filtered <- resultfiltered[["results_long_filtered"]]

				upsetslistdf <- results_long_filtered  %>% dplyr::group_by(ProjectID, test) %>%
				dplyr::summarise(labelids = list(unique(labelid)), .groups = "drop")

				vennlist <- list()
				for (row in 1:nrow(upsetslistdf)) {
					if (length(unique(upsetslistdf$ProjectID))  == 1)
					Upsetname <- upsetslistdf[row, "test"] %>% as.character()
					else {
						ProjectID <- upsetslistdf[row, "ProjectID"] %>% unique() %>% as.character()
						ShortName <- DataInSets[[ProjectID]]$ShortName
						Upsetname <- paste(ShortName, "->", upsetslistdf[row, "test"], sep="")
						#Upsetname <- paste(upsetslistdf[row, "ProjectID"], upsetslistdf[row, "test"],sep = "-")
					}
					vennlist[[Upsetname]] <- upsetslistdf[row, "labelids"][[1]][[1]]
				}
				#saveRDS(vennlist, "vennlist.rds")
				return(vennlist)
			})

			upsets_out <- reactive({
				vennlist <- DataUpsetsReactive()
				combinationmode <- input$combinationmode
				m = ComplexHeatmap::make_comb_mat(vennlist, mode = combinationmode)

				UpSetlist <- attributes(m)$data %>% as.data.frame() %>%
				dplyr::mutate(across(everything(), as.character)) %>%
				unite("Code", 1:ncol(.), sep="", na.rm = TRUE, remove = FALSE) %>%
				tibble::rownames_to_column(., var = "labelid") %>%
				dplyr::group_by(Code, .drop = FALSE) %>%
				dplyr::mutate(count = n(), labelids = list(unique(labelid))) %>%
				#dplyr::summarise(labelids = list(unique(labelid)), .groups = "drop") %>%
				dplyr::select(-labelid) %>%
				dplyr::ungroup() %>%
				dplyr::distinct()

				if (input$combinationtranspose == "yes")
				m <- t(m)

				p <- ComplexHeatmap::UpSet(m,
					comb_order = order(comb_size(m), decreasing = TRUE),
					comb_col = "#0000FF",
					bg_col = c("#F0F0FF", "#FFF0F0"),
					bg_pt_col = "#CCCCFF",
					top_annotation = upset_top_annotation(m, ylim = c(0, max(comb_size(m))*1.1), height = unit(4, "cm"), add_numbers = TRUE),
					right_annotation = upset_right_annotation(m, ylim = c(0, max(set_size(m))*1.1), gp = gpar(fill = "black"),	width = unit(4, "cm"), add_numbers = TRUE)
				)
				return(list("p" = p, "UpSetlist" = UpSetlist))
			})

			output$upsetsplot <- renderPlot({
				upsets_out()$p
			})

			output$UpSetsList <- DT::renderDT(server=FALSE,{
				result_long_tmp <- upsets_out()$UpSetlist
				result_long_tmp[,sapply(result_long_tmp,is.numeric)] <- signif(result_long_tmp[,sapply(result_long_tmp,is.numeric)],3)
				
				# Convert a list of strings to a vector, and then collapse the vector into a single string
				# Otherwise, "copy" will copy a list of strings with quotation marks
				result_long_tmp$Label_IDs <- NA
				
				for (i in 1:nrow(result_long_tmp)) {
				  result_long_tmp$Label_IDs[i] = paste0(unlist(result_long_tmp$labelids[i]), collapse = ",")
				}
				
				result_long_tmp$Action <- vapply(1L:nrow(result_long_tmp), function(i){
				  as.character(
				    rclipButton(
				      paste0("clipbtn_", i), 
				      label = "Copy", 
				      clipText = result_long_tmp[i, "Label_IDs"], 
				      #icon = icon("clipboard"),
				      icon = icon("copy", lib = "glyphicon"),
				      class = "btn-primary btn-sm"
				    )
				  )
				}, character(1L))
				

				
				# re-arrange columns
				result_long_tmp <- result_long_tmp %>% 
				  dplyr::select(-labelids) %>% 
				  dplyr::relocate(Action, .before = Label_IDs)
				
				DT::datatable(result_long_tmp,  extensions = 'Buttons', escape = FALSE, selection = "none", options = list(
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

			observeEvent(input$upsetsplot, {
				saved.num <- length(saved_plots$upsetsplot) + 1
				saved_plots$upsetsplot[[saved.num]] <- upsets_out()$p
			})

			######################################
			output$SvennDiagram <- renderPlot({
				venndata <- DataVennReactive()
				vennlist <- venndata$vennlist
				venn(vennlist, show.plot = TRUE, intersections = FALSE)
			})

			output$vennHTML <- renderText({
				req(length(working_project()) > 0)
				ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
				venndata <- DataVennReactive()
				vennlist <- venndata$vennlist
				v.table <- venn(vennlist,show.plot = FALSE, intersections = TRUE)
				intersect <- attr(v.table,"intersections")
				htmlstr <- "  <br>"
				for (i in 1:length(intersect)) {
					#if(input$vennlistname == "Gene.Name"){
					#	genes <- ProteinGeneName%>%dplyr::filter(UniqueID %in% intersect[[i]])%>%dplyr::select(Gene.Name)%>%unlist()
					#	intersectlist <- toString(genes)
					#	#intersectlist <- toString(sapply(strsplit(intersect[[i]],split= "\\_"),'[',1))
					#} else if (input$vennlistname == "Protein.ID") {
					#	intersectlist <- toString(sapply(strsplit(intersect[[i]],split= "\\_"),'[',2))
					#} else if (input$vennlistname == "UniqueID") {
					intersectlist <- toString(intersect[[i]])
					#}
					htmlstr <-  paste(htmlstr,"<p><b><font color=red>", names(intersect[i]),"</font></b>:",intersectlist, sep="")
				}
				htmlstr
			})
		}
	)
}
