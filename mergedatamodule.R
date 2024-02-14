###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: mergedatamodule.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 07/23/2023
##@version 1.0
###########################################################################################################
#pkgs:"DT", "shinyjqui", "dplyr","stringr", "rlang"
#data req: tests_order, ProteinGeneNameHeader, results_long

fisher_pvalue = function(marker_p_values){
	marker_p_values = data.frame(marker_p_values[!is.na(marker_p_values)])
	new_T = -2*sum(log(marker_p_values), na.rm = T)
	new_p = pchisq(new_T, df=2*nrow(marker_p_values), lower.tail=FALSE)
	if(new_p==0){
		new_p=1e-323
	}
	new_p
}

minP_pvalue = function(marker_p_values){
	T_min = min(marker_p_values, na.rm = T)
	p_min = pbeta(T_min, 1, length(marker_p_values))
	p_min
}

simes_pvalue = function(marker_p_values){
	marker_p_values = data.frame(marker_p_values[!is.na(marker_p_values)])
	r=rank(marker_p_values)
	pval=min(nrow(marker_p_values)*marker_p_values/r)
}

stouffer_pvalue = function(marker_p_values){
	marker_p_values = (marker_p_values[!is.na(marker_p_values)])
	w = rep(1, length(marker_p_values))
	zp <- qnorm(marker_p_values, lower.tail = F, log.p = F) %*%
	w / sqrt(sum(w^2))
	p.val = pnorm(zp, lower.tail = F, log.p = F)
	return(p.val)
}

mergedata_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
		column(3,
			wellPanel(
				conditionalPanel(ns = ns, "input.tabset =='DEG Counts'",
					uiOutput(ns('loadedprojects')),
					column(width=6, numericInput(ns("FCcut"), label= "Fold Change Cutoff", value = 1.2, min=1, step=0.1)),
					column(width=6, numericInput(ns("pvalcut"), label= "P Value Cutoff", value=0.01, min=0, step=0.001)),
					radioButtons(ns("psel"), label= "P value or P.adj Value?", choices= c("Pval"="Pval", "Padj"="Padj"), inline = TRUE)
				),
				conditionalPanel(ns = ns, "input.tabset =='Merge Data'",
					uiOutput(ns('loaddatasets')),
					span("Merging Option", style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
					radioButtons(ns("merged_by"),label="Merged By", inline = TRUE, choices=""),
					textInput(ns("merged_name"), "Merged Data Name", value = "MergedData"),
					numericInput(ns("overlapnum"), label= "Overlapped Datasets", value=1, min=1, max = 2, step=1),
					radioButtons(ns("intersect"), label= "Union, Overlap or Distinct?", choices= c("Union"="Union", "Overlap"="Overlap", "Distinct"="Distinct", "Append" = "Append"), inline = TRUE, selected = "Union"),
					radioButtons(ns("pmergemethod"), label= "For common IDs, select P value merged methods", choices= c("fisher"="fisher_Pvalue", "minP"="minP_pvalue", "simes"="simes_pvalue","stouffer"="stouffer_pvalue"), inline = TRUE, selected = "fisher_Pvalue"),
					radioButtons(ns("fcmergemethod"), label= "For common IDs, select Fold Change merged methods", choices= c("Average"="Average", "Multiplication"="Multiplication", "Formula"="Formula"), inline = TRUE, selected = "Average"),
					conditionalPanel(ns = ns, "input.fcmergemethod=='Multiplication'",
						numericInput(ns("multiplicationfactor"), label= "Multiplication Factor", value=1, min=1, max = 2, step=0.1)
					),
					conditionalPanel(ns = ns, "input.fcmergemethod=='Formula'",
						numericInput(ns("intercept"), label= "intercept", value=0),
						numericInput(ns("slope"), label= "slope",value=1)
					)
				)
			)
		),
		column(9,
			tabsetPanel(id=ns("tabset"),
				tabPanel(title="DEG Counts", value ="DEG Counts", tags$hr(style="border-color: black;"), DT::dataTableOutput(ns("deg_counts"))
				),
				tabPanel(title="Merge Data", value = "Merge Data",
					actionButton(ns("merge_data"), "Merge Selected Data Sets"),
					shiny::downloadButton(outputId = ns("download_mergeddata_button"),  label = "Download All as CSV file"),
					DT::dataTableOutput(ns("mergedatatable"))
				),
				tabPanel(title="Data Source", value = "Data Source",
					shiny::downloadButton(outputId = ns("download_datasource_button"),  label = "Download All as CSV file"),
					DT::dataTableOutput(ns("datasourcetable"))
				),
				tabPanel(title="Help", value = "Help", htmlOutput("help_MergedData")
				)
			)
		)
	)
}

mergedata_server <- function(id) {
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

			output$loaddatasets <- renderUI({
				req(length(working_project()) > 0)
				projectlist <- list()
				for (project in names(DataInSets)) {
					projectlist <-	append(projectlist, paste(project, DataInSets[[project]]$tests_order, sep="->"))
				}
				tagList(
					shinyjqui::sortableCheckboxGroupInput(ns("add_dataset"), label = "Add Dataset to Compare", choices = projectlist)
				)
			})

			maxcomparison <- reactiveVal(NULL)
			observeEvent(input$tabset, {
				maxcomparison(8)
			})

			observe({
				req(length(working_project()) > 0)
				req(input$tabset)
				req(input$add_dataset)
				updateNumericInput(session, "overlapnum", max = length(input$add_dataset))

				projectlist <- list()
				for (project in names(DataInSets)) {
					projectlist <-	append(projectlist, paste(project, DataInSets[[project]]$tests_order, sep="->"))
				}
				if(length(input$add_dataset) > maxcomparison()){
					updateCheckboxGroupInput(session, "add_dataset", selected= head(input$add_dataset, maxcomparison()))
				}
			})

			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$tests_order)
				tests = DataInSets[[working_project()]]$tests_order
				ProteinGeneNameHeader = DataInSets[[working_project()]]$ProteinGeneNameHeader
				updateRadioButtons(session,'merged_by', inline = TRUE, choices=ProteinGeneNameHeader, selected="Gene.Name")
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
				DT::datatable(deg_counts_data(), extensions = 'Buttons',
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

			###merge data
			DataMergeReactive <- eventReactive(input$merge_data,{
				withProgress(message = 'Processing. It may take a while...', value = 0, {
					req(length(working_project()) > 0)

					selectdatasets <- input$add_dataset
					validate(need(length(input$add_dataset) >= 2, message = "Please Select >=2 sets"))

					selectdatasets <- input$add_dataset
					genelabel <- input$merged_by
					pmergemethod <- input$pmergemethod
					fcmergemethod <- input$fcmergemethod

					merged_name <- input$merged_name
					overlapnum <- as.numeric(input$overlapnum)

					selectdatasetsdf <- stringr::str_split(selectdatasets, "->", simplify = TRUE) %>%
					as.data.frame() %>%
					rlang::set_names(c("project", "tests")) %>%
					dplyr::group_by(project) %>%
					dplyr::summarise(tests = list(tests))

					results_long_sel <- list()
					ProteinGeneName_merge <- list()
					for (row in 1:nrow(selectdatasetsdf)) {
						ProjectID <- selectdatasetsdf[row, "project"] %>% as.character()
						tests <- selectdatasetsdf[row, "tests"][[1]][[1]]

						res = DataInSets[[ProjectID]]$results_long %>%
						dplyr::filter(!is.na(P.Value)) %>%
						dplyr::filter(test %in% tests) %>%
						#dplyr::mutate(labelid = !!sym(genelabel))  %>%
						dplyr::select(c(!!sym(genelabel), test, logFC, P.Value, Adj.P.Value)) %>%
						dplyr:: group_by(!!sym(genelabel), test) %>%
						dplyr::slice_min(order_by = P.Value, n = 1) %>%
						dplyr::mutate (ProjectID = ProjectID)  %>%
						dplyr::ungroup()

						results_long_sel[[row]] <- res
						ProteinGeneName_merge[[row]] <- DataInSets[[ProjectID]]$ProteinGeneName %>%
						dplyr::select(one_of(c("UniqueID",   "Gene.Name",  "Protein.ID")))
					}

					selectdataset2 <- stringr::str_split(selectdatasets[2], "->", simplify = TRUE) %>%
					as.data.frame() %>%
					rlang::set_names(c("project", "tests"))

					res_wide <- dplyr::bind_rows(results_long_sel) %>%
					tidyr::pivot_wider(id_cols= !!sym(genelabel),
						names_from = c(ProjectID,'test'),
						values_from = c(logFC, P.Value, Adj.P.Value),
					names_glue = '{.value}_({ProjectID}.{test})')

					if(fcmergemethod == "Average") {
						res <- dplyr::bind_rows(results_long_sel) %>%
						dplyr::select(-c(ProjectID, test))  %>%
						dplyr::filter(!is.na(P.Value)) %>%
						dplyr::group_by(!!sym(genelabel))  %>%
						dplyr::mutate(mean_logFC = log2(mean(2^logFC, na.rm=TRUE)), n = n())
					}

					if(fcmergemethod == "Multiplication") {
						multiplicationfactor <- as.numeric(input$multiplicationfactor)
						res <- dplyr::bind_rows(results_long_sel) %>%
						dplyr::mutate(logFC = ifelse((ProjectID == (selectdataset2 %>% dplyr::pull(project)) & test == (selectdataset2 %>% dplyr::pull(tests))), logFC * multiplicationfactor, logFC )) %>%
						dplyr::select(-c(ProjectID, test)) %>%
						dplyr::filter(!is.na(P.Value)) %>%
						dplyr::group_by(!!sym(genelabel)) %>%
						dplyr::mutate(mean_logFC = log2(mean(2^logFC, na.rm=TRUE)), n = n())
					}

					if(fcmergemethod == "Formula") {
						intercept <- as.numeric(input$intercept)
						slope <- as.numeric(input$slope)

						res <- dplyr::bind_rows(results_long_sel) %>%
						dplyr::mutate(logFC = ifelse((ProjectID == (selectdataset2 %>% dplyr::pull(project)) & test == (selectdataset2 %>% dplyr::pull(tests))), (intercept  + slope * logFC), logFC )) %>%
						dplyr::select(-c(ProjectID, test)) %>%
						dplyr::filter(!is.na(P.Value)) %>%
						dplyr::group_by(!!sym(genelabel)) %>%
						dplyr::mutate(mean_logFC = log2(mean(2^logFC, na.rm=TRUE)), n = n())
					}

					if(pmergemethod == "fisher_Pvalue")
					res <- res %>%	dplyr::mutate(merged_Pvalue = if_else(n >= 2, fisher_pvalue(P.Value), P.Value))
					if(pmergemethod == "minP_pvalue")
					res <- res %>%	dplyr::mutate(merged_Pvalue = if_else(n >= 2, minP_pvalue(P.Value), P.Value))
					if(pmergemethod == "simes_pvalue")
					res <- res %>%	dplyr::mutate(merged_Pvalue = if_else(n >= 2, simes_pvalue(P.Value), P.Value))
					if(pmergemethod == "stouffer_pvalue")
					res <- res %>%	dplyr::mutate(merged_Pvalue = if_else(n >= 2, stouffer_pvalue(P.Value), P.Value))

					res <- res %>% dplyr::select(-c(logFC, P.Value, Adj.P.Value)) %>%
					dplyr::filter(n >= overlapnum) %>%
					dplyr::distinct() %>%
					dplyr::ungroup() %>%
					dplyr::rename(P.Value = merged_Pvalue, logFC = mean_logFC) %>%
					dplyr::mutate(Adj.P.Value = p.adjust(P.Value, method = "BH", n = length(P.Value))) %>%
					dplyr::mutate(test = merged_name) %>%
					dplyr::mutate(UniqueID = !!sym(genelabel)) %>%
					dplyr::mutate(Protein.ID = !!sym(genelabel)) %>%
					dplyr::relocate(UniqueID, .before = !!sym(genelabel))

					datasource_res <- dplyr::full_join(res, res_wide, by = genelabel) #12122023 by bgao show data source

					ProteinGeneName <- dplyr::bind_rows(ProteinGeneName_merge)  %>% dplyr::distinct()

					if ("Merged Dataset" %in% names(DataInSets)){
						if 	(!(merged_name %in% DataInSets[["Merged Dataset"]]$tests_order)){
							results_long  <- DataInSets[["Merged Dataset"]]$results_long %>%
							#dplyr::filter(test != merged_name) %>%
							dplyr::bind_rows(., res)
							DataInSets[["Merged Dataset"]]$results_long = results_long
							DataInSets[["Merged Dataset"]]$tests_order = c(DataInSets[["Merged Dataset"]]$tests_order, merged_name)
						} else {
							results_long  <- DataInSets[["Merged Dataset"]]$results_long %>%
							dplyr::filter(test != merged_name) %>%
							dplyr::bind_rows(., res)
							DataInSets[["Merged Dataset"]]$results_long = results_long
						}
					} else {
						DataInSets[["Merged Dataset"]]$Name = "Merged Dataset"
						DataInSets[["Merged Dataset"]]$ShortName = "Merged Dataset"
						DataInSets[["Merged Dataset"]]$ProteinGeneName = ProteinGeneName
						DataInSets[["Merged Dataset"]]$ProteinGeneNameHeader = names(ProteinGeneName)
						DataInSets[["Merged Dataset"]]$results_long = res
						DataInSets[["Merged Dataset"]]$tests_order = merged_name
						DataInSets[["Merged Dataset"]]$Species = DataInSets[[working_project()]]$Species
						#change DS_names
						DataInSets_List<-reactiveValuesToList(DataInSets)
						DataInSets_List<-DataInSets_List[!sapply(DataInSets_List, is.null)]
						DS_names(names(DataInSets_List))
	
					}
					return(list("res" = res, "datasource_res" = datasource_res))
				})
			})

			output$mergedatatable <- DT::renderDataTable(server=TRUE,{
				mergeddata <- DataMergeReactive()$res
				mergeddata[,sapply(mergeddata,is.numeric)] <- signif(mergeddata[,sapply(mergeddata,is.numeric)],3)
				DT::datatable(mergeddata,  extensions = 'Buttons',
					options = list(	dom = 'lBfrtip', pageLength = 15,
						buttons = list(
							list(extend = "csv", text = "Download Current Page", filename = "Page_Results",	exportOptions = list(modifier = list(page = "current")))
						)
					),
				rownames= T)
			})

			output$datasourcetable <- DT::renderDataTable(server=TRUE,{
				datasource <- DataMergeReactive()$datasource_res
				datasource[,sapply(datasource,is.numeric)] <- signif(datasource[,sapply(datasource,is.numeric)],3)
				DT::datatable(datasource,  extensions = 'Buttons',
					options = list(	dom = 'lBfrtip', pageLength = 15,
						buttons = list(
							list(extend = "csv", text = "Download Current Page", filename = "Page_Results",	exportOptions = list(modifier = list(page = "current")))
						)
					),
				rownames= T)
			})

			## download big table
			output$download_mergeddata_button <- shiny::downloadHandler(
				filename = function() {
					paste("Results-", Sys.Date(), ".csv", sep="")
				},
				content = function(file) {
					mergeddata <- DataMergeReactive()$res
					mergeddata[,sapply(mergeddata,is.numeric)] <- signif(mergeddata[,sapply(mergeddata,is.numeric)],3)
					write.csv(mergeddata, file)
				}
			)

			output$download_datasource_button <- shiny::downloadHandler(
				filename = function() {
					paste("Data-", Sys.Date(), ".csv", sep="")
				},
				content = function(file) {
					write.csv(DataMergeReactive()$datasource_res, file)
				}
			)
			
		}
	)
}

