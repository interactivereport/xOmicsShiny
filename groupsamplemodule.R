###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: groupsample.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 08/23/2023
##@version 3.0
###########################################################################################################

# Need "MetaData_long","groups", "group_order", "samples", "sample_order", "tests", "tests_order"

groupsample_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
		column(3,
		),
		column(6,
			wellPanel(
				uiOutput(ns('loadedprojects')),
				tags$hr(style="border-color: black;"),
				uiOutput(ns("ui_sampletype")),
				uiOutput(ns("ui_source")),
				uiOutput(ns("ui_dest")),
				uiOutput(ns("reset")),
				tags$hr(style="border-color: black;"),
				uiOutput(ns("ui_source_s")),
				uiOutput(ns("ui_dest_s")),
				tags$hr(style="border-color: black;"),
				uiOutput(ns("ui_source_test")),
				uiOutput(ns("ui_dest_test")),
				uiOutput(ns("reset_test")),
				tags$hr(style="border-color: black;"),
				span(textOutput(ns("selectGroupSample")), style = "color:red; font-size:20px; font-family:arial; font-style:italic")
			)
		),
		column(3,
		)
	)
}

groupsample_server <- function(id) {
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

			output$ui_sampletype <- renderUI({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$MetaData_long)

				types <- DataInSets[[working_project()]]$MetaData_long %>%
				dplyr::pull(type) %>% unique()
				sampletype <- input$sampletype
				radioButtons(ns("sampletype"), label="Select Sample Type",  inline = TRUE, choices=types, selected=sampletype)
			})


			observeEvent(input$sampletype, {
				req(length(working_project()) > 0)
				sampletype <- input$sampletype
				groups <- DataInSets[[working_project()]]$MetaData_long %>%
				dplyr::filter(type == sampletype) %>%
				dplyr::pull(group) %>% unique() %>% as.character()

				DataInSets[[working_project()]]$groups <- groups
				DataInSets[[working_project()]]$group_order <- groups
				DataInSets[[working_project()]]$sample_order <- DataInSets[[working_project()]]$samples
			})

			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$groups)
				itemsall <- DataInSets[[working_project()]]$groups
				items <- DataInSets[[working_project()]]$group_order
				itemsremove <- itemsall[!(itemsall %in% items)]

				output$ui_source <- renderUI({
					shinyjqui::orderInput(ns('source'), 'Available Groups (Drag to Order):', items = items,  width = '100%',  item_class = 'primary', connect = ns('dest'))
				})
				output$ui_dest <- renderUI({
					shinyjqui::orderInput(ns('dest'), 'Drag to Remove Groups:', items = itemsremove,  width = '100%', placeholder = 'Drag items here...', item_class = 'primary', connect = ns('source'))
				})
				output$reset <- renderUI({
					actionButton(ns("reset"),"reset")
				})

				output$update <- renderUI({
					actionButton(ns("update"), "update")
				})
			})

			observeEvent(input$source_order, {
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$MetaData_long)

				MetaData = DataInSets[[working_project()]]$MetaData_long
				groups = input$source_order %>% pull(text)
				group_order <- DataInSets[[working_project()]]$group_order
				sampletype <- input$sampletype

				if (!identical(groups, group_order)) {
					DataInSets[[working_project()]]$group_order <- groups

					keep_samples <- MetaData %>%
					dplyr::filter(type == sampletype) %>%
					dplyr::filter(group %in% groups) %>%
					dplyr::arrange(factor(group, levels = groups)) %>%
					dplyr::pull(sampleid) %>% unique() %>% as.character()

					if(rlang::is_empty(keep_samples))
					keep_samples <- "No Sample Selected"
					DataInSets[[working_project()]]$sample_order <- keep_samples
				}
			})

			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$sample_order)

				itemsall <- DataInSets[[working_project()]]$samples
				items <- DataInSets[[working_project()]]$sample_order
				itemsremove <- itemsall[!(itemsall %in% items)]

				output$ui_source_s <- renderUI({
					shinyjqui::orderInput(ns('source_s'), 'Available Samples:', items = items,  width = '100%', item_class = 'success', connect = ns('dest_s'))
				})

				output$ui_dest_s <- renderUI({
					shinyjqui::orderInput(ns('dest_s'), 'Drag to Remove Samples:', items = itemsremove,  width = '100%', placeholder = 'Drag items here...', item_class = 'success', connect = ns('source_s'))
				})
			})

			observeEvent(input$source_s_order, {
				selecteditems <- input$source_s_order %>% pull(text)
				DataInSets[[working_project()]]$sample_order <- selecteditems
			})


			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$tests)

				itemsall <- DataInSets[[working_project()]]$tests
				items <- DataInSets[[working_project()]]$tests_order
				itemsremove <- itemsall[!(itemsall %in% items)]

				output$ui_source_test <- renderUI({
					shinyjqui::orderInput(ns('source_test'), 'Available Comprisons', items = items,  width = '100%', item_class = 'success', connect = ns('dest_test'))
				})

				output$ui_dest_test <- renderUI({
					shinyjqui::orderInput(ns('dest_test'), 'Drag to Remove Comprisons:', items = itemsremove,  width = '100%', placeholder = 'Drag items here...', item_class = 'success', connect = ns('source_test'))
				})
				output$reset_test <- renderUI({
					actionButton(ns("btn_test"),"reset test")
				})
			})

			observeEvent(input$source_test_order, {
				req(length(working_project()) > 0)
				selecteditems <- input$source_test_order %>% pull(text)
				DataInSets[[working_project()]]$tests_order <- selecteditems
			})

			observeEvent(input$reset, {
				req(length(working_project()) > 0)
				itemsall <- DataInSets[[working_project()]]$groups

				output$ui_source <- renderUI({
					shinyjqui::orderInput(ns('source'), 'Available Groups (Drag to Order):', items = itemsall,  width = '100%', item_class = 'primary', connect = ns('dest'))
				})
				output$ui_dest <- renderUI({
					shinyjqui::orderInput(ns('dest'), 'Drag to Remove Groups:', items = NULL,  width = '100%', placeholder = 'Drag items here...', item_class = 'primary', connect = ns('source'))
				})
			})

			observeEvent(input$btn_test,{
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$tests)
				itemsall <- DataInSets[[working_project()]]$tests
				items <- DataInSets[[working_project()]]$tests_order
				itemsremove <- itemsall[!(itemsall %in% items)]

				output$ui_source_test <- renderUI({
					shinyjqui::orderInput(ns('source_test'), 'Available Comprisons', items = itemsall,  width = '100%', item_class = 'success', connect = ns('dest_test'))
				})

				output$ui_dest_test <- renderUI({
					shinyjqui::orderInput(ns('dest_test'), 'Selected Comprisons:', items = NULL,  width = '100%', placeholder = 'Drag items here...', item_class = 'success', connect = ns('source_test'))
				})
			})

			output$selectGroupSample <- renderText({
				paste("Selected ",length(DataInSets[[working_project()]]$group_order), " out of ", length(DataInSets[[working_project()]]$groups), " Groups, ",	length(DataInSets[[working_project()]]$sample_order), " out of ", length(DataInSets[[working_project()]]$samples), " Samples.", sep="")
			})
		}
	)
}
