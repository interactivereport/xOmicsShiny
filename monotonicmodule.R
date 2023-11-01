
###########################################################################################################
##Omics Dose-Response/Time-Course Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: monotonicmodule.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 05/23/2023
##@version 3.0
###########################################################################################################

##########################################################################################################
## Monotonic Trend
##########################################################################################################
#pkgs:  "ggpmisc,"  "parallel", "DT", "dplyr", "purrr", "plyr"  

library(ggpmisc)
library(parallel)

corenumber  = parallel::detectCores(logical = FALSE)

monotonic_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
		column(2,
			wellPanel(
				conditionalPanel(ns = ns,"input.expression_tabset3=='Fitting Curve'",
					selectizeInput(ns("sel_gene3"),	label="Gene Name",	choices = NULL,	multiple=FALSE, options = list(placeholder =	'Type to search')),
					selectizeInput(ns("sel_treatment3"),	label="Treatment",	choices = NULL,	multiple=TRUE),
					radioButtons(ns("separateplot3"), label="Plot Separately ", inline = TRUE, choices = c("Yes" = "Yes","No" = "No"))
				),
				conditionalPanel(ns = ns,"input.expression_tabset3=='Browsing'",
					selectizeInput(ns("sel_treatment4"),	label="Treatment",	choices = NULL,	multiple=TRUE),
					radioButtons(ns("orderby3"), label="Sort by", inline = TRUE, choices =  c("p_value" = "p_value", "padjust" = "padjust","r.squared" = "r.squared")),
					column(width=6,selectInput(ns("sel_page3"),	label="Select Page",	choices = NULL,	selected=1)),
					column(width=6,selectInput(ns("numperpage3"), label= "Plot Number per Page", choices= c("4"=4,"6"=6,"9"=9), selected=6))
				),
				conditionalPanel(ns = ns,"input.expression_tabset3=='Fitting Curve'  || input.expression_tabset3=='Browsing'",
					textInput(ns("xlabel3"), "xlabel", value = "Conc"),
					textInput(ns("ylabel3"), "ylabel", value = "Response"),
					sliderInput(ns("basefontsize3"), "Font Size:", min = 10, max = 24, step = 2, value = 14)
				),
				conditionalPanel(ns = ns,"input.expression_tabset3=='Result Table (Isotone/Monotone)'",
					radioButtons(ns("parallel3"), label= "parallel processing?", choices= c("yes"="yes","no"="no"),inline = TRUE),
					numericInput(ns("core3"), "Core will be used:", value = 4, min = 2, max = 12),
					radioButtons(ns("psel3"), label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
					numericInput(ns("pvalcut3"), label= "Choose P-value Threshold",  value=0.01, min=0, step=0.001),
					numericInput(ns("datapoint3"), label= "Minimal Data Points",  value=4, min=4, step=1),
					uiOutput(ns("filteredgene3")),
					radioButtons(ns("saveproject3"), label="Save Result?", inline = TRUE,  choices=c("Yes"= 1, "No"= 0), selected=0),
					textInput(ns("projectname3"), "Project Name"),
					actionButton(ns("fitallButton3"), "Fit filtered data")
				)
			)
		),
		column(10,
			tabsetPanel(id=ns("expression_tabset3"),
				tabPanel(title="Fitting Curve",
					#actionButton("FittingCurve3", "Save to output"),
					fluidRow(
						plotOutput(ns("FittingCurve3"), height=800)
					),
					fluidRow(
						DT::dataTableOutput(ns("fitresult3"))
					)
				),
				tabPanel(title="Data Table",	DT::dataTableOutput(ns("fitdata3"))),
				tabPanel(title="Result Table (Isotone/Monotone)",
					#actionButton("results3", "Save to output"),
					#shinycssloaders::withSpinner(dataTableOutput(ns("results3")),type = 4, size = 2,color = "#0000FF")
					dataTableOutput(ns("results3"))
				),
				tabPanel(title="Browsing",
					fluidRow(
						plotOutput(ns("browsing"), height=800)
					),
					fluidRow(
						DT::dataTableOutput(ns("browsing_result"))
					)
				),
				tabPanel(title="Help", htmlOutput(ns('help_FittingCurveDRC3')))
			)
		)
	)
}

monotonic_server <- function(id) {
	shiny::moduleServer(id,
		function(input, output, session) {
			ns <- shiny::NS(id)
			observe({
				req(length(working_project()) > 0)
				validate(
					need(DataInSets[[working_project()]]$data_long, "Need upload data")
				)

				data_long <- DataInSets[[working_project()]]$data_long
				req("UniqueID" %in% colnames(data_long) & "group" %in% colnames(data_long))

				DataIngenes <-  data_long  %>% dplyr::pull(UniqueID) %>% unique() %>% as.character()
				updateSelectizeInput(session,'sel_gene3', choices= DataIngenes, server=TRUE)

				group <-  data_long %>% dplyr::pull(group) %>% unique() %>% as.character()
				updateSelectizeInput(session,'sel_treatment3', choices= group,  selected=group)
				updateSelectizeInput(session,'sel_treatment4', choices= group,  selected=group)
			})

			DataExpReactive3 <- reactive({
				req(length(working_project()) > 0)
				data_long <- DataInSets[[working_project()]]$data_long

				shiny::validate(need(input$sel_gene3 != "","Please select a gene."))
				sel_gene = input$sel_gene3
				data_tmp = dplyr::filter(data_long, UniqueID == sel_gene)

				return(list("data_tmp"=data_tmp))
			})

			fitting_out3 <- reactive({
				req(length(working_project()) > 0)
				sel_gene <- input$sel_gene3
				sel_treatment <- input$sel_treatment3
				npcx <- input$npcx3
				npcy <- input$npcy3
				labelfontsize <- input$labelfontsize3
				basefontsize <- input$basefontsize3
				xlabel = input$xlabel3
				ylabel = input$ylabel3

				DataExpReactive_tmp <- DataExpReactive3()
				data_tmp <- DataExpReactive_tmp[["data_tmp"]] %>%
				dplyr::filter(group %in% sel_treatment)

				dataStmp <- named_group_split(data_tmp,UniqueID, group)
				Res <- lapply(dataStmp, LinFitOne)
				results  <- purrr::map_df(Res, ~as.data.frame(.x), .id="id") %>%
				dplyr::filter(!is.na(direction))

				Linearformula <- y ~ x

				if (input$separateplot3 == "Yes") {
					sel_treatment = intersect(sel_treatment, unique(data_tmp[['group']]))
					plist <- list()
					for (onegroup in sel_treatment) {
						dfgene1 <- data_tmp %>% as.data.frame() %>% dplyr::filter(group == onegroup)

						p <- ggplot(data = dfgene1, aes(x = conc, y = response)) +
						geom_smooth(method = "lm", se=FALSE, color="black", formula = Linearformula) +
						stat_poly_eq(formula = Linearformula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE) +
						geom_point() +
						ggtitle(onegroup) +
						theme_bw(base_size = basefontsize) + xlab(xlabel) +  ylab(ylabel) +
						theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), axis.text.x = element_text(angle = 0),legend.title = element_blank(), legend.position="bottom")

						plist[[onegroup]]  <- p
					}

					if (length(plist) == 1)
					ml <- plist[[1]]
					if (length(plist) == 2)
					ml <- marrangeGrob(plist, nrow=1, ncol=2, top = sel_gene)
					if (length(plist) > 2 & length(plist) < 5)
					ml <- marrangeGrob(plist, nrow=2, ncol=2, top = sel_gene)
					if (length(plist) > 4) {
						nrow = ceiling(length(plist)/3)
						ml <- marrangeGrob(plist, nrow=nrow, ncol=3, top = sel_gene)
					}
				} else {
					sel_treatment = intersect(sel_treatment, unique(data_tmp[['group']]))
					dfgene1 <- data_tmp %>% as.data.frame() %>%
					dplyr::filter(group %in% sel_treatment)

					p <- ggplot(dfgene1, aes(x=conc, y=response, color=group)) +
					geom_smooth(method = "lm", se=FALSE, formula = Linearformula) +
					stat_poly_eq(formula = Linearformula,   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),  parse = TRUE) +
					geom_point() +
					theme_bw(base_size = basefontsize) + xlab(xlabel) +  ylab(ylabel) +
					theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), axis.text.x = element_text(angle = 0),legend.title = element_blank(), legend.position="bottom")

					ml <- p
				}

				return(list(plot=ml, result = results))
			})

			output$FittingCurve3 <- renderPlot({
				fitting_out3()[["plot"]]
			})

			output$fitresult3 <- DT::renderDataTable({
				fitresult <- fitting_out3()[["result"]]
				fitresult <- fitresult %>%
				mutate_if(is.numeric, round, digits = 4)
				DT::datatable(fitresult, options = list(pageLength = 15))
			})

			output$fitdata3 <- DT::renderDataTable({
				data_tmp <- DataExpReactive3()[["data_tmp"]]
				data_tmp <- data_tmp %>%
				mutate_if(is.numeric, round, digits = 2)
				DT::datatable(data_tmp, options = list(pageLength = 15))
			})

			###########################################################################################################
			#monotonic fit all data
			observe({
				req(length(working_project()) > 0)
				updateNumericInput(session, 'core3', label = "Core will be used:", value = corenumber,  min = 2, max = corenumber, step =1)
				pcutoff <- input$pvalcut3
				datapoint <- input$datapoint3
				psel <- input$psel3

				if (!is.null(psel) & !is.null(DataInSets[[working_project()]]$statresult)) {
					if (psel == "Pval") {
						tmpdat3 = DataInSets[[working_project()]]$statresult %>% dplyr::filter(pvalue < pcutoff & n >= datapoint)
					}
					if (psel == "Padj") {
						tmpdat3 =DataInSets[[working_project()]]$statresult %>% dplyr::filter(padjust < pcutoff & n >= datapoint)
					}
					UniqueIDnum <- nrow(tmpdat3)
				} else {
					data_long <- DataInSets[[working_project()]]$data_long
					UniqueIDnum <- length(unique(data_long$UniqueID))
				}
				output$filteredgene3 =	 renderText({paste("<font color=\'red\'><b>Total Genes: ", UniqueIDnum, "</b></font>",sep="")})
				if (UniqueIDnum  < 500)
				updateRadioButtons(session, "parallel3",  label= "parallel processing (not recommended)",  choices= c("yes"="yes","no"="no"),  selected = "no")


			})

			LinFitOne <- function(x) {
				datatmp <- x %>% dplyr::arrange(conc)
				x <- datatmp$conc
				y <-  datatmp$response

				if (length(unique(x)) > 3) {
					lmfitres <- lm(y~x)
					lmfitsummary <- summary(lmfitres)
					b <- lmfitsummary[['coefficients']][2,'Estimate']
					direction <- if (b <= 0) "d" else "u"
					angle <- atan(b) * (180 / pi)
					#adj.r.squared <- lmfitsummary[['adj.r.squared']]
					r.squared <- lmfitsummary[['r.squared']]
					p_value <- lmfitsummary[['coefficients']][2,'Pr(>|t|)']
					coedf <- data.frame(direction = direction, slope = b, angle = angle, p_value = p_value,  r.squared = r.squared)
				} else {
					coedf <- data.frame(direction = NA, slope = NA, angle = NA, p_value = NA,  r.squared = NA)
				}
				coedf <- coedf  %>%
				dplyr::mutate(UniqueID = unique(datatmp$UniqueID)) %>%
				dplyr::mutate(group = unique(datatmp$group)) %>%
				dplyr::relocate(c(UniqueID,group), .before = direction)

				return(coedf)
			}

			MultiCoreLinFit <- function(dataS_ChunkByCore, corenum) {
				cl <- makeCluster(corenum)
				clusterEvalQ(cl, {
					library(dplyr)
				}
			)

			clusterExport(cl=cl, varlist=c("dataS_ChunkByCore","LinFitOne"), envir=environment())

			out <- parLapply(cl, 1:length(dataS_ChunkByCore), function (k) {
				dataStmp <- dataS_ChunkByCore[[k]]
				Res <- lapply(dataStmp, LinFitOne)
				coedf <- purrr::map_df(Res, ~as.data.frame(.x), .id="id")
				return(coedf)
			} )
			stopCluster(cl)
			return (out)
		}

		Monotone_all <- eventReactive(input$fitallButton3, {
			pcutoff <- input$pvalcut3
			datapoint <- input$datapoint3
			psel <- input$psel3
			parallel <- input$parallel3
			corenum <- input$core3

			data_long <- DataInSets[[working_project()]]$data_long

			if (!("UniqueID" %in% colnames(data_long))) {
				results = data.frame("no ID" ="no result")
			} else {
				if (!is.null(DataInSets[[working_project()]]$statresult)) {
					if (psel == "Pval") {
						filterdata = DataInSets[[working_project()]]$statresult %>% dplyr::filter(pvalue < pcutoff & n >= datapoint)
					}
					if (psel == "Padj") {
						filterdata =DataInSets[[working_project()]]$statresult %>% dplyr::filter(padjust < pcutoff & n >= datapoint)
					}

					data_long_sub <- dplyr::semi_join(x=data_long, y=filterdata, by=c("UniqueID","group"))
					dataS <- named_group_split(data_long_sub, UniqueID, group)
				} else {
					dataS <- named_group_split(data_long,UniqueID,group)
				}
				if (parallel == "yes")  {
					dataS_ChunkByCore <-  split(dataS, cut(seq_along(dataS), corenum, labels = FALSE))
					out <- MultiCoreLinFit(dataS_ChunkByCore, corenum)
					results <- plyr::ldply(out, data.frame) %>%
					dplyr::filter(!is.na(direction))
				} else {
					Res <- lapply(dataS, LinFitOne)
					results  <- purrr::map_df(Res, ~as.data.frame(.x), .id="id") %>%
					dplyr::filter(!is.na(direction))
				}
				results$padjust <- p.adjust(results$p_value, "BH",  n = length(results$p_value))
				results <- results %>% dplyr::mutate_if(is.numeric, round, digits = 6)  %>%
				dplyr::relocate(padjust, .before = r.squared)

				return(results)
			}


			return(results)
		})

		output$results3 <- DT::renderDataTable({
			results3 <- as.data.frame("No Fitting Results")
			if (input$fitallButton3[1] == 0) {
				if (!is.null(DataInSets[[working_project()]]$results_lin))
				results3 <- DataInSets[[working_project()]]$results_lin
			} else {
				withProgress(message = 'Caculating...',  detail = 'This may take a while...',  {
					results3 <-  Monotone_all()
					DataInSets[[working_project()]]$results_lin <-  results3

					if (input$saveproject3 == 1) {
						shiny::validate(need(input$projectname3!= "","Please provide project name."))
						filename <- paste("data/",input$projectname3,".RData",sep="")
						save(data_long, results3, file=filename )
					}
				})
			}
			results3 <- results3 %>%  mutate_if(is.numeric, round, digits = 6)
			DT::datatable(results3, options = list(pageLength = 15), rownames= FALSE)
		})

		###########################################################################################################
		#browsing
		observe({
			req(length(working_project()) > 0)
			updateSelectInput(session,'sel_page3', choices= seq_len(100))
		})

		browsing_out <- reactive({
			req(length(working_project()) > 0)
			req(DataInSets[[working_project()]]$data_long)
			results_lin <- DataInSets[[working_project()]]$results_lin
			data_long <- DataInSets[[working_project()]]$data_long
			sel_treatment <- input$sel_treatment4
			labelfontsize <- input$labelfontsize3
			basefontsize <- input$basefontsize3
			xlabel <- input$xlabel3
			ylabel <- input$ylabel3
			numperpage <- as.numeric(input$numperpage3)
			sel_page <- as.numeric(input$sel_page3)-1
			orderby <- input$orderby3

			startslice = sel_page * numperpage  + 1
			endslice = startslice + numperpage -1

			if (orderby == "p_value") {
				field <- "p_value"
			}
			if (orderby == "padjust") {
				field <- "padjust"
			}
			if (orderby == "r.squared") {
				field <- "r.squared"
			}

			if (!is.null(results_lin)) {
				sel_gene <- results_lin %>%
				dplyr::group_by(UniqueID) %>%
				dplyr::slice(which.min(!!as.symbol(field))) %>%
				dplyr::ungroup() %>%
				dplyr::arrange(!!as.symbol(field)) %>%
				dplyr::slice(startslice:endslice) %>%
				dplyr::pull(UniqueID)
			} else {
				sel_gene <- data_long %>%
				dplyr::distinct(UniqueID) %>%
				dplyr::slice(startslice:endslice) %>%
				dplyr::pull(UniqueID)
			}

			data_long_tmp  <- dplyr::filter(data_long, UniqueID %in% sel_gene) %>%
			dplyr::filter(group %in% sel_treatment)

			dataStmp <- named_group_split(data_long_tmp, UniqueID, group)

			Res <- lapply(dataStmp, LinFitOne)
			results  <- purrr::map_df(Res, ~as.data.frame(.x), .id="id") %>%
			dplyr::filter(!is.na(direction))

			if(numperpage==4) {
				nrow = 2; ncol = 2
			} else if(numperpage==6) {
				nrow = 2; ncol = 3
			} else {
				nrow = 3; ncol = 3
			}
			Linearformula <- y ~ x
			p <- ggplot(data_long_tmp, aes(x=conc, y=response, color=group)) +
			facet_wrap(~ UniqueID, scales = "free", nrow = nrow, ncol = ncol) +
			geom_smooth(method = "lm", se=FALSE, formula = Linearformula) +
			stat_poly_eq(formula = Linearformula,   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),  parse = TRUE) +
			geom_point() +
			theme_bw(base_size = basefontsize) + xlab(xlabel) +  ylab(ylabel) +
			theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), axis.text.x = element_text(angle = 0),legend.title = element_blank(), legend.position="bottom")

			return(list(plot=p, result = results))
		})


		output$browsing <- renderPlot({
			browsing_out()[["plot"]]
		})


		output$browsing_result  <- DT::renderDataTable({
			fitresult <- browsing_out()[["result"]]
			fitresult <- fitresult %>%
			mutate_if(is.numeric, round, digits = 4)
			DT::datatable(fitresult, options = list(pageLength = 10))
		})
	}
)
}

