
###########################################################################################################
##Omics Dose-Response/Time-Course Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: monotonicmodule.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 04/23/2024
##@version 3.0
###########################################################################################################

##########################################################################################################
## Monotonic Trend
##########################################################################################################
#pkgs:  "ggpmisc,"  "parallel", "DT", "dplyr", "purrr", "plyr"

library(ggpmisc)
library(parallel)

#corenumber  = parallel::detectCores(logical = FALSE) - 2
corenumber  = 20 #fixed based on cluster node
monotonic_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
		column(3,
			wellPanel(
				uiOutput(ns('loadedprojects')),
				conditionalPanel(ns = ns,"input.upper_tabset=='Fitting Curve'",
					selectizeInput(ns("sel_gene3"),	label="Gene Name",	choices = NULL,	multiple=FALSE, options = list(placeholder =	'Type to search')),
					radioButtons(ns("separateplot3"), label="Plot Separately ", inline = TRUE, choices = c("Yes" = "Yes","No" = "No"))
				),
				conditionalPanel(ns = ns,"input.upper_tabset=='Browsing'",
					radioButtons(ns("subset"), label="Genes Used in Plot", choices=c("Browsing", "Upload Genes"), inline = TRUE, selected="Browsing"),
					conditionalPanel(ns = ns, "input.subset=='Upload Genes'", textAreaInput(ns("uploadlist"), "Enter Gene List", "", cols = 5, rows=6)),
					radioButtons(ns("sel_geneid"), label="Select Gene Label", inline = TRUE, choices=c("UniqueID", "Gene.Name","Protein.ID"), selected="UniqueID"),
					fluidRow(
						column(width=6, radioButtons(ns("psel"), label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"), inline = TRUE)),
						column(width=6, radioButtons(ns("updown"), label= "All, Up or Down?", choices= c("All"="All","Up"="Up","Down"="Down"), inline = TRUE))
					),
					fluidRow(
						column(width=6, numericInput(ns("fccut"), label= "Fold Change Threshold", value = 1.2, min=1, step=0.1)),
						column(width=6, numericInput(ns("pvalcut"), label= "P-value Threshold", value=0.01, min=0, step=0.01))
					),
					radioButtons(ns("orderby"), label="Sort by", inline = TRUE, choices =  c("p_value" = "p_value", "padjust" = "padjust","r.squared" = "r.squared", "alphabetically" = "alphabetically")),
					fluidRow(
						column(width=4,sliderInput(ns("plot_ncol"), label= "Column Number", min = 1, max = 6, step = 1, value = 3)),
						column(width=4,sliderInput(ns("plot_nrow"), label= "Row Number", min = 1, max = 9, step = 1, value = 3)),
						column(width=4,selectInput(ns("sel_page"), label="Select Page",	choices = NULL,	selected=1))
					)
				),
				conditionalPanel(ns = ns,"input.upper_tabset=='Fitting Curve'  || input.upper_tabset=='Browsing'",
					selectizeInput(ns("sel_treatment"),	label="Treatment",	choices = NULL,	multiple=TRUE),
					fluidRow(
						column(width=6,textInput(ns("xlabel"), "xlabel", value = "Conc")),
						column(width=6,textInput(ns("ylabel"), "ylabel", value = "Response"))
					),
					fluidRow(
						column(width=6,sliderInput(ns("labelfontsize"), "Label Font Size:", min = 10, max = 24, step = 2, value = 14)),
						column(width=6,sliderInput(ns("basefontsize"), "Font Size:", min = 10, max = 24, step = 2, value = 14))
					),
					radioButtons(ns("dotline"), label="Regression or Connecting Dot", choices=c("Regression", "Connecting Dot"), inline = TRUE, selected="Regression"),
					radioButtons(ns("IndividualPoint"), label="Individual Point?", inline = TRUE, choices = c("YES" = "YES","NO" = "NO"), selected = "YES"),
					radioButtons(ns("ShowErrorBar"), label="Error Bar?", inline = TRUE, choices = c("SD" = "SD","SEM" = "SEM","NO" = "NO"), selected = "SD")
				),
				conditionalPanel(ns = ns,"input.upper_tabset=='Result Table (all)'",
					radioButtons(ns("parallel"), label= "parallel processing?", choices= c("yes"="yes","no"="no"),inline = TRUE),
					numericInput(ns("core"), "Core will be used:", value = 4, min = 2, max = 12),
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
		column(9,
			tabsetPanel(id=ns("upper_tabset"),
				tabPanel(title="Fitting Curve", value ="Fitting Curve",
					tabsetPanel(id=ns("FittingCurve_tabset"),
						tabPanel(title="FittingCurve", value="FittingCurve", plotOutput(ns("FittingCurve"), height=800)
						),
						tabPanel(title="Result Table", value="Result Table", DT::dataTableOutput(ns("fitresult"))
						),
						tabPanel(title="Data Table", value ="Data Table",	DT::dataTableOutput(ns("fitdata")))
					)
				),
				tabPanel(title="Result Table (all)", value ="Result Table (all)",
					shiny::downloadButton(outputId = ns("download_results_button"),  label = "Download All as CSV file"),
					DT::dataTableOutput(ns("results"))
				),
				tabPanel(title="Browsing", value ="Browsing",
					tabsetPanel(id=ns("Browsing_tabset"),
						tabPanel(title="Plot", value="Plot", plotOutput(ns("browsing"), height=1200)
						),
						tabPanel(title="Table", value="Table", DT::dataTableOutput(ns("browsing_result"))
						)
					)
				),
				tabPanel(title="Help", htmlOutput(ns('help_monotonic')))
			)
		)
	)
}

monotonic_server <- function(id) {
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

			observe({
				req(length(working_project()) > 0)
				validate(
					need(DataInSets[[working_project()]]$data_long, "Need upload data")
				)

				data_long <- DataInSets[[working_project()]]$data_long
				req("UniqueID" %in% colnames(data_long) & "conc" %in% colnames(data_long))
				if (("group" %in% colnames(data_long)) && !("treatment" %in% colnames(data_long))) {
					data_long <- data_long %>% dplyr::rename("treatment"="group") %>%
					dplyr::rename("expr" = "response")
				}

				DataIngenes <-  data_long  %>% dplyr::pull(UniqueID) %>% unique() %>% as.character()
				updateSelectizeInput(session,'sel_gene3', choices= DataIngenes, server=TRUE)

				treatment <-  data_long %>% dplyr::pull(treatment) %>% unique() %>% as.character()
				updateSelectizeInput(session,'sel_treatment', choices= treatment,  selected=treatment)
			})

			DataExpReactive3 <- reactive({
				req(length(working_project()) > 0)
				data_long <- DataInSets[[working_project()]]$data_long
				if (("group" %in% colnames(data_long)) && !("treatment" %in% colnames(data_long))) {
					data_long <- data_long %>% dplyr::rename("treatment"="group") %>%
					dplyr::rename("expr" = "response")
				}
				shiny::validate(need(input$sel_gene3 != "","Please select a gene."))
				sel_gene = input$sel_gene3
				data_tmp = dplyr::filter(data_long, UniqueID == sel_gene)

				return(list("data_tmp"=data_tmp))
			})

			fitting_out3 <- reactive({
				req(length(working_project()) > 0)
				sel_gene <- input$sel_gene3
				sel_treatment <- input$sel_treatment
				npcx <- input$npcx3
				npcy <- input$npcy3
				labelfontsize <- input$labelfontsize
				basefontsize <- input$basefontsize
				xlabel = input$xlabel
				ylabel = input$ylabel

				DataExpReactive_tmp <- DataExpReactive3()
				data_tmp <- DataExpReactive_tmp[["data_tmp"]] %>%
				dplyr::filter(treatment %in% sel_treatment)

				dataStmp <- named_group_split(data_tmp, UniqueID, treatment)
				Res <- lapply(dataStmp, LinFitOne)
				results  <- purrr::map_df(Res, ~as.data.frame(.x), .id="id") %>%
				dplyr::filter(!is.na(direction))

				Linearformula <- y ~ x

				if (input$separateplot3 == "Yes") {
					sel_treatment = intersect(sel_treatment, unique(data_tmp[['treatment']]))
					plist <- list()
					for (onetreatment in sel_treatment) {
						dfgene1 <- data_tmp %>% as.data.frame() %>% dplyr::filter(treatment == onetreatment)

						df.summary <- dfgene1 %>%
						dplyr::group_by(conc, treatment) %>%
						dplyr::summarise(sd = sd(expr),	expr = mean(expr), n = n(), se = sd / sqrt(n), .groups = "drop") %>%
						dplyr::select(-n)

						p <- ggplot(data = dfgene1, aes(x = conc, y = expr,  group = treatment))

						if (input$ShowErrorBar == "SD") {
							p <- p + geom_errorbar(aes(ymin = expr-sd, ymax = expr+sd, color = treatment), data = df.summary, width = 0.2)
						}

						if (input$ShowErrorBar == "SEM") {
							p <- p + geom_errorbar(aes(ymin = expr-se, ymax = expr+se, color = treatment), data = df.summary, width = 0.2)
						}

						if (input$IndividualPoint == "YES") {
							p <- p + geom_jitter(aes(color = treatment), position = position_jitter(0.2))
						}

						if (input$dotline == "Connecting Dot") {
							p <- p +
							geom_line(aes(group = treatment, color = treatment), data = df.summary)
						} else { #"Regression"
							p <- p +
							geom_smooth(method = "lm", se=FALSE, formula = Linearformula) +
							stat_poly_eq(formula = Linearformula,   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),  parse = TRUE)
						}

						p <- p +
						ggtitle(onetreatment) +
						theme_bw(base_size = basefontsize) + xlab(xlabel) +  ylab(ylabel) +
						theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
							axis.text.x = element_text(angle = 0),
							legend.title = element_blank(),
							strip.text.x = element_text(size=input$labelfontsize),
						legend.position="bottom")

						plist[[onetreatment]]  <- p
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
					sel_treatment = intersect(sel_treatment, unique(data_tmp[['treatment']]))
					dfgene1 <- data_tmp %>% as.data.frame() %>%
					dplyr::filter(treatment %in% sel_treatment)

					df.summary <- dfgene1 %>%
					group_by(conc, treatment) %>%
					summarise(
						sd = sd(expr),
						expr = mean(expr), .groups = "drop"
					)

					p <- ggplot(dfgene1, aes(x=conc, y=expr, color=treatment))

					if (input$ShowErrorBar == "SD") {
						p <- p + geom_errorbar(aes(ymin = expr-sd, ymax = expr+sd, color = treatment), data = df.summary, width = 0.2)
					}

					if (input$ShowErrorBar == "SEM") {
						p <- p + geom_errorbar(aes(ymin = expr-se, ymax = expr+se, color = treatment), data = df.summary, width = 0.2)
					}

					if (input$IndividualPoint == "YES") {
						p <- p + geom_jitter(aes(color = treatment), position = position_jitter(0.2))
					}

					if (input$dotline == "Connecting Dot") {
						p <- p +
						geom_line(aes(group = treatment, color = treatment), data = df.summary)
					} else { #"Regression"
						p <- p +
						geom_smooth(method = "lm", se=FALSE, formula = Linearformula) +
						stat_poly_eq(formula = Linearformula,   aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),  parse = TRUE)
					}

					p <- p +
					theme_bw(base_size = basefontsize) + xlab(xlabel) +  ylab(ylabel) +
					theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), axis.text.x = element_text(angle = 0),legend.title = element_blank(), legend.position="bottom")

					ml <- p
				}

				return(list(plot=p, result = results))
			})

			output$FittingCurve <- renderPlot({
				fitting_out3()[["plot"]]
			})

			output$fitresult <- DT::renderDataTable({
				fitresult <- fitting_out3()[["result"]]
				fitresult <- fitresult %>%
				mutate_if(is.numeric, round, digits = 4)
				DT::datatable(fitresult, options = list(pageLength = 15))
			})

			output$fitdata <- DT::renderDataTable({
				data_tmp <- DataExpReactive3()[["data_tmp"]]
				data_tmp <- data_tmp %>%
				mutate_if(is.numeric, round, digits = 2)
				DT::datatable(data_tmp, options = list(pageLength = 15))
			})

			###########################################################################################################
			#monotonic fit all data
			observe({
				req(length(working_project()) > 0)
				updateNumericInput(session, 'core', label = "Core will be used:", value = corenumber,  min = 2, max = corenumber, step =1)
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
					req("UniqueID" %in% colnames(data_long) & "conc" %in% colnames(data_long))
					if (("group" %in% colnames(data_long)) && !("treatment" %in% colnames(data_long))) {
						data_long <- data_long %>% dplyr::rename("treatment"="group") %>%
						dplyr::rename("expr" = "response")
					}
					UniqueIDnum <- length(unique(data_long$UniqueID))
				}
				output$filteredgene3 =	 renderText({paste("<font color=\'red\'><b>Total Genes: ", UniqueIDnum, "</b></font>",sep="")})
				if (UniqueIDnum  < 500)
				updateRadioButtons(session, "parallel",  label= "parallel processing (not recommended)",  choices= c("yes"="yes","no"="no"),  selected = "no")


			})

			LinFitOne <- function(x) {
				datatmp <- x %>% dplyr::arrange(conc)
				x <- datatmp$conc
				y <-  datatmp$expr

				if (length(unique(x)) > 3) {
					lmfitres <- lm(y~x)
					lmfitsummary <- summary(lmfitres)
					b <- lmfitsummary[['coefficients']][2,'Estimate']
					direction <- if (b <= 0) "d" else "u"
					angle <- atan(b) * (180 / pi)
					r.squared <- lmfitsummary[['r.squared']]
					p_value <- lmfitsummary[['coefficients']][2,'Pr(>|t|)']
					coedf <- data.frame(direction = direction, slope = b, angle = angle, p_value = p_value,  r.squared = r.squared)
				} else {
					coedf <- data.frame(direction = NA, slope = NA, angle = NA, p_value = NA,  r.squared = NA)
				}
				coedf <- coedf  %>%
				dplyr::mutate(UniqueID = unique(datatmp$UniqueID)) %>%
				dplyr::mutate(treatment = unique(datatmp$treatment)) %>%
				dplyr::relocate(c(UniqueID,treatment), .before = direction)

				return(coedf)
			}

			MultiCoreLinFit <- function(dataS_ChunkByCore, corenum) {
				cl <- makeCluster(corenum)
				clusterEvalQ(cl, {
					library(dplyr)
				})

				clusterExport(cl=cl, varlist=c("dataS_ChunkByCore","LinFitOne"), envir=environment())

				out <- parLapply(cl, 1:length(dataS_ChunkByCore), function (k) {
					dataStmp <- dataS_ChunkByCore[[k]]
					Res <- lapply(dataStmp, LinFitOne)
					coedf <- purrr::map_df(Res, ~as.data.frame(.x), .id="id")
					return(coedf)
				})
				stopCluster(cl)
				return (out)
			}

			Monotone_all <- eventReactive(input$fitallButton3, {
				pcutoff <- input$pvalcut3
				datapoint <- input$datapoint3
				psel <- input$psel3
				parallel <- input$parallel
				corenum <- input$core

				data_long <- DataInSets[[working_project()]]$data_long
				if (("group" %in% colnames(data_long)) && !("treatment" %in% colnames(data_long))) {
					data_long <- data_long %>% dplyr::rename("treatment"="group") %>%
					dplyr::rename("expr" = "response")
				}
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

						data_long_sub <- dplyr::semi_join(x=data_long, y=filterdata, by=c("UniqueID","treatment"))
						dataS <- named_group_split(data_long_sub, UniqueID, treatment)
					} else {
						dataS <- named_group_split(data_long,UniqueID,treatment)
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

			output$results <- DT::renderDataTable({
				results <- as.data.frame("No Fitting Results")
				if (input$fitallButton3[1] == 0) {
					if (!is.null(DataInSets[[working_project()]]$results_lin))
					results <- DataInSets[[working_project()]]$results_lin
				} else {
					withProgress(message = 'Caculating...',  detail = 'This may take a while...',  {
						results <-  Monotone_all()
						DataInSets[[working_project()]]$results_lin <-  results

						if (input$saveproject3 == 1) {
							shiny::validate(need(input$projectname3!= "","Please provide project name."))
							filename <- paste("data/",input$projectname3,".RData",sep="")
							save(data_long, results, file=filename )
						}
					})
				}
				results <- results %>%  mutate_if(is.numeric, round, digits = 6)
				DT::datatable(results, options = list(pageLength = 15), rownames= FALSE, filter = 'top')
			})

			###########################################################################################################
			observe({
				req(length(working_project()) > 0)
				updateSelectInput(session,'sel_page', choices= seq_len(100))
			})

			browsing_out <- reactive({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$data_long)
				validate(
					need(DataInSets[[working_project()]]$results_lin, "Need fitting results")
				)
				validate(
					need(DataInSets[[working_project()]]$data_long, "Need data")
				)

				results_lin <- DataInSets[[working_project()]]$results_lin
				data_long <- DataInSets[[working_project()]]$data_long
				if (("group" %in% colnames(data_long)) && !("treatment" %in% colnames(data_long))) {
					data_long <- data_long %>% dplyr::rename("treatment"="group") %>%
					dplyr::rename("expr" = "response")
				}
				labelfontsize <- as.numeric(input$labelfontsize)
				basefontsize <- as.numeric(input$basefontsize)
				xlabel <- input$xlabel
				ylabel <- input$ylabel
				genelabel <- input$sel_geneid
				sel_treatment <- input$sel_treatment

				orderby <- input$orderby

				ncol = as.numeric(input$plot_ncol)
				nrow = as.numeric(input$plot_nrow)
				sel_page <- as.numeric(input$sel_page)-1
				numberpage = ncol * nrow
				startslice = sel_page * numberpage  + 1
				endslice = startslice + numberpage -1


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

					if (input$subset == "Upload Genes") {
						req(input$uploadlist)
						gene_list <- input$uploadlist
						gene_list <- ProcessUploadGeneList(gene_list)
						validate(need(length(gene_list) > 0, message = "Please input at least 1 valid gene."))
						if (!is.null(DataInSets[[working_project()]]$ProteinGeneName)) {
							ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
							ProteinGeneName_sel <- dplyr::filter(ProteinGeneName, (UniqueID %in% gene_list) | (Protein.ID %in% gene_list) | (toupper(Gene.Name) %in% toupper(gene_list)))
							validate(need(nrow(ProteinGeneName_sel) > 0, message = "Please input at least 1 matched gene."))
							gene_list <- ProteinGeneName_sel %>% dplyr::pull(UniqueID)

						}
						results_lin <- results_lin %>%
						dplyr::filter(UniqueID %in% gene_list)
					}

					sel_gene <- results_lin %>%
					dplyr::group_by(UniqueID) %>%
					dplyr::slice(which.min(!!as.symbol(field))) %>%
					dplyr::ungroup() %>%
					dplyr::arrange(!!as.symbol(field)) %>%
					dplyr::slice(startslice:endslice) %>%
					dplyr::pull(UniqueID)
				} else {
					if (input$subset == "Upload Genes") {
						req(input$uploadlist)
						gene_list <- input$uploadlist
						gene_list <- ProcessUploadGeneList(gene_list)
						validate(need(length(gene_list) > 0, message = "Please input at least 1 valid gene."))
						data_long <- data_long %>%
						dplyr::filter(UniqueID %in% gene_list)
					}

					sel_gene <- data_long %>%
					dplyr::distinct(UniqueID) %>%
					dplyr::slice(startslice:endslice) %>%
					dplyr::pull(UniqueID)
				}

				data_long_tmp  <- dplyr::filter(data_long, UniqueID %in% sel_gene) %>%
				dplyr::filter(treatment %in% sel_treatment) %>% as.data.frame() %>%
				dplyr::arrange(Gene.Name)
				data_long_tmp$labelgeneid = data_long_tmp[,match(genelabel,colnames(data_long_tmp))]
				data_long_tmp$treatment = factor(data_long_tmp$treatment, levels = sel_treatment)

				df.summary <- data_long_tmp %>%
				group_by(labelgeneid, conc, treatment) %>%
				dplyr::summarise(sd = sd(expr), expr = mean(expr),  n = n(), se = sd / sqrt(n), .groups = "drop") %>%
				dplyr::select(-n)

				dataStmp <- named_group_split(data_long_tmp, labelgeneid, treatment)

				Res <- lapply(dataStmp, LinFitOne)
				results  <- purrr::map_df(Res, ~as.data.frame(.x), .id="id") %>%
				dplyr::filter(!is.na(direction))

				Linearformula <- y ~ x
				colorpal <- UserColorPlalette(colpalette = "Dark2", items = unique(data_long_tmp$treatment))

				p <- ggplot(data_long_tmp, aes(x=conc, y=expr, color=treatment)) +
				facet_wrap(~labelgeneid, scales = "free", nrow = nrow, ncol = ncol)

				if (input$ShowErrorBar == "SD") {
					p <- p + geom_errorbar(aes(ymin = expr-sd, ymax = expr+sd, color = treatment), data = df.summary, width = 0.2)
				}

				if (input$ShowErrorBar == "SEM") {
					p <- p + geom_errorbar(aes(ymin = expr-se, ymax = expr+se, color = treatment), data = df.summary, width = 0.2)
				}

				if (input$IndividualPoint == "YES") {
					p <- p + geom_jitter(aes(color = treatment), position = position_jitter(0.2))
				}

				if (input$dotline == "Connecting Dot") {
					p <- p +
					geom_line(aes(group = treatment, color = treatment), data = df.summary)
				} else { #"Regression"
					p <- p +
					geom_smooth(method = "lm", se=FALSE, formula = Linearformula) +
					stat_poly_eq(formula = Linearformula, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),  parse = TRUE)
				}

				p <- p +
				scale_fill_manual(values = colorpal) +
				scale_colour_manual(values = colorpal) +
				theme_bw(base_size = basefontsize) + xlab(xlabel) +  ylab(ylabel) +
				theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), axis.text.x = element_text(angle = 0), legend.title = element_blank(), plot.title = element_text(size=labelfontsize), legend.position="bottom")

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

			output$download_results_button <- shiny::downloadHandler(
				filename = function() {
					paste("Results-", Sys.Date(), ".csv", sep="")
				},
				content = function(file) {
					fitresult <- browsing_out()[["result"]]
					fitresult <- fitresult %>%
					mutate_if(is.numeric, round, digits = 4)
					write.csv(results, file)
				}
			)
		}
	)
}

