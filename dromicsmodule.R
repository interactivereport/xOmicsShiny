###########################################################################################################
##Omics Dose-Response/Time-Course Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: dromicsmodule.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 09/06/2024
##@version 3.0

## fitting functions are modified from DRomics R package (https://cran.r-project.org/web/packages/DRomics/)
## util-basicandfitfunc.R from https://cran.r-project.org/web/packages/DRomics/
###########################################################################################################

##########################################################################################################
## Curve Fitting Plot (no IC50/EC50)
##########################################################################################################
#pkgs: "ggpmisc",  "parallel", "DT", "dplyr","scales","purrr", "tidyr"
library(ggpmisc)
library(parallel)

#corenumber  = parallel::detectCores(logical = FALSE) -2
corenumber  = 20 #fixed based on cluster node
dromics_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
		column(3,
			wellPanel(
				uiOutput(ns('loadedprojects')),
				conditionalPanel(ns = ns,"input.upper_tabset=='Fitting Curve'",
					selectizeInput(ns("sel_treatment2"),	label="Treatment", choices = NULL,	multiple=FALSE),
					selectizeInput(ns("sel_gene2"),	label="Gene Name",	choices = NULL,	multiple=FALSE, options = list(placeholder =	'Type to search'))
				),
				conditionalPanel(ns = ns,"input.upper_tabset=='Browsing'",
					selectizeInput(ns("sel_treatment2b"),	label="Treatment", choices = NULL,	multiple=TRUE),
					selectizeInput(ns("typology"),	label="Typology", choices = NULL,	multiple=TRUE),
					selectizeInput(ns("trend"),	label="Trend", choices = NULL,	multiple=TRUE),
					radioButtons(ns("subset"), label="Genes Used in Plot", choices=c("Browsing", "Upload Genes"), inline = TRUE, selected="Browsing"),
					conditionalPanel(ns = ns, "input.subset=='Upload Genes'", textAreaInput(ns("uploadlist"), "Enter Gene List", "", cols = 5, rows=6)),
					radioButtons(ns("sel_geneid"), label="Select Gene Label", inline = TRUE, choices=c("UniqueID", "Gene.Name","Protein.ID"), selected="UniqueID"),
					fluidRow(
						column(width=4,sliderInput(ns("plot_ncol"), label= "Column Number", min = 1, max = 6, step = 1, value = 3)),
						column(width=4,sliderInput(ns("plot_nrow"), label= "Row Number", min = 1, max = 9, step = 1, value = 3)),
						column(width=4,selectInput(ns("sel_page"), label="Select Page",	choices = NULL,	selected=1))
					)

				),
				conditionalPanel(ns = ns,"input.upper_tabset=='Fitting Curve' || input.upper_tabset=='Browsing'",
					radioButtons(ns("npcx"), label="Label X Position", inline = TRUE, choices =  c("auto" = "auto", "left" = "left","middle"="middle","right"="right")),
					radioButtons(ns("npcy"), label="Label Y Position", inline = TRUE, choices =  c("auto" = "auto", "top" = "top","center" = "center","bottom" = "bottom")),
					radioButtons(ns("logbase"), label="Log tansform", inline = TRUE,  choices=c("Non"= 1, "10"= 10), selected=1),
					fluidRow(
						column(width=6,textInput(ns("xlabel"), "xlabel", value = "Time")),
						column(width=6,textInput(ns("ylabel"), "ylabel", value = "expr"))
					),
					fluidRow(
						column(width=6,sliderInput(ns("labelfontsize"), "Label Font Size:", min = 10, max = 24, step = 2, value = 14)),
						column(width=6,sliderInput(ns("basefontsize"), "Font Size:", min = 10, max = 24, step = 2, value = 14))
					)
				),
				conditionalPanel(ns = ns,"input.upper_tabset=='Result Table (all)'",
					radioButtons(ns("parallel"), label= "parallel processing?", choices= c("yes"="yes","no"="no"),inline = TRUE),
					numericInput(ns("core"), "Core will be used:", value = 4, min = 2, max = 12),
					radioButtons(ns("psel2"), label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
					numericInput(ns("pvalcut2"), label= "Choose P-value Threshold",  value=0.01, min=0, step=0.001),
					numericInput(ns("datapoint2"), label= "Minimal Data Points",  value=5, min=5, step=1),
					radioButtons(ns("logbaseb"), label="Log tansform", inline = TRUE,  choices=c("Non"= 1, "10"= 10), selected=10),
					uiOutput(ns("filteredgene2")),
					radioButtons(ns("saveproject2"), label="Save Result?", inline = TRUE,  choices=c("Yes"= 1, "No"= 0), selected=0),
					textInput(ns("projectname2"), "Project Name"),
					actionButton(ns("fitallButton2"), "Fit filtered data")
				)
			)
		),
		column(9,
			tabsetPanel(id=ns("upper_tabset"),
				tabPanel(title="Fitting Curve", value ="Fitting Curve",
					tabsetPanel(id=ns("FittingCurve_tabset"),
						tabPanel(title="FittingCurve", value="FittingCurve", plotOutput(ns("ModelSelection"), height=800)
						),
						tabPanel(title="Result Table", value="Result Table", DT::dataTableOutput(ns("fitresult"))
						),
						tabPanel(title="Data Table",	value="Data Table", DT::dataTableOutput(ns("fitdata")))
					)
				),
				tabPanel(title="Result Table (all)", value="Result Table (all)",
					shiny::downloadButton(outputId = ns("download_results_button"),  label = "Download All as CSV file"),
					DT::dataTableOutput(ns("results"))
				),
				tabPanel(title="Browsing", value="Browsing",
					tabsetPanel(id=ns("Browsing_tabset"),
						tabPanel(title="Plot", value="Plot", plotOutput(ns("browsing"), height=1200)
						),
						tabPanel(title="Table", value="Table", DT::dataTableOutput(ns("browsing_result"))
						)
					)
				),
				tabPanel(title="Help", htmlOutput('help_dromics'))
			)
		)
	)
}

dromics_server <- function(id) {
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
				req(DataInSets[[working_project()]]$data_long)
				data_long <- DataInSets[[working_project()]]$data_long
				req("UniqueID" %in% colnames(data_long) & "conc" %in% colnames(data_long))
				if (("group" %in% colnames(data_long)) && !("treatment" %in% colnames(data_long))) {
					data_long <- data_long %>% dplyr::rename("treatment"="group") %>%
					dplyr::rename("expr" = "response")
				}
				
				DataIngenes <-  data_long  %>% dplyr::pull(UniqueID) %>% unique() %>% as.character()
				updateSelectizeInput(session,'sel_gene2', choices= DataIngenes, server=TRUE)
				treatment <-  data_long %>% dplyr::pull(treatment) %>% unique() %>% as.character()
				updateSelectizeInput(session,'sel_treatment2', choices= treatment,  selected=treatment[1])
				updateSelectizeInput(session,'sel_treatment2b', choices= treatment,  selected=treatment)
			})

			DataExpReactive2 <- reactive({
				req(length(working_project()) > 0)
				data_long <- DataInSets[[working_project()]]$data_long
				if (("group" %in% colnames(data_long)) && !("treatment" %in% colnames(data_long))) {
					data_long <- data_long %>% dplyr::rename("treatment"="group") %>%
					dplyr::rename("expr" = "response")
				}

				shiny::validate(need(input$sel_gene2 != "","Please select a gene"))
				sel_gene = input$sel_gene2
				data_tmp = dplyr::filter(data_long, UniqueID ==sel_gene)
				return(data_tmp)
			})

			ModelSelection_out2 <- reactive({
				sel_gene <- input$sel_gene2
				sel_treatment3 <- input$sel_treatment2
				npcx <- input$npcx
				npcy <- input$npcy
				labelfontsize <- input$labelfontsize
				basefontsize <- input$basefontsize
				logbase = as.numeric(input$logbase)
				xlabel = input$xlabel
				ylabel = input$ylabel

				data_tmp <- DataExpReactive2()

				sel_treatment3 = intersect(sel_treatment3, unique(data_tmp[['treatment']]))
				dfgene1 <- data_tmp %>% dplyr::filter(treatment == sel_treatment3)

				data_tmp <- dfgene1 %>% arrange(conc)

				dose <- data_tmp$conc
				doseranks <- as.numeric(as.factor(dose))
				signal <- data_tmp$expr

				signalm <- data_tmp %>%
				group_by(conc) %>%
				dplyr::summarize(Mean = mean(expr, na.rm=TRUE)) %>%
				arrange(conc)  %>%
				pull(Mean)

				doseu <-  data_tmp %>%
				group_by(conc) %>%
				dplyr::summarize(Mean = mean(expr, na.rm=TRUE)) %>%
				arrange(conc)  %>%
				pull(conc)

				# preparation of data for modelling with nls
				dset <- data.frame(signal = signal, dose = dose, doseranks = doseranks)
				shiny::validate(need(nrow(dset) >= 5, "Need data point >= 5"))

				# calculations for starting values and other uses
				dosemin <- min(dose)
				dosemax <- max(dose)
				dosemed <- median(dose[dose!=0])

				# number of points per dose-response curve
				npts <- length(dose)
				ndoses <- length(unique(dose))
				lessthan5doses <- ndoses < 5

				# Information criterion definition
				information.criterion = "AICc"
				AICdigits <- 2 # number of digits for rounding the AIC values
				information.criterion <- match.arg(information.criterion, c("AICc", "BIC", "AIC"))

				# kcrit gives the argument k to pass to function AIC()
				# dependeing of the number of parameters of the model
				# (1 to 5, corresponding to the index of the vector)
				if (information.criterion == "AIC"){
					kcrit <- rep(2, 5)
				} else if (information.criterion == "AICc") {
					nparwithsigma <- 1:5 + 1
					kcrit <- 2*npts /(npts - nparwithsigma - 1)
				} else { # BIC last choice
					lnpts <- log(npts)
					kcrit <- rep(lnpts, 5)
				}

				equalcdG <- FALSE # use to define the value of c equal to d in the Gauss4p model if needed
				equalcdLG <- FALSE # use to define the value of c equal to d in the LGauss4p model if needed

				if (any(!complete.cases(dset))){
					# removing lines with NA values for the signal
					dset <- dset[complete.cases(dset$signal), ]
				}

				# for choice of the linear trend (decreasing or increasing)
				modlin <- lm(signal ~ doseranks, data = dset)
				increaseranks <- coef(modlin)[2] >= 0
				increaseminmax <- dset$dose[which.min(dset$signal)] < dset$dose[which.max(dset$signal)]

				# for choice of the quadratic trend (Ushape or Umbrella shape)
				modquad <- lm(signal ~ doseranks + I(doseranks^2), data = dset)
				Ushape <- coef(modquad)[3] >= 0

				#
				#dose_log_transfo = TRUE
				npts = 50

				if (logbase == 10){
					minx <- min(dose[dose != 0])
					maxx <- max(dose)
					xplot <- c(0, 10^seq(log10(minx), log10(maxx), length.out = npts))
				} else {
					xplot <- seq(0, max(dose), length.out = npts)
				}

				resultlist <- list()
				predictedlist <- list()

				################### Gauss fit ########################
				# No 1
				startGauss4p <- startvalGauss4pnls(xm = doseu, ym = signalm, Ushape = Ushape)
				Gauss4p <- suppressWarnings(try(nls(formGauss4p, start = startGauss4p, data = dset, lower = c(0, -Inf, 0, -Inf), algorithm = "port"), silent = TRUE))

				if (lessthan5doses)	{
					if (!inherits(Gauss4p, "try-error")) {
						equalcdG <- TRUE
						Gauss <- Gauss4p
						AICGaussi <- round(AIC(Gauss4p, k = kcrit[4]), digits = AICdigits)
					} else {
						AICGaussi <- Inf
					}
				} else {
					startGauss5p <- startvalGauss5pnls(xm = doseu, ym = signalm, Ushape = Ushape)
					Gauss5p <- suppressWarnings(try(nls(formGauss5p, start = startGauss5p, data = dset,	lower = c(0, -Inf, -Inf, 0, -Inf), algorithm = "port"), silent = TRUE))

					#### convergence of both models
					if ((!inherits(Gauss4p, "try-error")) & (!inherits(Gauss5p, "try-error")))	{ ## both fitting successfully
						AICGauss4p <- round(AIC(Gauss4p, k = kcrit[4]), digits = AICdigits)
						AICGauss5p <- round(AIC(Gauss5p, k = kcrit[5]), digits = AICdigits)

						if (AICGauss5p < AICGauss4p) {
							Gauss <- Gauss5p
							AICGaussi <- AICGauss5p
						} else {
							Gauss <- Gauss4p
							equalcdG <- TRUE
							AICGaussi <- AICGauss4p
						}
					} else if ((!inherits(Gauss4p, "try-error")) & inherits(Gauss5p, "try-error")) { #### convergence only of Gauss4p
						equalcdG <- TRUE
						Gauss <- Gauss4p
						AICGaussi <- round(AIC(Gauss4p, k = kcrit[4]), digits = AICdigits)
					} else if ((!inherits(Gauss5p, "try-error")) & inherits(Gauss4p, "try-error")) { #### convergence only of Gauss5p
						equalcdG <- FALSE
						Gauss <- Gauss5p
						AICGaussi <- round(AIC(Gauss5p, k = kcrit[5]), digits = AICdigits)
					} else {
						AICGaussi <- Inf
					}
				}


				if (AICGaussi != Inf) {
					modelname = "Gauss-probit"
					AIC <- AICGaussi
					fit <- Gauss

					par <- coef(fit)

					b <- par[["b"]]
					c <- ifelse(equalcdG, par[["d"]], par[["c"]])
					d <- par[["d"]]
					e <- par[["e"]]
					f <- par[["f"]]
					SDres <- sigma(fit)
					nbpari <- ifelse(equalcdLG, 4, 5)
					datapred <- fGauss5p(x = xplot, c = c, d = d, b = b, e = e, f = f)

					predicted <- data.frame(modelname= modelname, x = xplot, predicted=datapred)
					predictedlist[[modelname]] <- predicted

					if (f < 0){
						typology <- "GP.U"
						trend <- "U"
					} else if (f >=0){
						typology <- "GP.bell"
						trend <- "bell"
					}
					result <- data.frame("model" = modelname, "AIC" = AIC, "b" = b, "c" = c, "d" =d, "e" = e, "f" = f, "SDres" =  SDres, "nbpar" = nbpari, "typology" = typology, "trend" = trend )
					resultlist[[modelname]] <- result
				}


				################# LGauss fit ####################
				#NO 2

				if (!lessthan5doses){
					startLGauss5p <- startvalLGauss5pnls(xm = doseu, ym = signalm,Ushape = Ushape)
					LGauss5p <- suppressWarnings(try(nls(formLGauss5p, start = startLGauss5p, data = dset, lower = c(0, -Inf, -Inf, 0, -Inf), algorithm = "port"), silent = TRUE))
				}
				startLGauss4p <- startvalLGauss4pnls(xm = doseu, ym = signalm, Ushape = Ushape)
				LGauss4p <- suppressWarnings(try(nls(formLGauss4p, start = startLGauss4p, data = dset,lower = c(0, -Inf, 0, -Inf), algorithm = "port"), silent = TRUE))
				if (lessthan5doses)	{
					if (!inherits(LGauss4p, "try-error"))	{
						equalcdLG <- TRUE
						LGauss <- LGauss4p
						AICLGaussi <- round(AIC(LGauss4p, k = kcrit[4]), digits = AICdigits)
					} else {
						AICLGaussi <- Inf
					}

				} else { # if (lessthan5doses)
					#### convergence of both models
					if ((!inherits(LGauss4p, "try-error")) & (!inherits(LGauss5p, "try-error")))	{
						AICLGauss4p <- round(AIC(LGauss4p, k = kcrit[4]), digits = AICdigits)
						AICLGauss5p <- round(AIC(LGauss5p, k = kcrit[5]), digits = AICdigits)
						if (AICLGauss5p < AICLGauss4p) {
							LGauss <- LGauss5p
							AICLGaussi <- AICLGauss5p
						} else {
							LGauss <- LGauss4p
							equalcdLG <- TRUE
							AICLGaussi <- AICLGauss4p
						}
					} else if (inherits(LGauss4p, "try-error") & inherits(LGauss5p, "try-error")) { #### no convergence of both models
						AICLGaussi <- Inf
						LGauss <- LGauss5p # we could have given LGauss4p
					} else if ((!inherits(LGauss4p, "try-error")) & inherits(LGauss5p, "try-error")) { #### convergence only of LGauss4p
						equalcdLG <- TRUE
						LGauss <- LGauss4p
						AICLGaussi <- round(AIC(LGauss4p, k = kcrit[4]), digits = AICdigits)
					} else if ((!inherits(LGauss5p, "try-error")) & inherits(LGauss4p, "try-error")) { #### convergence only of LGauss5p
						LGauss <- LGauss5p
						AICLGaussi <- round(AIC(LGauss5p, k = kcrit[5]), digits = AICdigits)
					} else {
						AICLGaussi <- Inf
					}
				}

				if (AICLGaussi != Inf) {
					modelname = "log-Gauss-probit"
					AIC <- AICLGaussi
					fit <- LGauss
					par <- coef(fit)
					b <- par[["b"]]
					c <- ifelse(equalcdLG, par[["d"]], par[["c"]])
					d <- par[["d"]]
					e <- par[["e"]]
					f <- par[["f"]]
					SDres <- sigma(fit)
					nbpari <- ifelse(equalcdLG, 4, 5)
					datapred <- fLGauss5p(x = xplot, c = c, d = d, b = b, e = e, f = f)

					predicted <- data.frame(modelname= modelname, x = xplot, predicted=datapred)
					predictedlist[[modelname]] <- predicted

					if (f < 0) {
						typology <- "lGP.U"
						trend <- "U"
					} else if (f >=0) {
						typology <- "lGP.bell"
						trend <- "bell"
					}

					result <- data.frame("model" = modelname,"AIC" = AIC, "b" = b, "c" = c, "d" =d, "e" = e, "f" = f, "SDres" =  SDres, "nbpar" = nbpari, "typology" = typology, "trend" = trend )
					resultlist[[modelname]] <- result
				}

				################## Hill fit (npar = 4) ##########################
				# No. 3
				startHill <- startvalHillnls2(x = dose, y = signal, xm = doseu, ym = signalm, increase = increaseminmax)
				Hill <- suppressWarnings(try(nls(formHill, start = startHill, data = dset, lower = c(0, -Inf, -Inf, 0), algorithm = "port"), silent = TRUE))
				if (!inherits(Hill, "try-error"))	{
					AICHilli <- round(AIC(Hill, k = kcrit[4]), digits = AICdigits)
				} else {
					AICHilli <- Inf
				}

				if (AICHilli != Inf) {
					modelname = "Hill"
					AIC <- AICHilli
					fit <- Hill
					par <- coef(fit)
					b <- par[["b"]]
					c <- par[["c"]]
					d <- par[["d"]]
					e <- par[["e"]]
					f <- NA
					SDres <- sigma(fit)
					nbpari <- 4
					datapred <- fHill(x = xplot, c = c, d = d, b = b, e = e)
					predicted <- data.frame(modelname= modelname, x = xplot, predicted=datapred)
					predictedlist[[modelname]] <- predicted

					if (c > d) {
						typology <- "H.inc"
						trend <- "inc"
					} else if (c <= d) {
						typology <- "H.dec"
						trend <- "dec"
					}

					result <- data.frame("model" = modelname,"AIC" = AIC, "b" = b, "c" = c, "d" =d, "e" = e, "f" = f, "SDres" =  SDres, "nbpar" = nbpari, "typology" = typology, "trend" = trend )
					resultlist[[modelname]] <- result
				}

				############### Lprobit fit #################
				# NO 4
				startLprobit <- startvalLprobitnls2(x = dose, y = signal, xm = doseu, ym = signalm,	increase = increaseminmax)
				Lprobit <- suppressWarnings(try(nls(formLprobit, start = startLprobit, data = dset,	lower = c(0, -Inf, -Inf, 0), algorithm = "port"), silent = TRUE))
				if (!inherits(Lprobit, "try-error")){
					AICLprobiti <- round(AIC(Lprobit, k = kcrit[4]), digits = AICdigits)
				} else {
					AICLprobiti <- Inf
				}

				if (AICLprobiti != Inf) {
					modelname = "log-probit"
					AIC <- AICLprobiti
					fit <- Lprobit
					par <- coef(fit)
					b <- par[["b"]]
					c <- par[["c"]]
					d <- par[["d"]]
					e <- par[["e"]]
					f <- 0 # to enable the use of the LGauss function to plot the model and calculate the BMD
					SDres <- sigma(fit)
					nbpari <- 4
					datapred <- fLGauss5p(x = xplot, c = c, d = d, b = b, e = e, f = f)
					predicted <- data.frame(modelname= modelname, x = xplot, predicted=datapred)
					predictedlist[[modelname]] <- predicted

					if (c > d) {
						typology <- "lP.inc"
						trend <- "inc"
					} else if (c <= d) {
						typology <- "lP.dec"
						trend <- "dec"
					}

					result <- data.frame("model" = modelname,"AIC" = AIC, "b" = b, "c" = c, "d" =d, "e" = e, "f" = f, "SDres" =  SDres, "nbpar" = nbpari, "typology" = typology, "trend" = trend )
					resultlist[[modelname]] <- result
				}

				################ Expo fit (npar = 3) ###############################
				#NO 5
				# fit of the exponential model with two starting values for abs(e)
				# 0.1 x max(dose) or max(dose)
				startExpo3p.1 <- startvalExp3pnls.1(xm = doseu, ym = signalm,	increase = increaseranks,Ushape = Ushape)
				startExpo3p.2 <- startvalExp3pnls.2(xm = doseu, ym = signalm,	increase = increaseranks,Ushape = Ushape)

				if ((increaseranks & !Ushape) | (!increaseranks & Ushape)){ # e < 0
					# Fit of the 3 par model
					Expo3p.1 <- suppressWarnings(try(nls(formExp3p, start = startExpo3p.1, data = dset,	lower = c(-Inf, -Inf, -Inf), upper = c(Inf, Inf, 0), algorithm = "port"), silent = TRUE))
					Expo3p.2 <- suppressWarnings(try(nls(formExp3p, start = startExpo3p.2, data = dset,	lower = c(-Inf, -Inf, -Inf), upper = c(Inf, Inf, 0), algorithm = "port"), silent = TRUE))
				} else { # e > 0
					# Fit of the 3 par model
					Expo3p.1 <- suppressWarnings(try(nls(formExp3p, start = startExpo3p.1, data = dset,	lower = c(-Inf, -Inf, 0),	algorithm = "port"), silent = TRUE))
					Expo3p.2 <- suppressWarnings(try(nls(formExp3p, start = startExpo3p.2, data = dset,	lower = c(-Inf, -Inf, 0),	algorithm = "port"), silent = TRUE))
				}

				#### convergence of both models
				if ((!inherits(Expo3p.1, "try-error")) & (!inherits(Expo3p.2, "try-error"))){
					AICExpo3p.1 <- round(AIC(Expo3p.1, k = kcrit[3]), digits = AICdigits)
					AICExpo3p.2 <- round(AIC(Expo3p.2, k = kcrit[3]), digits = AICdigits)
					if (AICExpo3p.1 < AICExpo3p.2) {
						Expo <- Expo3p.1
						AICExpoi <- AICExpo3p.1
					} else {
						Expo <- Expo3p.2
						AICExpoi <- AICExpo3p.2
					}
				} else if (inherits(Expo3p.1, "try-error") & inherits(Expo3p.2, "try-error")) { #### no convergence of both models
					AICExpoi <- Inf
					Expo <- Expo3p.1 # we could have given Expo3p.2
				} else if ((!inherits(Expo3p.2, "try-error")) & inherits(Expo3p.1, "try-error")){ #### convergence only of Expo3p.2
					Expo <- Expo3p.2
					AICExpoi <- round(AIC(Expo3p.2, k = kcrit[3]), digits = AICdigits)
				} else if ((!inherits(Expo3p.1, "try-error")) & inherits(Expo3p.2, "try-error")){ #### convergence only of Expo3p.1
					Expo <- Expo3p.1
					AICExpoi <- round(AIC(Expo3p.1, k = kcrit[3]), digits = AICdigits)
				}

				if (AICExpoi != Inf) {
					modelname <- "exponential"
					AIC <- AICExpoi
					fit <- Expo
					par <- coef(fit)
					b <- par[["b"]]
					c <- NA
					d <- par[["d"]]
					e <- par[["e"]]
					f <- NA
					SDres <- sigma(fit)
					nbpari <- 3
					datapred <- fExpo(x = xplot, d = d, b = b, e = e)
					predicted <- data.frame(modelname= modelname, x = xplot, predicted=datapred)
					predictedlist[[modelname]] <- predicted

					if (e > 0 & b > 0) {
						typology <- "E.inc.convex"
						trend <- "inc"
					} else if (e <= 0 & b > 0) {
						typology <- "E.dec.convex"
						trend <- "dec"
					} else if (e <= 0 & b <= 0) {
						typology <- "E.inc.concave"
						trend <- "inc"
					} else if (e > 0 & b <= 0){
						typology <- "E.dec.concave"
						trend <- "dec"
					}

					result <- data.frame("model" = modelname,"AIC" = AIC, "b" = b, "c" = c, "d" =d, "e" = e, "f" = f, "SDres" =  SDres, "nbpar" = nbpari, "typology" = typology, "trend" = trend )
					resultlist[[modelname]] <- result
				}

				######### Fit of the linear model ############################
				# NO 6
				lin <- lm(signal ~ dose, data = dset)
				AIClini <- round(AIC(lin, k = kcrit[2]), digits = AICdigits)

				if (AIClini != Inf) {
					modelname= "linear"
					AIC <- AIClini
					fit <- lin
					par <- coef(fit)
					b <- par[[2]]
					c <- NA
					d <- par[[1]]
					e <- NA
					f <- NA
					SDres <- sigma(fit)
					nbpari <- 2
					datapred <- xplot * b + d
					predicted <- data.frame(modelname= modelname, x = xplot, predicted=datapred)
					predictedlist[[modelname]] <- predicted

					if (b > 0) {
						typology <- "L.inc"
						trend <- "inc"
					} else if (b <= 0) {
						typology <- "L.dec"
						trend <- "dec"
					}

					result <- data.frame("model" = modelname,"AIC" = AIC, "b" = b, "c" = c, "d" =d, "e" = e, "f" = f, "SDres" =  SDres, "nbpar" = nbpari, "typology" = typology, "trend" = trend )
					resultlist[[modelname]] <- result
				}


				######## Fit of the null model (constant) ###########################
				# NO 7
				if (length(resultlist) == 0) {
					constmodel <- lm(signal ~ 1, data = dset)
					AICconsti <-  round(AIC(constmodel, k = kcrit[1]), digits = AICdigits)
					modelname= "constant"
					AIC <- AICconsti
					fit <- constmodel
					par <- NA
					nbpari <- 1
					b <- NA
					c <- mean(dset$signal)
					d <- NA
					e <- NA
					f <- NA
					SDres <- sigma(constmodel)

					result <- data.frame("model" = modelname,"AIC" = AIC, "b" = b, "c" = c, "d" =d, "e" = e, "f" = f, "SDres" =  SDres, "nbpar" = nbpari, "typology" = "NA", "trend" = "NA")
					datapred <- rep(mean(dset$signal), length(xplot))
					resultlist[[modelname]] <- result
					predictedlist[[modelname]] <- predicted

				}

				predicteddf <- do.call(rbind, predictedlist) %>% as.data.frame()
				resultdf <- do.call(rbind, resultlist) %>% as.data.frame() %>%
				dplyr::arrange(AIC) %>%
				dplyr::select(-model) %>%
				dplyr::mutate_if(is.numeric, round, digits = 4)


				fitDAT <- data.frame(conc = dose, expr = signal)

				p <- ggplot(fitDAT, aes(x=conc, y=expr)) +
				geom_point() +
				ylim(min(predicteddf['predicted']), max(predicteddf['predicted'])) +
				geom_line(data = predicteddf, aes(x=x, y=predicted, group=modelname, colour=modelname))+
				ggtitle("Fit Models") +
				theme_bw(base_size = basefontsize) + xlab(xlabel) +  ylab(ylabel) +
				theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), axis.text.x = element_text(angle = 0), legend.title = element_blank(), legend.position="bottom")

				if (logbase != 1){
					p <- p +  scale_x_continuous(trans=scales::pseudo_log_trans(base = logbase))
				}

				if (nrow(resultdf) >0){
					if (npcy == "auto")
					npcy = "top"
					if (npcx == "auto") {
						npcx = "right"
					}

					npc_table <-  geom_table_npc(data = resultdf, label = list(resultdf), table.colnames = TRUE, table.rownames = TRUE,
						npcx = npcx, npcy = npcy,
						table.theme = ttheme_default(base_size = basefontsize, padding = unit(c(1, 1), "mm"))
					)
					p <- p + npc_table
				}

				return(list(plot=p, result=resultdf))
			})

			output$ModelSelection <- renderPlot({
				ModelSelection_out2()[["plot"]]
			})

			output$fitdata <- DT::renderDataTable({
				data_tmp <- DataExpReactive2() %>%
				mutate_if(is.numeric, round, digits = 2)
				DT::datatable(data_tmp, options = list(pageLength = 15))
			})
			
			session$onFlushed(function() {
			  if (!public_dataset) {
			    hideTab(inputId = ns("FittingCurve_tabset"), target = "Data Table")
			  }
			}, once = TRUE)
			
			output$fitresult <- DT::renderDataTable({
				fitresult <- 	ModelSelection_out2()[["result"]]
				fitresult <- fitresult %>%
				mutate_if(is.numeric, round, digits = 4)
				DT::datatable(fitresult)
			})

			###########################################################################################################
			#DRomics fit all data
			observe({
				req(length(working_project()) > 0)
				updateNumericInput(session, 'core', label = "Core will be used:", value = corenumber,  min = 2, max = corenumber, step =1)
				pcutoff <- input$pvalcut2
				datapoint <- input$datapoint2
				psel <- input$psel2

				if (!is.null(psel) & !is.null(DataInSets[[working_project()]]$statresult)) {
					if (psel == "Pval") {
						tmpdat = DataInSets[[working_project()]]$statresult %>% filter(pvalue < pcutoff & n >= datapoint)
					}
					if (psel == "Padj") {
						tmpdat =DataInSets[[working_project()]]$statresult %>% filter(padjust < pcutoff & n >= datapoint)
					}
					UniqueIDnum <- nrow(tmpdat)
				} else {
					data_long <- DataInSets[[working_project()]]$data_long
					req("UniqueID" %in% colnames(data_long) & "conc" %in% colnames(data_long))
					if (("group" %in% colnames(data_long)) && !("treatment" %in% colnames(data_long))) {
						data_long <- data_long %>% dplyr::rename("treatment"="group") %>%
						dplyr::rename("expr" = "response")
					}

					UniqueIDnum <- length(unique(data_long$UniqueID))
				}

				output$filteredgene2 =	 renderText({paste("<font color=\'red\'><b>Total Genes: ", UniqueIDnum, "</b></font>",sep="")})
				if (UniqueIDnum  < 200)
				updateRadioButtons(session, "parallel",  label= "parallel processing (not recommended)",  choices= c("yes"="yes","no"="no"),  selected = "no")
			})

			DROmicsFitOne <- function(x, logbase) {
				data_tmp <- x %>% arrange(conc)
				dose <- data_tmp$conc
				doseranks <- as.numeric(as.factor(dose))
				signal <- data_tmp$expr

				signalm <- data_tmp %>%
				group_by(conc) %>%
				dplyr::summarize(Mean = mean(expr, na.rm=TRUE)) %>%
				arrange(conc)  %>%
				pull(Mean)

				doseu <-  data_tmp %>%
				group_by(conc) %>%
				dplyr::summarize(Mean = mean(expr, na.rm=TRUE)) %>%
				arrange(conc)  %>%
				pull(conc)


				# preparation of data for modelling with nls
				dset <- data.frame(signal = signal, dose = dose, doseranks = doseranks)
				# calculations for starting values and other uses
				dosemin <- min(dose)
				dosemax <- max(dose)
				dosemed <- median(dose[dose!=0])

				# number of points per dose-response curve
				npts <- length(dose)
				ndoses <- length(unique(dose))
				lessthan5doses <- ndoses < 5

				# Information criterion definition
				information.criterion = "AICc"
				AICdigits <- 2 # number of digits for rounding the AIC values
				information.criterion <- match.arg(information.criterion, c("AICc", "BIC", "AIC"))

				# kcrit gives the argument k to pass to function AIC()
				# dependeing of the number of parameters of the model
				# (1 to 5, corresponding to the index of the vector)
				if (information.criterion == "AIC"){
					kcrit <- rep(2, 5)
				} else if (information.criterion == "AICc") {
					nparwithsigma <- 1:5 + 1
					kcrit <- 2*npts /(npts - nparwithsigma - 1)
				} else { # BIC last choice
					lnpts <- log(npts)
					kcrit <- rep(lnpts, 5)
				}

				equalcdG <- FALSE # use to define the value of c equal to d in the Gauss4p model if needed
				equalcdLG <- FALSE # use to define the value of c equal to d in the LGauss4p model if needed

				if (any(!complete.cases(dset))){
					# removing lines with NA values for the signal
					dset <- dset[complete.cases(dset$signal), ]
				}

				# for choice of the linear trend (decreasing or increasing)
				modlin <- lm(signal ~ doseranks, data = dset)
				increaseranks <- coef(modlin)[2] >= 0
				increaseminmax <- dset$dose[which.min(dset$signal)] < dset$dose[which.max(dset$signal)]

				# for choice of the quadratic trend (Ushape or Umbrella shape)
				modquad <- lm(signal ~ doseranks + I(doseranks^2), data = dset)
				Ushape <- coef(modquad)[3] >= 0

				#
				dose_log_transfo = TRUE
				npts = 50

				if (logbase == 10){
					minx <- min(dose[dose != 0])
					maxx <- max(dose)
					xplot <- c(0, 10^seq(log10(minx), log10(maxx), length.out = npts))
				} else {
					xplot <- seq(0, max(dose), length.out = npts)
				}

				resultlist <- list()

				################### Gauss fit ########################
				# No 1
				startGauss4p <- startvalGauss4pnls(xm = doseu, ym = signalm, Ushape = Ushape)
				Gauss4p <- suppressWarnings(try(nls(formGauss4p, start = startGauss4p, data = dset, lower = c(0, -Inf, 0, -Inf), algorithm = "port"), silent = TRUE))

				if (lessthan5doses)	{
					if (!inherits(Gauss4p, "try-error")) {
						equalcdG <- TRUE
						Gauss <- Gauss4p
						AICGaussi <- round(AIC(Gauss4p, k = kcrit[4]), digits = AICdigits)
					} else {
						AICGaussi <- Inf
					}
				} else {
					startGauss5p <- startvalGauss5pnls(xm = doseu, ym = signalm, Ushape = Ushape)
					Gauss5p <- suppressWarnings(try(nls(formGauss5p, start = startGauss5p, data = dset,	lower = c(0, -Inf, -Inf, 0, -Inf), algorithm = "port"), silent = TRUE))

					#### convergence of both models
					if ((!inherits(Gauss4p, "try-error")) & (!inherits(Gauss5p, "try-error")))	{ ## both fitting successfully
						AICGauss4p <- round(AIC(Gauss4p, k = kcrit[4]), digits = AICdigits)
						AICGauss5p <- round(AIC(Gauss5p, k = kcrit[5]), digits = AICdigits)

						if (AICGauss5p < AICGauss4p) {
							Gauss <- Gauss5p
							AICGaussi <- AICGauss5p
						} else {
							Gauss <- Gauss4p
							equalcdG <- TRUE
							AICGaussi <- AICGauss4p
						}
					} else if ((!inherits(Gauss4p, "try-error")) & inherits(Gauss5p, "try-error")) { #### convergence only of Gauss4p
						equalcdG <- TRUE
						Gauss <- Gauss4p
						AICGaussi <- round(AIC(Gauss4p, k = kcrit[4]), digits = AICdigits)
					} else if ((!inherits(Gauss5p, "try-error")) & inherits(Gauss4p, "try-error")) { #### convergence only of Gauss5p
						equalcdG <- FALSE
						Gauss <- Gauss5p
						AICGaussi <- round(AIC(Gauss5p, k = kcrit[5]), digits = AICdigits)
					} else {
						AICGaussi <- Inf
					}
				}


				if (AICGaussi != Inf) {
					modelname = "Gauss-probit"
					AIC <- AICGaussi
					fit <- Gauss
					par <- coef(fit)
					b <- par[["b"]]
					c <- ifelse(equalcdG, par[["d"]], par[["c"]])
					d <- par[["d"]]
					e <- par[["e"]]
					f <- par[["f"]]
					SDres <- sigma(fit)
					nbpari <- ifelse(equalcdLG, 4, 5)

					if (f < 0){
						typology <- "GP.U"
						trend <- "U"
					} else if (f >=0){
						typology <- "GP.bell"
						trend <- "bell"
					}
					result <- data.frame("model" = modelname, "AIC" = AIC, "b" = b, "c" = c, "d" =d, "e" = e, "f" = f, "SDres" =  SDres, "nbpar" = nbpari, "typology" = typology, "trend" = trend )
					resultlist[[modelname]] <- result
				}


				################# LGauss fit ####################
				#NO 2

				if (!lessthan5doses){
					startLGauss5p <- startvalLGauss5pnls(xm = doseu, ym = signalm,Ushape = Ushape)
					LGauss5p <- suppressWarnings(try(nls(formLGauss5p, start = startLGauss5p, data = dset, lower = c(0, -Inf, -Inf, 0, -Inf), algorithm = "port"), silent = TRUE))
				}
				startLGauss4p <- startvalLGauss4pnls(xm = doseu, ym = signalm, Ushape = Ushape)
				LGauss4p <- suppressWarnings(try(nls(formLGauss4p, start = startLGauss4p, data = dset,lower = c(0, -Inf, 0, -Inf), algorithm = "port"), silent = TRUE))
				if (lessthan5doses)	{
					if (!inherits(LGauss4p, "try-error"))	{
						equalcdLG <- TRUE
						LGauss <- LGauss4p
						AICLGaussi <- round(AIC(LGauss4p, k = kcrit[4]), digits = AICdigits)
					} else {
						AICLGaussi <- Inf
					}

				} else { # if (lessthan5doses)
					#### convergence of both models
					if ((!inherits(LGauss4p, "try-error")) & (!inherits(LGauss5p, "try-error")))	{
						AICLGauss4p <- round(AIC(LGauss4p, k = kcrit[4]), digits = AICdigits)
						AICLGauss5p <- round(AIC(LGauss5p, k = kcrit[5]), digits = AICdigits)
						if (AICLGauss5p < AICLGauss4p) {
							LGauss <- LGauss5p
							AICLGaussi <- AICLGauss5p
						} else {
							LGauss <- LGauss4p
							equalcdLG <- TRUE
							AICLGaussi <- AICLGauss4p
						}
					} else if (inherits(LGauss4p, "try-error") & inherits(LGauss5p, "try-error")) { #### no convergence of both models
						AICLGaussi <- Inf
						LGauss <- LGauss5p # we could have given LGauss4p
					} else if ((!inherits(LGauss4p, "try-error")) & inherits(LGauss5p, "try-error")) { #### convergence only of LGauss4p
						equalcdLG <- TRUE
						LGauss <- LGauss4p
						AICLGaussi <- round(AIC(LGauss4p, k = kcrit[4]), digits = AICdigits)
					} else if ((!inherits(LGauss5p, "try-error")) & inherits(LGauss4p, "try-error")) { #### convergence only of LGauss5p
						LGauss <- LGauss5p
						AICLGaussi <- round(AIC(LGauss5p, k = kcrit[5]), digits = AICdigits)
					} else {
						AICLGaussi <- Inf
					}
				}

				if (AICLGaussi != Inf) {
					modelname = "log-Gauss-probit"
					AIC <- AICLGaussi
					fit <- LGauss
					par <- coef(fit)
					b <- par[["b"]]
					c <- ifelse(equalcdLG, par[["d"]], par[["c"]])
					d <- par[["d"]]
					e <- par[["e"]]
					f <- par[["f"]]
					SDres <- sigma(fit)
					nbpari <- ifelse(equalcdLG, 4, 5)

					if (f < 0) {
						typology <- "lGP.U"
						trend <- "U"
					} else if (f >=0) {
						typology <- "lGP.bell"
						trend <- "bell"
					}

					result <- data.frame("model" = modelname,"AIC" = AIC, "b" = b, "c" = c, "d" =d, "e" = e, "f" = f, "SDres" =  SDres, "nbpar" = nbpari, "typology" = typology, "trend" = trend )
					resultlist[[modelname]] <- result
				}

				################## Hill fit (npar = 4) ##########################
				# No. 3
				startHill <- startvalHillnls2(x = dose, y = signal, xm = doseu, ym = signalm, increase = increaseminmax)
				Hill <- suppressWarnings(try(nls(formHill, start = startHill, data = dset, lower = c(0, -Inf, -Inf, 0), algorithm = "port"), silent = TRUE))
				if (!inherits(Hill, "try-error"))	{
					AICHilli <- round(AIC(Hill, k = kcrit[4]), digits = AICdigits)
				} else {
					AICHilli <- Inf
				}

				if (AICHilli != Inf) {
					modelname = "Hill"
					AIC <- AICHilli
					fit <- Hill
					par <- coef(fit)
					b <- par[["b"]]
					c <- par[["c"]]
					d <- par[["d"]]
					e <- par[["e"]]
					f <- NA
					SDres <- sigma(fit)
					nbpari <- 4

					if (c > d) {
						typology <- "H.inc"
						trend <- "inc"
					} else if (c <= d) {
						typology <- "H.dec"
						trend <- "dec"
					}

					result <- data.frame("model" = modelname,"AIC" = AIC, "b" = b, "c" = c, "d" =d, "e" = e, "f" = f, "SDres" =  SDres, "nbpar" = nbpari, "typology" = typology, "trend" = trend )
					resultlist[[modelname]] <- result
				}

				############### Lprobit fit #################
				# NO 4
				startLprobit <- startvalLprobitnls2(x = dose, y = signal, xm = doseu, ym = signalm,	increase = increaseminmax)
				Lprobit <- suppressWarnings(try(nls(formLprobit, start = startLprobit, data = dset,	lower = c(0, -Inf, -Inf, 0), algorithm = "port"), silent = TRUE))
				if (!inherits(Lprobit, "try-error")){
					AICLprobiti <- round(AIC(Lprobit, k = kcrit[4]), digits = AICdigits)
				} else {
					AICLprobiti <- Inf
				}

				if (AICLprobiti != Inf) {
					modelname = "log-probit"
					AIC <- AICLprobiti
					fit <- Lprobit
					par <- coef(fit)
					b <- par[["b"]]
					c <- par[["c"]]
					d <- par[["d"]]
					e <- par[["e"]]
					f <- 0 # to enable the use of the LGauss function to plot the model and calculate the BMD
					SDres <- sigma(fit)
					nbpari <- 4

					if (c > d) {
						typology <- "lP.inc"
						trend <- "inc"
					} else if (c <= d) {
						typology <- "lP.dec"
						trend <- "dec"
					}

					result <- data.frame("model" = modelname,"AIC" = AIC, "b" = b, "c" = c, "d" =d, "e" = e, "f" = f, "SDres" =  SDres, "nbpar" = nbpari, "typology" = typology, "trend" = trend )
					resultlist[[modelname]] <- result
				}

				################ Expo fit (npar = 3) ###############################
				#NO 5
				# fit of the exponential model with two starting values for abs(e)
				# 0.1*max(dose) or max(dose)
				startExpo3p.1 <- startvalExp3pnls.1(xm = doseu, ym = signalm,	increase = increaseranks,Ushape = Ushape)
				startExpo3p.2 <- startvalExp3pnls.2(xm = doseu, ym = signalm,	increase = increaseranks,Ushape = Ushape)

				if ((increaseranks & !Ushape) | (!increaseranks & Ushape)){ # e < 0
					# Fit of the 3 par model
					Expo3p.1 <- suppressWarnings(try(nls(formExp3p, start = startExpo3p.1, data = dset,	lower = c(-Inf, -Inf, -Inf), upper = c(Inf, Inf, 0), algorithm = "port"), silent = TRUE))
					Expo3p.2 <- suppressWarnings(try(nls(formExp3p, start = startExpo3p.2, data = dset,	lower = c(-Inf, -Inf, -Inf), upper = c(Inf, Inf, 0), algorithm = "port"), silent = TRUE))
				} else { # e > 0
					# Fit of the 3 par model
					Expo3p.1 <- suppressWarnings(try(nls(formExp3p, start = startExpo3p.1, data = dset,	lower = c(-Inf, -Inf, 0),	algorithm = "port"), silent = TRUE))
					Expo3p.2 <- suppressWarnings(try(nls(formExp3p, start = startExpo3p.2, data = dset,	lower = c(-Inf, -Inf, 0),	algorithm = "port"), silent = TRUE))
				}

				#### convergence of both models
				if ((!inherits(Expo3p.1, "try-error")) & (!inherits(Expo3p.2, "try-error"))){
					AICExpo3p.1 <- round(AIC(Expo3p.1, k = kcrit[3]), digits = AICdigits)
					AICExpo3p.2 <- round(AIC(Expo3p.2, k = kcrit[3]), digits = AICdigits)
					if (AICExpo3p.1 < AICExpo3p.2) {
						Expo <- Expo3p.1
						AICExpoi <- AICExpo3p.1
					} else {
						Expo <- Expo3p.2
						AICExpoi <- AICExpo3p.2
					}
				} else if (inherits(Expo3p.1, "try-error") & inherits(Expo3p.2, "try-error")) { #### no convergence of both models
					AICExpoi <- Inf
					Expo <- Expo3p.1 # we could have given Expo3p.2
				} else if ((!inherits(Expo3p.2, "try-error")) & inherits(Expo3p.1, "try-error")){ #### convergence only of Expo3p.2
					Expo <- Expo3p.2
					AICExpoi <- round(AIC(Expo3p.2, k = kcrit[3]), digits = AICdigits)
				} else if ((!inherits(Expo3p.1, "try-error")) & inherits(Expo3p.2, "try-error")){ #### convergence only of Expo3p.1
					Expo <- Expo3p.1
					AICExpoi <- round(AIC(Expo3p.1, k = kcrit[3]), digits = AICdigits)
				}

				if (AICExpoi != Inf) {
					modelname <- "exponential"
					AIC <- AICExpoi
					fit <- Expo
					par <- coef(fit)
					b <- par[["b"]]
					c <- NA
					d <- par[["d"]]
					e <- par[["e"]]
					f <- NA
					SDres <- sigma(fit)
					nbpari <- 3

					if (e > 0 & b > 0) {
						typology <- "E.inc.convex"
						trend <- "inc"
					} else if (e <= 0 & b > 0) {
						typology <- "E.dec.convex"
						trend <- "dec"
					} else if (e <= 0 & b <= 0) {
						typology <- "E.inc.concave"
						trend <- "inc"
					} else if (e > 0 & b <= 0){
						typology <- "E.dec.concave"
						trend <- "dec"
					}

					result <- data.frame("model" = modelname,"AIC" = AIC, "b" = b, "c" = c, "d" =d, "e" = e, "f" = f, "SDres" =  SDres, "nbpar" = nbpari, "typology" = typology, "trend" = trend )
					resultlist[[modelname]] <- result
				}

				######### Fit of the linear model ############################
				# NO 6
				lin <- lm(signal ~ dose, data = dset)
				AIClini <- round(AIC(lin, k = kcrit[2]), digits = AICdigits)

				if (AIClini != Inf) {
					modelname= "linear"
					AIC <- AIClini
					fit <- lin
					par <- coef(fit)
					b <- par[[2]]
					c <- NA
					d <- par[[1]]
					e <- NA
					f <- NA
					SDres <- sigma(fit)
					nbpari <- 2

					if (b > 0) {
						typology <- "L.inc"
						trend <- "inc"
					} else if (b <= 0) {
						typology <- "L.dec"
						trend <- "dec"
					}

					result <- data.frame("model" = modelname,"AIC" = AIC, "b" = b, "c" = c, "d" =d, "e" = e, "f" = f, "SDres" =  SDres, "nbpar" = nbpari, "typology" = typology, "trend" = trend )
					resultlist[[modelname]] <- result
				}


				######## Fit of the null model (constant) ###########################
				# NO 7
				if (length(resultlist) == 0) {
					constmodel <- lm(signal ~ 1, data = dset)
					AICconsti <-  round(AIC(constmodel, k = kcrit[1]), digits = AICdigits)
					modelname= "constant"
					AIC <- AICconsti
					fit <- constmodel
					par <- NA
					nbpari <- 1
					b <- NA
					c <- mean(dset$signal)
					d <- NA
					e <- NA
					f <- NA
					SDres <- sigma(constmodel)

					result <- data.frame("model" = modelname,"AIC" = AIC, "b" = b, "c" = c, "d" =d, "e" = e, "f" = f, "SDres" =  SDres, "nbpar" = nbpari, "typology" = "NA", "trend" = "NA")
					resultlist[[modelname]] <- result

				}


				resultdf <- do.call(rbind, resultlist) %>%
				as.data.frame() %>%
				dplyr::arrange(AIC) %>%
				#top_n(1) %>%
				dplyr::filter(row_number()==1) %>%
				dplyr::mutate_if(is.numeric, round, digits = 4) %>%
				dplyr::mutate(UniqueID = unique(data_tmp$UniqueID)) %>%
				dplyr::mutate(treatment = unique(data_tmp$treatment)) %>%
				dplyr::relocate(c(UniqueID,treatment), .before = model)

				return(resultdf)
			}

			MultiCoreFitDROmics <- function(dataS_ChunkByCore,corenum,logbase) {
				cl <- parallel::makeCluster(corenum)
				parallel::clusterEvalQ(cl, {
					library(dplyr)
					source("util-basicandfitfunc.R")
				})
				parallel::clusterExport(cl=cl, varlist=c("dataS_ChunkByCore","DROmicsFitOne","logbase"), envir=environment())

				out <- parallel::parLapply(cl, 1:length(dataS_ChunkByCore), function (k) {
					dataStmp <- dataS_ChunkByCore[[k]]
					Res <- lapply(dataStmp, DROmicsFitOne,logbase)
					coedf <- purrr::map_df(Res, ~as.data.frame(.x), .id="id")
					return(coedf)
				})
				stopCluster(cl)
				return (out)
			}

			fit_all2 <- eventReactive(input$fitallButton2, {
				pcutoff <- input$pvalcut2
				datapoint <- input$datapoint2
				psel <- input$psel2
				parallel <- input$parallel
				corenum <- input$core
				logbase = as.numeric(input$logbaseb)

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
						dataS <- named_group_split(data_long, UniqueID, treatment)
					}

					if (parallel == "yes")  {
						dataS_ChunkByCore <-  split(dataS, cut(seq_along(dataS), corenum, labels = FALSE))
						out <- MultiCoreFitDROmics(dataS_ChunkByCore,corenum,logbase)
						results <- ldply(out, data.frame) %>%
						#dplyr::filter(!is.na(direction_fit)) %>%
						dplyr::mutate_if(is.numeric, round, digits = 4)
					} else {
						Res <- lapply(dataS, DROmicsFitOne,logbase)
						results  <- purrr::map_df(Res, ~as.data.frame(.x), .id="id") %>%
						#dplyr::filter(!is.na(direction))  %>%
						dplyr::mutate_if(is.numeric, round, digits = 4)
					}
				}
				return(results)
			})

			output$results <- DT::renderDataTable({
				results <- as.data.frame("No Fitting Results")
				if (input$fitallButton2[1] == 0) {
					if (!is.null(DataInSets[[working_project()]]$results_omics))
					results <- DataInSets[[working_project()]]$results_omics
				} else {
					withProgress(message = 'Caculating...',  detail = 'This may take a while...',  {
						results <- fit_all2()
						DataInSets[[working_project()]]$results_omics <-  results

						if (input$saveproject2 == 1) {
							shiny::validate(need(input$projectname2!= "","Please provide project name."))
							filename <- paste("data/",input$projectname2,".RData",sep="")
							save(data_long, results, file=filename )
						}

					})
				}
				results <- results %>%  mutate_if(is.numeric, round, digits = 2)
				DT::datatable(results, options = list(pageLength = 15), rownames= FALSE, filter = 'top')
			})

			###########################################################################################################
			observe({
				req(length(working_project()) > 0)
				updateSelectInput(session,'sel_page2', choices= seq_len(100))
				validate(
					need(DataInSets[[working_project()]]$results_omics, "Need fitting results")
				)
				sel_treatment <- input$sel_treatment2b

				typology <- DataInSets[[working_project()]]$results_omics %>%
				dplyr::filter(!is.na(typology)) %>%
				dplyr::filter(typology != "NA") %>%
				dplyr::filter(treatment %in% sel_treatment) %>%
				dplyr::pull(typology) %>% unique()

				trend <- DataInSets[[working_project()]]$results_omics %>%
				dplyr::filter(!is.na(trend)) %>%
				dplyr::filter(trend != "NA") %>%
				dplyr::filter(treatment %in% sel_treatment) %>%
				dplyr::pull(trend) %>% unique()

				updateSelectInput(session,'typology', choices=typology,selected=typology)
				updateSelectInput(session,'trend', choices=trend,selected=trend)
			})

			observe({
				req(length(working_project()) > 0)
				updateSelectInput(session,'sel_page', choices= seq_len(100))
			})

			browsing_out <- reactive({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$data_long)

				validate(
					need(DataInSets[[working_project()]]$results_omics, "Need fitting results")
				)
				validate(
					need(DataInSets[[working_project()]]$data_long, "Need data")
				)

				results_omics <- DataInSets[[working_project()]]$results_omics
				
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
				  results_omics <- results_omics %>%
				    dplyr::filter(UniqueID %in% gene_list)
				}
				
				
				data_long <- DataInSets[[working_project()]]$data_long
				if (("group" %in% colnames(data_long)) && !("treatment" %in% colnames(data_long))) {
					data_long <- data_long %>% dplyr::rename("treatment"="group") %>%
					dplyr::rename("expr" = "response")
				}
				

				sel_treatment <- input$sel_treatment2b
				labelfontsize <- as.numeric(input$labelfontsize)
				basefontsize <- as.numeric(input$basefontsize)
				xlabel <- input$xlabel
				ylabel <- input$ylabel
				logbase = as.numeric(input$logbase)

				sel_typology <- input$typology
				sel_trend <- input$trend
				genelabel <- input$sel_geneid

				ncol = as.numeric(input$plot_ncol)
				nrow = as.numeric(input$plot_nrow)
				sel_page <- as.numeric(input$sel_page)-1
				numberpage = ncol * nrow
				startslice = sel_page * numberpage  + 1
				endslice = startslice + numberpage -1

				sliced_id <- results_omics %>%
				dplyr::filter((treatment %in% sel_treatment) & (typology %in% sel_typology) & (trend %in% sel_trend))  %>%
				dplyr::distinct(UniqueID) %>%
				dplyr::slice(startslice:endslice)

				sliced_df <- results_omics%>%
				dplyr::filter((treatment %in% sel_treatment) & (typology %in% sel_typology) & (trend %in% sel_trend))  %>%
				dplyr::filter(UniqueID %in% sliced_id$UniqueID) %>%
				tidyr::unite(id, c("UniqueID","treatment"), remove = FALSE, sep = "-")

				data_long_tmp  <- data_long %>% as.data.frame() %>%
				tidyr::unite(id, c("UniqueID","treatment"), remove = FALSE, sep = "-") %>%
				dplyr::filter(id %in% sliced_df$id)

				###################
				npts = 50
				dose = unique(data_long_tmp$conc)
				if (logbase == 10){
					minx <- min(dose[dose != 0])
					maxx <- max(dose)
					xplot <- c(0, 10^seq(log10(minx), log10(maxx), length.out = npts))
				} else {
					xplot <- seq(0, max(dose), length.out = npts)
				}
				####################
				predictedlist <- list()

				for (row in 1:nrow(sliced_df)) {
					modelname <- sliced_df[row, "model"]
					UniqueID <- sliced_df[row, "UniqueID"]
					treatment <- sliced_df[row, "treatment"]

					c <- sliced_df[row, "c"]
					d <- sliced_df[row, "d"]
					b <- sliced_df[row, "b"]
					e <- sliced_df[row, "e"]
					f <- sliced_df[row, "f"]
					### Gauss fit
					if (modelname ==  "Gauss-probit")
					datapred <- fGauss5p(x = xplot, c = c, d = d, b = b, e = e, f = f)
					### LGauss fit
					if (modelname == "log-Gauss-probit")
					datapred <- fLGauss5p(x = xplot, c = c, d = d, b = b, e = e, f = f)
					### Hill fit (npar = 4)
					if (modelname == "Hill")
					datapred <- fHill(x = xplot, c = c, d = d, b = b, e = e)
					### Lprobit fit
					if (modelname == "log-probit")
					datapred <- fLGauss5p(x = xplot, c = c, d = d, b = b, e = e, f = f)
					### Expo fit (npar = 3)
					if (modelname == "exponential")
					datapred <- fExpo(x = xplot, d = d, b = b, e = e)
					### Fit of the linear model
					if (modelname == "linear")
					datapred <- xplot * b + d
					### Fit of the null model (constant)
					if(modelname == "constant")
					datapred <- rep(mean(dset$signal), length(xplot))

					predicteddf <- data.frame(UniqueID = UniqueID, treatment = treatment,  x = xplot, predicted = datapred)
					predictedlist[[row]]  = predicteddf
				}

				predicteddf2 <- do.call(rbind,predictedlist)

				data_long_tmp <- data_long_tmp %>%
				dplyr::select(UniqueID, treatment, conc, expr)%>% as.data.frame()
				colnames(predicteddf2) <- colnames(data_long_tmp)

				gg.df <- rbind(cbind(geom="pt", data_long_tmp), cbind(geom="ln",predicteddf2)) %>%
				dplyr::rename(GroupName = treatment)

				p <- ggplot(gg.df, aes(x=conc, y=expr, color=GroupName)) +
				geom_point(data=gg.df[gg.df$geom=="pt",], shape=4) +
				geom_line(data=gg.df[gg.df$geom=="ln",]) +
				facet_wrap(~ UniqueID, scales = "free", nrow = nrow, ncol = ncol) +
				theme_bw(base_size = basefontsize) + xlab(xlabel) +  ylab(ylabel) +
				theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), 
					axis.text.x = element_text(angle = 0),
					legend.title = element_blank(), 
					strip.text.x = element_text(size=input$labelfontsize),
					legend.position="bottom")

				if (logbase != 1){
					p <- p +  scale_x_continuous(trans=scales::pseudo_log_trans(base = logbase))
				}

				return(list(plot=p, result = sliced_df))

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
