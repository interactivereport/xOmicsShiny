###########################################################################################################
##Omics Dose-Response/Time-Course Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: drcmodule.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 10/24/2024
##@version 3.0
###########################################################################################################

##########################################################################################################
## Curve Fitting Plot
##########################################################################################################
#pkgs: "drc", "ggpmisc", "parallel", "dplyr", "scales", "DT", "purrr", "modelr", "tibble"

library(drc)
library(ggpmisc)
library(parallel)
library(gridExtra)

#corenumber  = parallel::detectCores(logical = FALSE) - 2
corenumber  = 20 #fixed based on cluster node

fctList <- list(LL.2(), LL.3(), LL.3u(), LL.4(), LL.5(),
	W1.2(), W1.3(), W1.4(), W2.2(), W2.3(), W2.4()
	#BC.4(),
	#BC.5(),
	#LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5(),
	#AR.2(), AR.3(),
	#MM.2(), MM.3()
)

modelnames = c("LL.2:Log-logistic(lower=0; upper=1)"="LL.2",
	"LL.3:Log-logistic(lower=0)"="LL.3",
	"LL.3u:Log-logistic(upper=1)"="LL.3u",
	"LL.4:Log-logistic"="LL.4",
	"LL.5:Generalized log-logistic"="LL.5",
	"W1.2:Weibull(type 1; lower=0; upper=1)"="W1.2",
	"W1.3:Weibull(type 1; lower=0)"="W1.3",
	"W1.4:Weibull(type 1)"="W1.4",
	"W2.2:Weibull(type 2; lower=0; upper=1"="W2.2",
	"W2.3:Weibull(type 2; lower=0)"="W2.3",
	"W2.4:Weibull(type 2)"="W2.4"
	#"BC.4:Brain-Cousens (hormesis; lower=0)"="BC.4",
	#"BC.5:Brain-Cousens (hormesis)"="BC.5",
	#"LL2.2:Log-logistic(log(ED50);lower=0; upper=1)"="LL2.2",
	#"LL2.3:Log-logistic(log(ED50);lower=0)"="LL2.3",
	#"LL2.3u:Log-logistic(log(ED50); upper=1)"="LL2.3u",
	#"LL2.4:Log-logistic(log(ED50))"="LL2.4",
	#"LL2.5:Generalised log-logistic(log(ED50))"="LL2.5",
	#"AR.2:Asymptotic regression(lower=0)"="AR.2",
	#"AR.3:Shifted asymptotic regression"="AR.3",
	#"MM.2:Michaelis-Menten"="MM.2",
	#"MM.3:Shifted Michaelis-Menten"="MM.3"
)

names(fctList) <- modelnames


plotfitting <- function(model, treatment, TimeDose="conc",  npcx="auto", npcy="auto", label = "ED50",xlabel="Log10(conc)", ylabel="Response", basefontsize=14, logbase=10) {
	parameterdf <- data.frame(coename = c("b:(Intercept)", "c:(Intercept)", "d:(Intercept)","e:(Intercept)","f:(Intercept)"), parameter = c("Slope", "Lower Limit", "Upper Limit", "ED50","f") )
	coedf <- data.frame(coename = names(model$coefficients), val = unname(model$coefficients)) %>% left_join(parameterdf, by = "coename") %>% dplyr::select(one_of("parameter","val"))

	rSquared <- cor(model$predres[,1], model$data$expr)^2
	#Pvalue <- noEffect(model)[3]
	#Pvalue <- modelFit(model)$"p value"[2]
	coedf <-  rbind(coedf, data.frame(parameter ="r2",val = rSquared)) %>% mutate_if(is.numeric, round, digits = 2)
	coedft <-  setNames(data.frame(t(coedf[,-1])), coedf[,1])

	fitDAT <- model$origData
	xvalue <- model$dataList$dose
	newdata <- data.frame(DOSE = 0:max(xvalue), CURVE = rep(1,length(0:max(xvalue))))
	predicted <- data.frame(x = c(0:max(xvalue)), predicted= predict(model,newdata))

	p <- ggplot(fitDAT, aes(x=conc, y=expr)) +
	geom_point() +
	ylim(0, max(fitDAT['expr'])) +
	geom_line(color='red', data = predicted, aes(x=x, y=predicted)) +
	#geom_hline(yintercept = 0.5, linetype="dashed",  color = "blue", size=1) +
	ggtitle(treatment) +
	theme_bw(base_size = basefontsize) + xlab(xlabel) +  ylab(ylabel) +
	theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), axis.text.x = element_text(angle = 0), legend.title = element_blank(), legend.position="bottom")

	if (logbase != 1){
		p <- p +  scale_x_continuous(trans=scales::pseudo_log_trans(base = logbase))
	}

	if (npcy == "auto")
	npcy = "top"

	if (npcx == "auto") {
		if (coedf[(coedf$parameter == "Slope"),"val"] > 0)
		npcx = "right"
		else
		npcx = "left"
	}

	npc_table <-  geom_table_npc(data = coedf, label = list(coedf), table.colnames = FALSE,
		npcx = npcx, npcy = npcy,	table.theme = ttheme_default(base_size = basefontsize, padding = unit(c(1, 1), "mm"))
	)
	p <- p + npc_table

	return(list(plot=p, result=	coedft))
}

######
drc_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
		column(3,
			wellPanel(
				uiOutput(ns('loadedprojects')),
				conditionalPanel(ns = ns, "input.upper_tabset=='Fitting Curve'",
					selectizeInput(ns("sel_gene1"),	label="Gene Name",	choices = NULL,	multiple=FALSE, options = list(placeholder =	'Type to search')),
					conditionalPanel(ns = ns, "input.FittingCurve_tabset=='Fitting Curve'",
						radioButtons(ns("separateplot"), label="Plot Separately ", inline = TRUE, choices = c("Yes" = "Yes","No" = "No")),
						selectizeInput(ns("sel_model"),	label="Model",	choices = NULL,	multiple=FALSE)
					),
					conditionalPanel(ns = ns, "input.FittingCurve_tabset=='Model Selection'",
						selectizeInput(ns("sel_treatment1b"),	label="Treatment",	choices = NULL,	multiple=FALSE),
						sliderInput(ns("sel_topn"), "Top Model Number:", min = 1, max = 11, step = 1, value = 5)
					)
				),
				conditionalPanel(ns = ns,"input.upper_tabset=='Browsing'",
					radioButtons(ns("subset"), label="Genes Used in Plot", choices=c("Browsing", "Upload Genes"), inline = TRUE, selected="Browsing"),
					conditionalPanel(ns = ns, "input.subset=='Upload Genes'",	textAreaInput(ns("uploadlist"), "Enter Gene List", "", cols = 5, rows=6)),
					radioButtons(ns("sel_geneid"), label="Select Gene Label", inline = TRUE, choices=c("UniqueID", "Gene.Name","Protein.ID"), selected="UniqueID"),
					radioButtons(ns("orderby"), label="Sort by", inline = TRUE, choices =  c("p_value" = "p_value", "padjust" = "padjust","r.squared" = "r.squared")),
					fluidRow(
						column(width=4,sliderInput(ns("plot_ncol"), label= "Column Number", min = 1, max = 6, step = 1, value = 3)),
						column(width=4,sliderInput(ns("plot_nrow"), label= "Row Number", min = 1, max = 9, step = 1, value = 3)),
						column(width=4,selectInput(ns("sel_page"), label="Select Page",	choices = NULL,	selected=1))
					)
				),
				conditionalPanel(ns = ns,"input.upper_tabset=='Fitting Curve' || input.upper_tabset=='Browsing'",
					conditionalPanel(ns = ns, "input.FittingCurve_tabset!=='Model Selection'",
						selectizeInput(ns("sel_treatment1"),	label="Treatment",	choices = NULL,	multiple=TRUE)
					),
					numericInput(ns("upper"), "Upper:", 1, min = 1, max = 100),
					radioButtons(ns("npcx"), label="Label X Position", inline = TRUE, choices =  c("auto" = "auto", "left" = "left","middle"="middle","right"="right")),
					radioButtons(ns("npcy"), label="Label Y Position", inline = TRUE, choices =  c("auto" = "auto", "top" = "top","center" = "center","bottom" = "bottom")),
					radioButtons(ns("logbase"), label="Log tansform", inline = TRUE,  choices=c("Non"= 1, "10"= 10), selected=10),
					fluidRow(
						column(width=6,textInput(ns("xlabel"), "xlabel", value = "Log10(Conc.)")),
						column(width=6,textInput(ns("ylabel"), "ylabel", value = "Response"))
					),
					fluidRow(
						column(width=6,sliderInput(ns("labelfontsize"), "Label Font Size:", min = 10, max = 24, step = 2, value = 14)),
						column(width=6,sliderInput(ns("basefontsize"), "Font Size:", min = 10, max = 24, step = 2, value = 14))
					),
					sliderInput(ns("refreshrate"), "Plot Refresh Rate:", min = 0, max = 8000, step = 1000, value = 1000)
				),
				conditionalPanel(ns = ns,"input.upper_tabset=='Result Table (all)'",
					radioButtons(ns("parallel"), label= "parallel processing?", choices= c("yes"="yes","no"="no"),inline = TRUE),
					numericInput(ns("core"), "Core will be used:", value = 4, min = 2, max = 12),
					radioButtons(ns("psel1"), label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
					numericInput(ns("pvalcut1"), label= "Choose P-value Threshold",  value=0.01, min=0, step=0.001),
					numericInput(ns("datapoint1"), label= "Minimal Data Points",  value=5, min=5, step=1),
					selectizeInput(ns("sel_model2"),	label="Model",	choices = NULL,	multiple=FALSE),
					numericInput(ns("upper1"), "Upper:", 1, min = 1, max = 100),
					uiOutput(ns("filteredgene1")),
					radioButtons(ns("saveproject1"), label="Save Result?", inline = TRUE,  choices=c("Yes"= 1, "No"= 0), selected=0),
					textInput(ns("projectname1"), "Project Name"),
					actionButton(ns("fitallButton1"), "Fit filtered data")
				)
			)
		),
		column(9,
			tabsetPanel(id=ns("upper_tabset"),
				tabPanel(title="Fitting Curve", value ="Fitting Curve",
					tabsetPanel(id=ns("FittingCurve_tabset"),
						tabPanel(title="Fitting Curve", value="Fitting Curve", plotOutput(ns("FittingCurve"), height=800)
						),
						tabPanel(title="Result Table", value="Result Table", DT::dataTableOutput(ns("fitresult"))
						),
						tabPanel(title="Data Table", value ="Data Table", DT::dataTableOutput(ns("fitdata"))),
						tabPanel(title="Model Selection", value ="Model Selection",
							fluidRow(
								plotOutput(ns("ModelSelection"), height=800)
							),
							fluidRow(
								DT::dataTableOutput(ns("fitresultb"))
							)
						)
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
				tabPanel(title="Help", htmlOutput('help_drc'))
			)
		)
	)
}

drc_server <- function(id) {
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
				data_long <- DataInSets[[working_project()]]$data_long
				req("UniqueID" %in% colnames(data_long) & "conc" %in% colnames(data_long))

				if (("group" %in% colnames(data_long)) && !("treatment" %in% colnames(data_long))) {
					data_long <- data_long %>% dplyr::rename("treatment"="group") %>%
					dplyr::rename("expr" = "response")
				}

				DataIngenes <-  data_long  %>% dplyr::pull(UniqueID) %>% unique() %>% as.character()
				updateSelectizeInput(session,'sel_gene1', choices= DataIngenes, server=TRUE)
				treatment <-  data_long %>% dplyr::pull(treatment) %>% unique() %>% as.character()
				updateSelectizeInput(session,'sel_treatment1', choices= treatment,  selected=treatment)
				updateSelectizeInput(session,'sel_treatment1b', choices= treatment,  selected=treatment[1])
				updateSelectizeInput(session,'sel_model', choices= modelnames,  selected="LL.4")
				updateSelectizeInput(session,'sel_model2', choices= modelnames,  selected="LL.4")
				modelnames2 <- c("Gauss-probit", "log-Gauss-probit", "Hill", "log-probit", "exponential", "linear")
				updateSelectizeInput(session,'sel_model2b', choices= modelnames2,  selected="Hill")
			})

			InputReactive <- reactive({
				sel_gene <- input$sel_gene1
				sel_treatment <- input$sel_treatment1
				sel_model <- input$sel_model
				npcx <- input$npcx
				npcy <- input$npcy
				labelfontsize <- input$labelfontsize
				basefontsize <- input$basefontsize
				logbase = as.numeric(input$logbase)
				xlabel = input$xlabel
				ylabel = input$ylabel
				upperlimit = input$upper
				lowerlimit = input$lowerr

				return(list("sel_gene" =sel_gene,
					"sel_treatment" =sel_treatment,
					"sel_model" =sel_model,
					"npcx" =npcx,
					"npcy" =npcy,
					"labelfontsize" =labelfontsize,
					"basefontsize" =basefontsize,
					"logbase" = logbase,
					"xlabel" = xlabel,
					"ylabel" = ylabel,
				"upperlimit" = upperlimit))
			})

			lookupDRC <- c("Lower.Limit" = "Lower Limit", "Upper.Limit" = "Upper Limit")

			DataExpReactive <- reactive({
				req(length(working_project()) > 0)
				data_long <- DataInSets[[working_project()]]$data_long
				if (("group" %in% colnames(data_long)) && !("treatment" %in% colnames(data_long))) {
					data_long <- data_long %>% dplyr::rename("treatment"="group") %>%
					dplyr::rename("expr" = "response")
				}

				shiny::validate(need(input$sel_gene1 != "","Please select a gene."))
				refreshrate <- as.numeric(input$refreshrate)
				InputData <- InputReactive %>% debounce(refreshrate)
				sel_gene <- InputData()$sel_gene
				data_tmp = dplyr::filter(data_long, UniqueID ==sel_gene)
				return(list("data_tmp"=data_tmp))
			})

			fitting_out <- reactive({
				refreshrate <- as.numeric(input$refreshrate)
				InputData <- InputReactive %>% debounce(refreshrate)

				sel_gene = InputData()[["sel_gene"]]
				sel_treatment = InputData()[["sel_treatment"]]
				sel_model = InputData()[["sel_model"]]
				npcx = InputData()[["npcx"]]
				npcy = InputData()[["npcy"]]
				labelfontsize = InputData()[["labelfontsize"]]
				basefontsize = InputData()[["basefontsize"]]
				logbase = InputData()[["logbase"]]
				xlabel = InputData()[["xlabel"]]
				ylabel = InputData()[["ylabel"]]
				upperlimit = InputData()[["upperlimit"]]

				if (upperlimit != 1) {
					fctList <- list(LL.2(upper = upperlimit), LL.3(), LL.3u(upper = upperlimit), LL.4(), LL.5(),
						W1.2(upper = upperlimit), W1.3(), W1.4(), W2.2(), W2.3(), W2.4()
					)

					names(fctList) <- c("LL.2:Log-logistic(lower=0; upper=1)"="LL.2",
						"LL.3:Log-logistic(lower=0)"="LL.3",
						"LL.3u:Log-logistic(upper=1)"="LL.3u",
						"LL.4:Log-logistic"="LL.4",
						"LL.5:Generalized log-logistic"="LL.5",
						"W1.2:Weibull(type 1; lower=0; upper=1)"="W1.2",
						"W1.3:Weibull(type 1; lower=0)"="W1.3",
						"W1.4:Weibull(type 1)"="W1.4",
						"W2.2:Weibull(type 2; lower=0; upper=1"="W2.2",
						"W2.3:Weibull(type 2; lower=0)"="W2.3",
						"W2.4:Weibull(type 2)"="W2.4"
					)
				}

				DataExpReactive_tmp <- DataExpReactive()
				data_tmp <- DataExpReactive_tmp[["data_tmp"]]

				sel_treatment = intersect(sel_treatment, unique(data_tmp[['treatment']]))
				#sel_treatment = unique(data_tmp[['treatment']])
				out <- list()
				plist <- list()
				#models <- list()
				fitDATlist <- list()
				predictedlist <- list()

				for (onetreatment in sel_treatment) {
					dfgene1 <- data_tmp %>% as.data.frame() %>% dplyr::filter(treatment == onetreatment)
					model <- try(drm(expr~conc, data = dfgene1, fct = fctList[[sel_model]]), silent = TRUE)
					df.summary <- dfgene1 %>%
					group_by(conc, treatment) %>%
					summarise(
						sd = sd(expr),
						expr = mean(expr), .groups = "drop"
					)
					if (class(model) == "drc") {
						#models[[onetreatment]] <- model
						fitDATlist[[onetreatment]] = model$origData
						xvalue <- model$dataList$dose
						newdata <- data.frame(DOSE = 0:max(xvalue), CURVE = rep(1,length(0:max(xvalue))))
						predicted <- data.frame(x = c(0:max(xvalue)), predicted= predict(model,newdata))
						######
						#predslm = predict(model, newdata, interval = "confidence")
						#######
						predictedlist[[onetreatment]] <- predicted
						fitting_PR <- plotfitting(model, treatment=onetreatment, TimeDose="conc", label = "ED50", npcx, npcy, xlabel, ylabel, basefontsize, logbase)
						p <- fitting_PR$plot
						result <- fitting_PR$result
						out[[onetreatment]] <- result
					}	else {
						p <- ggplot(dfgene1, aes(x=conc, y=expr)) +
						geom_point() +
						ylim(0, max(dfgene1['expr'])) +
						ggtitle(onetreatment) +
						theme_bw(base_size = basefontsize) + xlab(xlabel) +  ylab(ylabel) +
						theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), axis.text.x = element_text(angle = 0),legend.title = element_blank(), legend.position="bottom")
						if (logbase != 1){
							p <- p +  scale_x_continuous(trans=scales::pseudo_log_trans(base = logbase))
						}
					}
					plist[[onetreatment]]  <- p
				}
				outdf <- do.call(rbind, out)

				if (input$separateplot == "Yes") {
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
					predicteddf <- cbind(treatment=rep(names(predictedlist), sapply(predictedlist,nrow)),do.call(rbind,predictedlist))
					fitDATdf <- do.call(rbind, fitDATlist)

					p <- ggplot(fitDATdf, aes(x=conc, y=expr, color=treatment)) +
					#geom_point() +
					geom_jitter(aes(color = treatment),position = position_jitter(0.2)) +
					ylim(0, max(fitDATdf['expr']))

					p <- p +
					geom_line(data = predicteddf, aes(x=x, y=predicted, color=treatment))

					p <- p +
					theme_bw(base_size = basefontsize) + xlab(xlabel) +  ylab(ylabel) +
					theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), axis.text.x = element_text(angle = 0), legend.title = element_blank(), legend.position="bottom")
					if (logbase != 1) {
						p <- p +  scale_x_continuous(trans=scales::pseudo_log_trans(base = logbase))
					}

					if (npcy == "auto")
					npcy = "top"

					if (npcx == "auto") {
						if (colMeans(outdf["Slope"], na.rm = FALSE, dims = 1) > 0)
						npcx = "right"
						else
						npcx = "left"
					}
					npc_table <-  geom_table_npc(data = outdf, label = list(outdf), table.rownames = TRUE,	npcx = npcx, npcy = npcy,
						table.theme = ttheme_default(base_size = basefontsize, padding = unit(c(1, 1), "mm"))
					)
					ml <- p + npc_table
				}
				return(list(plot=ml, result=outdf))

			})

			ModelSelection_out <- reactive({
				sel_gene <- input$sel_gene1
				sel_treatment2 <- input$sel_treatment1b

				npcx <- input$npcx
				npcy <- input$npcy
				labelfontsize <- input$labelfontsize
				basefontsize <- input$basefontsize
				logbase = as.numeric(input$logbase)
				xlabel = input$xlabel
				ylabel = input$ylabel
				topn = as.numeric(input$sel_topn)
				upperlimit = as.numeric(input$upper)

				DataExpReactive_tmp <- DataExpReactive()
				data_tmp <- DataExpReactive_tmp[["data_tmp"]]

				if (is.null(data_tmp$treatment)) {
					dfgene <- data_tmp %>%
					dplyr::select(expr, conc) %>%
					as.data.frame()
				} else {
					sel_treatment2 = intersect(sel_treatment2, unique(data_tmp[['treatment']]))
					dfgene <- data_tmp %>% dplyr::filter(treatment == sel_treatment2) %>%
					dplyr::select(expr, conc) %>%
					as.data.frame()
				}

				modelLL4 <- try(drm(expr~conc, data = dfgene, fct = LL.4()), silent = TRUE)

				if (upperlimit != 1) {
					fctList <- list(LL.2(upper = upperlimit), LL.3(), LL.3u(upper = upperlimit), LL.4(), LL.5(),
						W1.2(upper = upperlimit), W1.3(), W1.4(), W2.2(), W2.3(), W2.4()
					)

					names(fctList) <- c("LL.2:Log-logistic(lower=0; upper=1)"="LL.2",
						"LL.3:Log-logistic(lower=0)"="LL.3",
						"LL.3u:Log-logistic(upper=1)"="LL.3u",
						"LL.4:Log-logistic"="LL.4",
						"LL.5:Generalized log-logistic"="LL.5",
						"W1.2:Weibull(type 1; lower=0; upper=1)"="W1.2",
						"W1.3:Weibull(type 1; lower=0)"="W1.3",
						"W1.4:Weibull(type 1)"="W1.4",
						"W2.2:Weibull(type 2; lower=0; upper=1"="W2.2",
						"W2.3:Weibull(type 2; lower=0)"="W2.3",
						"W2.4:Weibull(type 2)"="W2.4"
					)
				}


				if (!inherits(modelLL4, "try-error")){
					if (npcy == "auto")
					npcy = "top"

					if (npcx == "auto") {
						if (modelLL4$coefficients["b:(Intercept)"] > 0)
						npcx = "right"
						else
						npcx = "left"
					}

					fctList2 <- fctList
					fctList2[['LL.4']] <- NULL

					sel <- try(drc::mselect(modelLL4, fctList2, sorted = "IC", linreg = FALSE), silent = TRUE)
					sel <- sel %>% as.data.frame() %>% top_n(., -topn, IC) %>% mutate_if(is.numeric, round, digits = 6)
					fitDAT = modelLL4$origData
					xvalue <- modelLL4$dataList$dose
					randomdose <- 0:max(xvalue)
					newdata <- data.frame(DOSE = randomdose, CURVE = rep(1,length(randomdose)))
					predictedlist <- list()
					for (modelname in dimnames(sel)[[1]]){
						model <- try(drm(expr~conc, data = dfgene, fct = fctList[[modelname]]), silent = TRUE)
						if (!inherits(model, "try-error")) {
							predicted <- data.frame(modelname=rep(modelname,length(randomdose)), x = randomdose, predicted= predict(model, newdata))
							predictedlist[[modelname]] <- predicted
						}
					}
					predicteddf <- do.call(rbind, predictedlist)

					p <- ggplot(fitDAT, aes(x=conc, y=expr)) +
					geom_point() +
					ylim(0, max(predicteddf['predicted'])) +
					geom_line(data = predicteddf, aes(x=x, y=predicted, group=modelname, colour=modelname))+
					ggtitle("Top Models") +
					theme_bw(base_size = basefontsize) + xlab(xlabel) +  ylab(ylabel) +
					theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), axis.text.x = element_text(angle = 0), legend.title = element_blank(), legend.position="bottom")

					if (logbase != 1){
						p <- p +  scale_x_continuous(trans=scales::pseudo_log_trans(base = logbase))
					}

					npc_table <-  geom_table_npc(data = sel, label = list(sel), table.colnames = TRUE, table.rownames = TRUE,
						npcx = npcx, npcy = npcy,
						table.theme = ttheme_default(base_size = basefontsize, padding = unit(c(1, 1), "mm"))
					)
					p <- p + npc_table
					return(list(plot=p, result=sel))
				}	else {
					p <- ggplot(dfgene, aes(x=conc, y=expr)) +
					geom_point() +
					ylim(0, max(dfgene['expr'])) +
					#ggtitle("") +
					theme_bw(base_size = basefontsize) + xlab(xlabel) +  ylab(ylabel) +
					theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), axis.text.x = element_text(angle = 0),legend.title = element_blank(), legend.position="bottom")
					if (logbase != 1){
						p <- p +  scale_x_continuous(trans=scales::pseudo_log_trans(base = logbase))
					}
					sel <- data.frame()
					return(list(plot=p, result=sel))
				}
			})

			output$FittingCurve <- renderPlot({
				fitting_out()[["plot"]]
			})

			output$fitresult <- DT::renderDataTable({
				fitresult <- fitting_out()[["result"]]
				fitresult <- fitresult %>%
				mutate_if(is.numeric, round, digits = 4)
				DT::datatable(fitresult, options = list(pageLength = 15))
			})

			observeEvent(input$FittingCurve, {
				saved.num <- length(saved_plots$FittingCurve) + 1
				saved_plots$FittingCurve[[saved.num]] <- fitting_out()[["plot"]]
			})

			output$ModelSelection <- renderPlot({
				ModelSelection_out()[["plot"]]
			})

			output$fitresultb <- DT::renderDataTable({
				fitresult <- ModelSelection_out()[["result"]]
				fitresult <- fitresult %>%
				mutate_if(is.numeric, round, digits = 4)
				DT::datatable(fitresult, options = list(pageLength = 15))
			})

			output$fitdata <- DT::renderDataTable({
				data_tmp <- DataExpReactive()[["data_tmp"]]
				data_tmp <- data_tmp %>%
				mutate_if(is.numeric, round, digits = 2)
				DT::datatable(data_tmp, options = list(pageLength = 15))
			})
			
			session$onFlushed(function() {
			  if (!public_dataset) {
			    hideTab(inputId = ns("FittingCurve_tabset"), target = "Data Table")
			  }
			}, once = TRUE)
			
			###########################################################################################################
			#DRC fit all data
			observe({
				req(length(working_project()) > 0)

				updateNumericInput(session, 'core', label = "Core will be used:", value = corenumber,  min = 2, max = corenumber, step =1)
				pcutoff <- input$pvalcut1
				datapoint <- input$datapoint1
				psel <- input$psel1
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
				output$filteredgene1 =	 renderText({paste("<font color=\'red\'><b>Total Genes: ", UniqueIDnum, "</b></font>",sep="")})
				if (UniqueIDnum  < 200)
				updateRadioButtons(session, "parallel",  label= "parallel processing (not recommended)",  choices= c("yes"="yes","no"="no"),  selected = "no")
			})

			DRCFitOneS <- function(x,sel_model) {
				data_tmp <- x %>% arrange(conc)
				model <- try(drm(expr~conc, data = data_tmp, fct = fctList[[sel_model]]), silent = TRUE)
				if (class(model) == "drc") {
					parameterdf <- data.frame(coename = c("b:(Intercept)", "c:(Intercept)", "d:(Intercept)","e:(Intercept)","f:(Intercept)"), parameter = c("Slope", "Lower Limit", "Upper Limit", "ED50","f") )
					coedf <- data.frame(coename = names(model$coefficients), val = unname(model$coefficients)) %>% left_join(parameterdf, by = "coename") %>% dplyr::select(one_of("parameter","val"))
					rSquared <- cor(model$predres[,1], model$data$expr)^2
					Pvalue <- noEffect(model)[3]
					#Pvalue <- modelFit(model)$"p value"[2]
					coedf <-  rbind(coedf, data.frame(parameter = "r.squared", val = rSquared))
					coedf <-  rbind(coedf, data.frame(parameter = "p_value", val = Pvalue))
					coedft <- setNames(data.frame(t(coedf[,-1])), coedf[,1])
					coedft <- coedft %>%
					dplyr::mutate(UniqueID = unique(data_tmp$UniqueID)) %>%
					dplyr::mutate(treatment = unique(data_tmp$treatment)) %>%
					dplyr::relocate(c(UniqueID, treatment), .before = Slope)
					return(coedft)
				} #else {
				#coedft <- data.frame("Slope" = NA, "Lower.Limit" = NA, "Upper.Limit" = NA, "ED50" = NA, "r.squared" = NA, "p_value"= NA)
				#}

			}

			MultiCoreFit <- function(dataS_ChunkByCore, DRCFitOneS, fctList, sel_model, corenum) {
				cl <- makeCluster(corenum)
				clusterEvalQ(cl, {
					library(drc)
					library(dplyr)
				})
				clusterExport(cl=cl, varlist=c("dataS_ChunkByCore","DRCFitOneS","fctList","sel_model"), envir=environment())

				out <- parallel::parLapply(cl, 1:length(dataS_ChunkByCore), function (k) {
					dataStmp <- dataS_ChunkByCore[[k]]
					Res <- lapply(dataStmp, DRCFitOneS, sel_model)
					coedf <- purrr::map_df(Res, ~as.data.frame(.x), .id="id")
					return(coedf)
				})

				stopCluster(cl)
				return (out)
			}

			fit_all <- eventReactive(input$fitallButton1, {
				pcutoff <- input$pvalcut1
				datapoint <- input$datapoint1
				psel <- input$psel1
				parallel <- input$parallel
				corenum <- input$core
				sel_model <- input$sel_model2

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

					upperlimit2 = input$upper
					if (upperlimit2 != 1) {
						fctList <- list(LL.2(upper = upperlimit), LL.3(), LL.3u(upper = upperlimit), LL.4(), LL.5(),
							W1.2(upper = upperlimit), W1.3(), W1.4(), W2.2(), W2.3(), W2.4()
						)

						names(fctList) <- c("LL.2:Log-logistic(lower=0; upper=1)"="LL.2",
							"LL.3:Log-logistic(lower=0)"="LL.3",
							"LL.3u:Log-logistic(upper=1)"="LL.3u",
							"LL.4:Log-logistic"="LL.4",
							"LL.5:Generalized log-logistic"="LL.5",
							"W1.2:Weibull(type 1; lower=0; upper=1)"="W1.2",
							"W1.3:Weibull(type 1; lower=0)"="W1.3",
							"W1.4:Weibull(type 1)"="W1.4",
							"W2.2:Weibull(type 2; lower=0; upper=1"="W2.2",
							"W2.3:Weibull(type 2; lower=0)"="W2.3",
							"W2.4:Weibull(type 2)"="W2.4"
						)
					}

					if (parallel == "yes")  {
						dataS_ChunkByCore <-  split(dataS, cut(seq_along(dataS), corenum, labels = FALSE))
						out <- MultiCoreFit(dataS_ChunkByCore, DRCFitOneS, fctList, sel_model= sel_model, corenum)
						results <- ldply(out, data.frame) %>%
						#dplyr::filter(!is.na(ED50)) %>%
						dplyr::mutate_if(is.numeric, round, digits = 6)
					} else {
						Res <- lapply(dataS, DRCFitOneS, sel_model)
						results  <- purrr::map_df(Res, ~as.data.frame(.x), .id="id") %>%
						#dplyr::filter(!is.na(ED50))  %>%
						dplyr::mutate_if(is.numeric, round, digits = 6)
					}
					results$padjust <- p.adjust(results$p_value, "BH",  n = length(results$p_value))
				}
				return(results)
			})

			output$results <- DT::renderDataTable({
				results <- as.data.frame("No Fitting Results")

				if (input$fitallButton1[1] == 0) {
					if (!is.null(DataInSets[[working_project()]]$results_drc))
					results <- DataInSets[[working_project()]]$results_drc
				} else {
					withProgress(message = 'Caculating...',  detail = 'This may take a while...',  {
						results <- fit_all()
						DataInSets[[working_project()]]$results_drc <-  results
					})
				}
				if (input$saveproject1 == 1) {
					shiny::validate(need(input$projectname1!= "","Please provide project name."))
					filename <- paste("data/",input$projectname1,".RData",sep="")
					save(data_long, results, file=filename )
				}

				results <- results %>%
				dplyr::rename(any_of(lookupDRC)) %>%
				dplyr::mutate(across(c("Slope", "Lower.Limit", "Upper.Limit", "ED50", "r.squared"), round, 3)) #%>%
				#dplyr::mutate(across(c("p_value", "padjust"), round, 6))

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
					need(DataInSets[[working_project()]]$results_drc, "Need fitting results")
				)
				validate(
					need(DataInSets[[working_project()]]$data_long, "Need data")
				)

				results_drc <- DataInSets[[working_project()]]$results_drc
				data_long <- DataInSets[[working_project()]]$data_long

				if (("group" %in% colnames(data_long)) && !("treatment" %in% colnames(data_long))) {
					data_long <- data_long %>% dplyr::rename("treatment"="group") %>%
					dplyr::rename("expr" = "response")
				}

				sel_treatment <- input$sel_treatment1
				sel_model <- input$sel_model
				logbase = as.numeric(input$logbase)
				labelfontsize <- as.numeric(input$labelfontsize)
				basefontsize <- as.numeric(input$basefontsize)
				xlabel <- input$xlabel
				ylabel <- input$ylabel
				orderby <- input$orderby
				genelabel <- input$sel_geneid

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

				if (!is.null(results_drc)) {

					if (field == "r.squared") {
						sel_gene <- results_drc %>%
						dplyr::group_by(UniqueID) %>%
						dplyr::slice(which.max(!!as.symbol(field))) %>%
						dplyr::ungroup() %>%
						dplyr::arrange(desc(!!as.symbol(field))) %>%
						dplyr::slice(startslice:endslice) %>%
						dplyr::pull(UniqueID)
					} else {
						sel_gene <- results_drc %>%
						dplyr::group_by(UniqueID) %>%
						dplyr::slice(which.min(!!as.symbol(field))) %>%
						dplyr::ungroup() %>%
						dplyr::arrange(!!as.symbol(field)) %>%
						dplyr::slice(startslice:endslice) %>%
						dplyr::pull(UniqueID)
					}
				} else {
					sel_gene <- data_long %>%
					dplyr::distinct(UniqueID) %>%
					dplyr::slice(startslice:endslice) %>%
					dplyr::pull(UniqueID)
				}

				if (input$subset == "Upload Genes") {
					req(input$uploadlist)
					gene_list <- input$uploadlist
					gene_list <- ProcessUploadGeneList(gene_list)
					validate(need(length(gene_list) > 0, message = "Please input at least 1 valid gene."))
					sel_gene = gene_list
				}

				data_long_tmp  <- dplyr::filter(data_long, UniqueID %in% sel_gene) %>%
				dplyr::filter(treatment %in% sel_treatment) %>% as.data.frame()
				data_long_tmp$labelgeneid = data_long_tmp[, match(genelabel,colnames(data_long_tmp))]
				data_long_tmp$treatment = factor(data_long_tmp$treatment, levels = sel_treatment)

				data_long_tmp  <- data_long_tmp %>%
				dplyr::select(labelgeneid, treatment, conc, expr) %>%
				as.data.frame()

				df.summary <- data_long_tmp %>%
				group_by(labelgeneid, conc, treatment) %>%
				dplyr::summarise(sd = sd(expr), expr = mean(expr),  n = n(), se = sd / sqrt(n), .groups = "drop") %>%
				dplyr::select(-n)

				sel_treatment = intersect(sel_treatment, unique(data_long_tmp[['treatment']]))
				coedftlist1 <- list()
				predictedlist1 <- list()
				for (gene in unique(data_long_tmp[['labelgeneid']])) {
					dfgene1 <- data_long_tmp %>% as.data.frame() %>% dplyr::filter(labelgeneid == gene)

					predictedlist2 <- list()
					coedftlist2 <- list()
					for (onetreatment in unique(dfgene1$treatment)) {
						df <- dfgene1  %>% as.data.frame() %>% dplyr::filter(treatment == onetreatment)
						model <- try(drm(expr~conc, data = df, fct = fctList[[sel_model]]), silent = TRUE)
						if (class(model) == "drc") {
							predicted <- modelr::add_predictions(data.frame(Dose = seq(0,max(df$conc))), model)
							predictedlist2[[onetreatment]] <- predicted  %>%
							dplyr::mutate(labelgeneid = gene) %>%
							dplyr::mutate(treatment = onetreatment)

							parameterdf <- data.frame(coename = c("b:(Intercept)", "c:(Intercept)", "d:(Intercept)","e:(Intercept)","f:(Intercept)"), parameter = c("Slope", "Lower Limit", "Upper Limit", "ED50","f") )
							coedf <- data.frame(coename = names(model$coefficients), val = unname(model$coefficients)) %>% left_join(parameterdf, by = "coename") %>% dplyr::select(one_of("parameter","val"))
							rSquared <- cor(model$predres[,1], model$data$expr)^2
							Pvalue <- noEffect(model)[3]

							coedf <-  rbind(coedf, data.frame(parameter = "r.squared", val = rSquared))
							coedf <-  rbind(coedf, data.frame(parameter = "p_value", val = Pvalue))
							coedft <- setNames(data.frame(t(coedf[,-1])), coedf[,1])

							coedftlist2[[onetreatment]]  <- coedft %>%
							dplyr::mutate(labelgeneid = gene) %>%
							dplyr::mutate(treatment = onetreatment) %>%
							dplyr::relocate(c(labelgeneid, treatment), .before = Slope)

						}
					}
					predicteddf2 <- do.call(rbind,predictedlist2)
					coedftdf2 <- do.call(rbind,coedftlist2)
					predictedlist1[[gene]] <- predicteddf2
					coedftlist1[[gene]] <- coedftdf2
				}

				results  <- do.call(rbind,coedftlist1)

				predicteddf1 <- do.call(rbind, predictedlist1) %>%
				dplyr::select(c(labelgeneid, treatment, Dose, pred)) %>%
				tibble::remove_rownames()  %>%
				dplyr::rename(conc = Dose,  expr = pred)


				gg.df <- rbind(cbind(geom="pt", data_long_tmp), cbind(geom="ln", predicteddf1))%>%
				dplyr::rename(GroupName = treatment)

				p <- ggplot(gg.df, aes(x=conc, y=expr, color=GroupName)) +
				facet_wrap(~ labelgeneid, scales = "free", nrow = nrow, ncol = ncol)

				p <- p +
				geom_point(data=gg.df[gg.df$geom=="pt",], shape=4)+
				geom_line(data=gg.df[gg.df$geom=="ln",])

				p <- p +
				geom_errorbar(aes(ymin = expr-sd, ymax = expr+sd, color = treatment), data = df.summary, width = 0.2) +
				theme_bw(base_size = basefontsize) +
				xlab(xlabel) +  ylab(ylabel) +
				theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
					axis.text.x = element_text(angle = 0),
					legend.title = element_blank(),
					strip.text.x = element_text(size=input$labelfontsize),
				legend.position="bottom")

				if (logbase != 1){
					p <- p +  scale_x_continuous(trans=scales::pseudo_log_trans(base = logbase))
				}

				return(list(plot=p, result = results))
			})

			output$browsing <- renderPlot({
				browsing_out()[["plot"]]
			})

			output$browsing_result  <- DT::renderDataTable({
				fitresult <- browsing_out()[["result"]]
				fitresult <- fitresult %>%
				dplyr::rename(any_of(lookupDRC)) %>%
				dplyr::mutate(across(c("Slope", "Lower.Limit", "Upper.Limit", "ED50", "r.squared"), round, 3)) #%>%
				#dplyr::mutate(across(c("p_value"), round, 6))
				DT::datatable(fitresult, options = list(pageLength = 10)) 
			})

			output$download_results_button <- shiny::downloadHandler(
				filename = function() {
					paste("Results-", Sys.Date(), ".csv", sep="")
				},
				content = function(file) {
					results <- DataInSets[[working_project()]]$results_drc
					results[,sapply(results,is.numeric)] <- signif(results[,sapply(results,is.numeric)],3)
					write.csv(results, file)
				}
			)
		}
	)
}
