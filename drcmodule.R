###########################################################################################################
##Omics Dose-Response/Time-Course Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: drcmodule.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 12/16/2021
##@version 2.0
###########################################################################################################


##########################################################################################################
## Curve Fitting Plot
##########################################################################################################
#pkgs: "drc", "ggpmisc", "parallel", "dplyr", "scales", "DT", "purrr", "modelr", "tibble"  

library(drc)
library(ggpmisc)
library(parallel)
library(gridExtra)

corenumber  = parallel::detectCores(logical = FALSE)

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
fctList2 <- fctList

plotfitting <- function(model, group, TimeDose="conc",  npcx="auto", npcy="auto", label = "ED50",xlabel="Log10(conc)", ylabel="Response", basefontsize=14, logbase=10) {
	parameterdf <- data.frame(coename = c("b:(Intercept)", "c:(Intercept)", "d:(Intercept)","e:(Intercept)","f:(Intercept)"), parameter = c("Slope", "Lower Limit", "Upper Limit", "ED50","f") )
	coedf <- data.frame(coename = names(model$coefficients), val = unname(model$coefficients)) %>% left_join(parameterdf, by = "coename") %>% dplyr::select(one_of("parameter","val"))

	rSquared <- cor(model$predres[,1], model$data$response)^2
	#Pvalue <- noEffect(model)[3]
	#Pvalue <- modelFit(model)$"p value"[2]
	coedf <-  rbind(coedf, data.frame(parameter ="r2",val = rSquared)) %>% mutate_if(is.numeric, round, digits = 2)
	coedft <-  setNames(data.frame(t(coedf[,-1])), coedf[,1])

	fitDAT <- model$origData
	xvalue <- model$dataList$dose
	newdata <- data.frame(DOSE = 0:max(xvalue), CURVE = rep(1,length(0:max(xvalue))))
	predicted <- data.frame(x = c(0:max(xvalue)), predicted= predict(model,newdata))

	p <- ggplot(fitDAT, aes(x=conc, y=response)) +
	geom_point() +
	ylim(0, max(fitDAT['response'])) +
	geom_line(color='red', data = predicted, aes(x=x, y=predicted)) +
	#geom_hline(yintercept = 0.5, linetype="dashed",  color = "blue", size=1) +
	ggtitle(group) +
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
		column(2,
			wellPanel(
				conditionalPanel(ns = ns,"input.expression_tabset=='Model Selection' || input.expression_tabset=='Fitting Curve' || input.expression_tabset=='Data Table' || input.expression_tabset=='Result Table' ",
					selectizeInput(ns("sel_gene1"),	label="Gene Name",	choices = NULL,	multiple=FALSE, options = list(placeholder =	'Type to search'))
				),
				conditionalPanel(ns = ns,"input.expression_tabset=='Fitting Curve' || input.expression_tabset=='Browsing' || input.expression_tabset=='Data Table' || input.expression_tabset=='Result Table' ",
					selectizeInput(ns("sel_treatment1"),	label="Treatment",	choices = NULL,	multiple=TRUE),
				),
				conditionalPanel(ns = ns,"input.expression_tabset=='Fitting Curve' || input.expression_tabset=='Browsing'",
					selectizeInput(ns("sel_model"),	label="Model",	choices = NULL,	multiple=FALSE),
					numericInput(ns("upper"), "Upper:", 1, min = 1, max = 100),
				),
				conditionalPanel(ns = ns,"input.expression_tabset=='Fitting Curve'",
					radioButtons(ns("separateplot"), label="Plot Separately ", inline = TRUE, choices = c("Yes" = "Yes","No" = "No"))
				),
				conditionalPanel(ns = ns,"input.expression_tabset=='Model Selection'",
					selectizeInput(ns("sel_treatment1b"),	label="Treatment",	choices = NULL,	multiple=FALSE),
					sliderInput(ns("sel_topn"), "Top Model Number:", min = 1, max = 11, step = 1, value = 5)
				),
				conditionalPanel(ns = ns,"input.expression_tabset=='Browsing'",
					radioButtons(ns("orderby1"), label="Sort by", inline = TRUE, choices =  c("p_value" = "p_value", "padjust" = "padjust","r.squared" = "r.squared")),
					column(width=6,selectInput(ns("sel_page1"),	label="Select Page",	choices = NULL,	selected=1)),
					column(width=6,selectInput(ns("numperpage1"), label= "Plot Number per Page", choices= c("4"=4,"6"=6,"9"=9), selected=6))
				),
				conditionalPanel(ns = ns,"input.expression_tabset=='Model Selection' || input.expression_tabset=='Fitting Curve' || input.expression_tabset=='Browsing'",
					radioButtons(ns("npcx"), label="Label X Position", inline = TRUE, choices =  c("auto" = "auto", "left" = "left","middle"="middle","right"="right")),
					radioButtons(ns("npcy"), label="Label Y Position", inline = TRUE, choices =  c("auto" = "auto", "top" = "top","center" = "center","bottom" = "bottom")),
					radioButtons(ns("logbase"), label="Log tansform", inline = TRUE,  choices=c("Non"= 1, "10"= 10), selected=10),
					textInput(ns("xlabel"), "xlabel", value = "Log10(Conc.)"),
					textInput(ns("ylabel"), "ylabel", value = "Response"),
					sliderInput(ns("labelfontsize"), "Label Font Size:", min = 10, max = 24, step = 2, value = 14),
					sliderInput(ns("basefontsize"), "Font Size:", min = 10, max = 24, step = 2, value = 14),
					sliderInput(ns("refreshrate"), "Plot Refresh Rate:", min = 0, max = 8000, step = 1000, value = 1000)
				),
				conditionalPanel(ns = ns,"input.expression_tabset=='Result Table (IC50/EC50)'",
					radioButtons(ns("parallel1"), label= "parallel processing?", choices= c("yes"="yes","no"="no"),inline = TRUE),
					numericInput(ns("core1"), "Core will be used:", value = 4, min = 2, max = 12),
					radioButtons(ns("psel1"), label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
					numericInput(ns("pvalcut1"), label= "Choose P-value Threshold",  value=0.01, min=0, step=0.001),
					numericInput(ns("datapoint1"), label= "Minimal Data Points",  value=5, min=5, step=1),
					selectizeInput(ns("sel_model2"),	label="Model",	choices = NULL,	multiple=FALSE),
					numericInput(ns("upper1"), "Upper:", 1, min = 1, max = 100),
					uiOutput(ns("filteredgene1")),
					radioButtons(ns("saveproject1"), label="Save Result?", inline = TRUE,  choices=c("Yes"= 1, "No"= 0), selected=0),
					textInput(ns("projectname1"), "Project Name"),
					actionButton(ns("fitallButton1"), "Fit filtered data")
				),
			)
		),
		column(10,
			tabsetPanel(id=ns("expression_tabset"),
				tabPanel(title="Fitting Curve",
					#actionButton(ns("FittingCurve"), "Save to output"),
					fluidRow(
						plotOutput(ns("FittingCurve"), height=800)
					),
					fluidRow(
						DT::dataTableOutput(ns("fitresult"))
					)
				),
				tabPanel(title="Model Selection",
					#actionButton(ns("ModelSelection"), "Save to output"),
					fluidRow(
						plotOutput(ns("ModelSelection"), height=800)
					),
					fluidRow(
						DT::dataTableOutput(ns("fitresultb"))
					)
				),
				tabPanel(title="Data Table",	DT::dataTableOutput(ns("fitdata"))),
				tabPanel(title="Result Table (IC50/EC50)",
					#actionButton("results", "Save to output"),
					#shinycssloaders::withSpinner(dataTableOutput(ns("results1")),type = 4, size = 2,color = "#0000FF")
					dataTableOutput(ns("results1"))
				),
				tabPanel(title="Browsing",
					fluidRow(
						plotOutput(ns("browsing1"), height=800)
					),
					fluidRow(
						DT::dataTableOutput(ns("browsing1_result"))
					)
				),
				tabPanel(title="Help", htmlOutput('help_FittingCurveDRC'))
			)
		)
	)
}

drc_server <- function(id) {
	shiny::moduleServer(id,
		function(input, output, session) {
			ns <- shiny::NS(id)
			observe({
				req(length(working_project()) > 0)
				#req(DataInSets[[working_project()]]$data_long)
				data_long <- DataInSets[[working_project()]]$data_long
				req("UniqueID" %in% colnames(data_long) & "group" %in% colnames(data_long))
				DataIngenes <-  data_long  %>% dplyr::pull(UniqueID) %>% unique() %>% as.character()
					updateSelectizeInput(session,'sel_gene1', choices= DataIngenes, server=TRUE)
				group <-  data_long %>% dplyr::pull(group) %>% unique() %>% as.character()
				updateSelectizeInput(session,'sel_treatment1', choices= group,  selected=group)
				updateSelectizeInput(session,'sel_treatment1b', choices= group,  selected=group[1])
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

			DataExpReactive <- reactive({
				req(length(working_project()) > 0)
				data_long <- DataInSets[[working_project()]]$data_long
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

				sel_treatment = intersect(sel_treatment, unique(data_tmp[['group']]))
				#sel_treatment = unique(data_tmp[['group']])
				out <- list()
				plist <- list()
				#models <- list()
				fitDATlist <- list()
				predictedlist <- list()

				for (onegroup in sel_treatment) {
					dfgene1 <- data_tmp %>% as.data.frame() %>% dplyr::filter(group == onegroup)
					model <- try(drm(response~conc, data = dfgene1, fct = fctList[[sel_model]]), silent = TRUE)

					if (class(model) == "drc") {
						#models[[onegroup]] <- model
						fitDATlist[[onegroup]] = model$origData
						xvalue <- model$dataList$dose
						newdata <- data.frame(DOSE = 0:max(xvalue), CURVE = rep(1,length(0:max(xvalue))))
						predicted <- data.frame(x = c(0:max(xvalue)), predicted= predict(model,newdata))
						######
						#predslm = predict(model, newdata, interval = "confidence")
						#######
						predictedlist[[onegroup]] <- predicted
						fitting_PR <- plotfitting(model,group=onegroup, TimeDose="conc", label = "ED50", npcx, npcy, xlabel, ylabel, basefontsize, logbase)
						p <- fitting_PR$plot
						result <- fitting_PR$result
						#rownames(result) <- paste(sel_gene, group,sep="_")
						out[[onegroup]] <- result
					}	else {
						p <- ggplot(dfgene1, aes(x=conc, y=response)) +
						geom_point() +
						ylim(0, max(dfgene1['response'])) +
						ggtitle(onegroup) +
						theme_bw(base_size = basefontsize) + xlab(xlabel) +  ylab(ylabel) +
						theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), axis.text.x = element_text(angle = 0),legend.title = element_blank(), legend.position="bottom")
						if (logbase != 1){
							p <- p +  scale_x_continuous(trans=scales::pseudo_log_trans(base = logbase))
						}
					}
					plist[[onegroup]]  <- p
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
					predicteddf <- cbind(group=rep(names(predictedlist), sapply(predictedlist,nrow)),do.call(rbind,predictedlist))
					fitDATdf <- do.call(rbind, fitDATlist)

					p <- ggplot(fitDATdf, aes(x=conc, y=response, color=group)) +
					geom_point() +
					ylim(0, max(fitDATdf['response'])) +
					geom_line(data = predicteddf, aes(x=x, y=predicted, color=group)) +
					#geom_ribbon(data = predicteddf)+
					#ggtitle(group) +
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
				topn = input$sel_topn
				upperlimit = input$upper

				DataExpReactive_tmp <- DataExpReactive()
				data_tmp <- DataExpReactive_tmp[["data_tmp"]]

				if (is.null(data_tmp$group)) {
					dfgene1 <- data_tmp
				} else {
					sel_treatment2 = intersect(sel_treatment2, unique(data_tmp[['group']]))
					dfgene1 <- data_tmp %>% dplyr::filter(group == sel_treatment2)
				}

				model <- try(drm(response~conc, data = dfgene1, fct = LL.4()), silent = TRUE)

				if (!inherits(model, "try-error")){
					if (npcy == "auto")
					npcy = "top"

					if (npcx == "auto") {
						if (model$coefficients["b:(Intercept)"] > 0)
						npcx = "right"
						else
						npcx = "left"
					}

					fctList2[['LL.4']] <- NULL

					sel <- try(drc::mselect(model, fctList2, sorted = "IC", linreg = FALSE), silent = TRUE)
					sel <- sel %>% as.data.frame() %>% top_n(., -topn, IC) %>% mutate_if(is.numeric, round, digits = 6)

					fitDAT = model$origData
					xvalue <- model$dataList$dose
					randomdose <- 0:max(xvalue)
					newdata <- data.frame(DOSE = randomdose, CURVE = rep(1,length(randomdose)))
					predictedlist <- list()
					for (modelname in dimnames(sel)[[1]]){
						model <- try(drm(response~conc, data = dfgene1, fct = fctList[[modelname]]), silent = TRUE)
						if (!inherits(model, "try-error")) {
							predicted <- data.frame(modelname=rep(modelname,length(randomdose)), x = randomdose, predicted= predict(model,newdata))
							predictedlist[[modelname]] <- predicted
						}
					}
					predicteddf <- do.call(rbind, predictedlist)

					p <- ggplot(fitDAT, aes(x=conc, y=response)) +
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
					p <- ggplot(dfgene1, aes(x=conc, y=response)) +
					geom_point() +
					ylim(0, max(dfgene1['response'])) +
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

			###########################################################################################################
			#DRC fit all data

			observe({
				req(length(working_project()) > 0)

				updateNumericInput(session, 'core1', label = "Core will be used:", value = corenumber,  min = 2, max = corenumber, step =1)
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
					UniqueIDnum <- length(unique(data_long$UniqueID))
				}
				output$filteredgene1 =	 renderText({paste("<font color=\'red\'><b>Total Genes: ", UniqueIDnum, "</b></font>",sep="")})
				if (UniqueIDnum  < 200)
				updateRadioButtons(session, "parallel1",  label= "parallel processing (not recommended)",  choices= c("yes"="yes","no"="no"),  selected = "no")
			})


			DRCFitOneS <- function(x,sel_model) {
				data_tmp <- x %>% arrange(conc)
				model <- try(drm(response~conc, data = data_tmp, fct = fctList[[sel_model]]), silent = TRUE)
				if (class(model) == "drc") {
					parameterdf <- data.frame(coename = c("b:(Intercept)", "c:(Intercept)", "d:(Intercept)","e:(Intercept)","f:(Intercept)"), parameter = c("Slope", "Lower Limit", "Upper Limit", "ED50","f") )
					coedf <- data.frame(coename = names(model$coefficients), val = unname(model$coefficients)) %>% left_join(parameterdf, by = "coename") %>% dplyr::select(one_of("parameter","val"))
					rSquared <- cor(model$predres[,1], model$data$response)^2
					Pvalue <- noEffect(model)[3]
					#Pvalue <- modelFit(model)$"p value"[2]
					coedf <-  rbind(coedf, data.frame(parameter = "r.squared", val = rSquared))
					coedf <-  rbind(coedf, data.frame(parameter = "p_value", val = Pvalue))
					coedft <- setNames(data.frame(t(coedf[,-1])), coedf[,1])
					coedft <- coedft %>%
					dplyr::mutate(UniqueID = unique(data_tmp$UniqueID)) %>%
					dplyr::mutate(group = unique(data_tmp$group)) %>%
					dplyr::relocate(c(UniqueID, group), .before = Slope)
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
				}
			)
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
			parallel <- input$parallel1
			corenum <- input$core1
			sel_model <- input$sel_model2

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
					dataS <- named_group_split(data_long_sub,UniqueID,group)
				} else {
					dataS <- named_group_split(data_long,UniqueID,group)
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

		output$results1 <- DT::renderDataTable({
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

			results <- results %>%  mutate_if(is.numeric, round, digits = 2)
			DT::datatable(results, options = list(pageLength = 15), rownames= FALSE)
		})

		###########################################################################################################
		#browsing
		observe({
			updateSelectInput(session,'sel_page1', choices= seq_len(100))
		})

		browsing_out_drc <- reactive({
			req(DataInSets[[working_project()]]$data_long)
			results_drc <- DataInSets[[working_project()]]$results_drc
			data_long <- DataInSets[[working_project()]]$data_long
			sel_treatment <- input$sel_treatment1
			sel_model <- input$sel_model
			logbase = as.numeric(input$logbase)
			labelfontsize <- input$labelfontsize
			basefontsize <- input$basefontsize
			xlabel <- input$xlabel
			ylabel <- input$ylabel
			numperpage <- as.numeric(input$numperpage1)
			sel_page <- as.numeric(input$sel_page1)-1
			orderby <- input$orderby1

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

			if (!is.null(results_drc)) {
				sel_gene <- results_drc %>%
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

			#sel_gene <- unique(data_long$UniqueID)[1:10]
			#sel_model  <- "LL.3"
			data_long_tmp  <- data_long %>%
			dplyr::filter((UniqueID %in% sel_gene) & (group %in% sel_treatment)) %>%
			dplyr::select(UniqueID, group, conc, response) %>%
			as.data.frame()

			if(numperpage==4) {
				nrow = 2; ncol = 2
			} else if(numperpage==6) {
				nrow = 2; ncol = 3
			} else {
				nrow = 3; ncol = 3
			}

			sel_treatment = intersect(sel_treatment, unique(data_long_tmp[['group']]))
			coedftlist1 <- list()
			predictedlist1 <- list()
			for (gene in  sel_gene) {
				dfgene1 <- data_long_tmp %>% as.data.frame() %>% dplyr::filter(UniqueID == gene)
				predictedlist2 <- list()
				coedftlist2 <- list()
				for (onegroup in unique(dfgene1$group)) {
					df <- dfgene1  %>% as.data.frame() %>% dplyr::filter(group == onegroup)
					model <- try(drm(response~conc, data = df, fct = fctList[[sel_model]]), silent = TRUE)
					if (class(model) == "drc") {
						predicted <- modelr::add_predictions(data.frame(Dose = seq(0,max(df$conc))), model)
						predictedlist2[[onegroup]] <- predicted  %>%
						dplyr::mutate(UniqueID = gene) %>%
						dplyr::mutate(group = onegroup)

						parameterdf <- data.frame(coename = c("b:(Intercept)", "c:(Intercept)", "d:(Intercept)","e:(Intercept)","f:(Intercept)"), parameter = c("Slope", "Lower Limit", "Upper Limit", "ED50","f") )
						coedf <- data.frame(coename = names(model$coefficients), val = unname(model$coefficients)) %>% left_join(parameterdf, by = "coename") %>% dplyr::select(one_of("parameter","val"))
						rSquared <- cor(model$predres[,1], model$data$response)^2
						Pvalue <- noEffect(model)[3]

						coedf <-  rbind(coedf, data.frame(parameter = "r.squared", val = rSquared))
						coedf <-  rbind(coedf, data.frame(parameter = "p_value", val = Pvalue))
						coedft <- setNames(data.frame(t(coedf[,-1])), coedf[,1])

						coedftlist2[[onegroup]]  <- coedft %>%
						dplyr::mutate(UniqueID = gene) %>%
						dplyr::mutate(group = onegroup) %>%
						dplyr::relocate(c(UniqueID, group), .before = Slope)

					} #else {
					#coedft <- data.frame("Slope" = NA, "Lower.Limit" = NA, "Upper.Limit" = NA, "ED50" = NA, "r.squared" = NA, "p_value"= NA)
					#}
				}
				predicteddf2 <- do.call(rbind,predictedlist2)
				coedftdf2 <- do.call(rbind,coedftlist2)
				predictedlist1[[gene]] <- predicteddf2
				coedftlist1[[gene]] <- coedftdf2
			}

			results  <- do.call(rbind,coedftlist1)

			predicteddf1 <- do.call(rbind,predictedlist1) %>%
			dplyr::select(c( UniqueID, group, Dose, pred)) %>%
			tibble::remove_rownames()  %>%
			dplyr::rename(conc = Dose,  response = pred)

			gg.df <- rbind(cbind(geom="pt", data_long_tmp),cbind(geom="ln",predicteddf1))%>%
			dplyr::rename(GroupName = group)

			p <- ggplot(gg.df, aes(x=conc, y=response, color=GroupName)) +
			geom_point(data=gg.df[gg.df$geom=="pt",], shape=4)+
			geom_line(data=gg.df[gg.df$geom=="ln",])+
			facet_wrap(~ UniqueID, scales = "free", nrow = nrow, ncol = ncol)	+
			theme_bw(base_size = basefontsize) + xlab(xlabel) +  ylab(ylabel) +
			theme (plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), axis.text.x = element_text(angle = 0),legend.title = element_blank(), legend.position="bottom")


			if (logbase != 1){
				p <- p +  scale_x_continuous(trans=scales::pseudo_log_trans(base = logbase))
			}

			return(list(plot=p, result = results))
		})

		output$browsing1 <- renderPlot({
			browsing_out_drc()[["plot"]]
		})

		output$browsing1_result  <- DT::renderDataTable({
			fitresult <- browsing_out_drc()[["result"]]
			fitresult <- fitresult %>%
			mutate_if(is.numeric, round, digits = 4)
			DT::datatable(fitresult, options = list(pageLength = 10))
		})

	}
)
}
