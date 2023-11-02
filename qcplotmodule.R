###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: qcplotmodule.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 1/6/2022
##@version 3.0
###########################################################################################################
# pkgs: "rgl","car","factoextra","ComplexHeatmap", "DT", "dplyr", "tibble", "ggpubr"
# req data:MetaData, ProteinGeneNameHeader, exp_unit, data_long, data_wide, group_order, sample_order
library(rgl)
library(car)
library(factoextra)
library(ComplexHeatmap)
library(cowplot) #a simple add-on to ggplot

qcplot_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
		column(3,
			wellPanel(
				uiOutput(ns('loadedprojects')),
				tags$style(mycss),
				uiOutput(ns("selectGroupSample")),

				conditionalPanel(ns = ns, "input.tabset=='PCA Plot'",
					numericInput(ns("MaxPCANum"), label= "maximal number of principal components", value=10, min=2, max = 20, step=1),
					conditionalPanel(ns = ns, "input.PCA_tabset=='PCA Plot' || input.PCA_tabset=='PCA 3D Interactive' || input.PCA_tabset=='PCA 3D Plot'",
						selectInput(ns("PCAcolorby"), label="Color (Group) By", choices=NULL)
					),
					conditionalPanel(ns = ns, "input.PCA_tabset=='PCA Plot' || input.PCA_tabset=='PCA 3D Interactive'",
						selectInput(ns("PCAshapeby"), label="Shape By", choices=NULL)
					),
					conditionalPanel(ns = ns, "input.PCA_tabset=='PCA Plot'",
						#selectInput(ns("PCAsizeby"), label="Size By", choices=NULL),
						selectizeInput(ns("pcnum"),label="Select Principal Components", choices=1:10, multiple=TRUE, selected=1:2, options = list(maxItems = 2)),
						radioButtons(ns("ellipsoid"), label="Plot Ellipsoid (>3 per Group)", inline = TRUE, choices = c("No" = FALSE,"Yes" = TRUE)),
						radioButtons(ns("mean_point"), label="Show Mean Point", inline = TRUE, choices = c("No" = FALSE, "Yes" = TRUE)),
						radioButtons(ns("loading"), label="Show PCA Loading", inline = TRUE, choices = c("No" = "No", "Yes" = "Yes")),
						conditionalPanel(ns = ns,"input.loading == 'Yes'",
							radioButtons(ns("genelabel"), label="Select Gene Label", inline = TRUE, choices="Gene.Name"),
							numericInput(ns("contrib_factors"), "Select the number of the top contributing factors",10, min = 5, max = 100)
						),
						radioButtons(ns("convex"), label="Show Convex Polygon", inline = TRUE, choices = c("No" = "No","Yes" =  "Yes")),
						radioButtons(ns("rug"), label="Show Marginal Rugs", inline = TRUE, choices = c("No" = FALSE,"Yes" = TRUE)),
						fluidRow(
							column(width=6, sliderInput(ns("PCAdotsize"), "Dot Size:", min = 1, max = 20, step = 1, value = 4)),
							column(width=6, sliderInput(ns("PCAfontsize"), "Label Font Size:", min = 1, max = 20, step = 1, value = 10))
						),
						radioButtons(ns("PCA_subsample"), label="Label Samples:", inline = TRUE, choices = c("All","None", "Subset"), selected = "All"),
						conditionalPanel("input.PCA_subsample!='None'",
							radioButtons(ns("PCA_label"),label="Select Sample Label", inline = TRUE, choices="")
						),
						conditionalPanel(ns = ns, "input.PCA_subsample=='Subset'",
							actionButton(ns("PCA_refresh_sample"), "Reload Sample IDs"),
							textAreaInput(ns("PCA_list"), "List of Samples to Label", "", cols = 5, rows=6)
						)
					),
					conditionalPanel(ns = ns, "input.PCA_tabset=='PCA 3D Plot'",
						radioButtons(ns("ellipsoid3d"), label="Plot Ellipsoid (>3 per Group)", inline = TRUE, choices = c("No" = "No","Yes" = "Yes")),
						radioButtons(ns("dotlabel"), label="Dot Label", inline = TRUE, choices =  c("No" = "No","Yes" = "Yes"))
					),
					conditionalPanel(ns = ns, "input.PCA_tabset=='PCA Loading Table'",
						radioButtons(ns("abs"), label="Absolute value", inline = TRUE, choices = c("No" = "No", "Yes" = "Yes"))
					)
				),
				conditionalPanel(ns = ns, "input.tabset=='PCA Plot'  || input.tabset=='Box Plot' || input.tabset=='CV Distribution' || input.tabset=='Dendrograms'",
					selectInput(ns("colpalette"), label= "Select palette", choices="")
				),
				conditionalPanel(ns = ns, "input.tabset=='Covariates'",
					selectizeInput(ns("covar_variates"), label="Select Covariates:", choices=NULL, multiple = TRUE),
					numericInput(ns("covar_PC_cutoff"), label= "Principle Component (PC) Cutoff (% explained variance)",  value = 5, min=0, max=100, step=1),
					numericInput(ns("covar_FDR_cutoff"), label= "Choose FDR Value Cutoff", value=0.1, min=0, max=1, step=0.001),
					sliderInput(ns("covar_ncol"), label= "Number of Columns for Plots", min = 1, max = 6, step = 1, value = 3)
				),
				conditionalPanel(ns = ns, "input.tabset=='Dendrograms'",
					fluidRow(
						column(width=6, sliderInput(ns("DendroCut"), label="tree cut number:", min = 2, max = 10, step = 1, value = 4)),
						column(width=6, sliderInput(ns("DendroFont"), label= "Label Font Size:", min = 0.5, max = 4, step = 0.5, value = 1))
					),
					radioButtons(ns("dendroformat"), label="Select Plot Format", inline = TRUE, choices = c("tree" = "tree","horizontal" = "horiz", "circular" = "circular"), selected="circular")
				),
				conditionalPanel(ns = ns, "input.tabset=='Box Plot'",
					sliderInput(ns("axisfontsize"), "Axis Font Size:", min = 10, max = 28, step = 1, value = 16),
					textInput(ns("Ylab"), "Y label", width = "100%"),
					textInput(ns("Xlab"), "X label", width = "100%"),
					sliderInput(ns("Xangle"), label= "X Angle", min = 0, max = 90, step = 5, value = 90)
				),
				conditionalPanel(ns = ns, "input.tabset=='AlignQC'",
					tags$hr(style="border-color: black;") ,
					conditionalPanel(ns = ns, "input.AlignQC_tabset=='Plot Other Variables'",
						selectizeInput(ns("alignQC_var"), label="Select Numerica Variable to Plot", choices=NULL,multiple=FALSE)
					),
					conditionalPanel(ns = ns, "input.AlignQC_tabset!='Top Gene List'",
						sliderInput(ns("alignQC_height"), "Plot Height:", min = 200, max = 3000, step = 50, value = 700),
						sliderInput(ns("alignQC_width"), "Plot Width:", min = 200, max = 3000, step = 50, value = 1000)
					),
					conditionalPanel(ns = ns, "input.AlignQC_tabset=='Top Gene List'",
						sliderInput(ns("alignQC_TG_Ngene"), "Max # of Top Genes Per Sample:", min = 1, max = 100, step = 1, value = 30),
						sliderInput(ns("alignQC_TG_Ntotal"), "Max # of Genes to Show", min = 10, max = 200, step = 5, value = 90),
						checkboxInput(ns("alignQC_TG_list"), "Show Table Instead of Graph?",  FALSE, width="90%"),
						conditionalPanel(ns = ns,"input.alignQC_TG_list==0",
							sliderInput(ns("alignQC_TG_height"), "Gene Plot Height:", min = 200, max = 3000, step = 50, value = 800),
							sliderInput(ns("alignQC_TG_width"), "Gene Plot Width:", min = 200, max = 3000, step = 50, value = 800)
						)
					),
					conditionalPanel(ns = ns, "input.AlignQC_tabset!='Top Gene List' || input.alignQC_TG_list==0",
						sliderInput(ns("alignQC_fontsize"), "Label Font Size:", min = 1, max = 18, step = 1, value = 12)
					),
					conditionalPanel(ns = ns, "input.AlignQC_tabset=='Top Gene Ratio' || input.AlignQC_tabset=='Top Gene List'",
						radioButtons(ns("convert_exp"), label="Convert Expression Data", inline = TRUE, choices = c("None", "log to linear"), selected = "log to linear"),
						conditionalPanel(ns = ns, "input.convert_exp=='log to linear'",
							fluidRow(
								column(width=6, numericInput(ns("convert_logbase"), label= "log base",  value = 2, min=1, step=1)),
								column(width=6, numericInput(ns("convert_small"), label= "small value",  value=1, min=0, step=0.1))
							)
						),
						conditionalPanel(ns = ns, "input.AlignQC_tabset!='Top Gene List' | input.alignQC_TG_list==0",
							textInput(ns("alignQC_TGR_Y"), "Label for Top Gene %", value="% of Total TPM", width = "100%")
						)
					)
				)
			)
		),
		column(9,
			tabsetPanel(id=ns("tabset"),
				tabPanel(title="PCA Plot", value="PCA Plot",
					tabsetPanel(id=ns("PCA_tabset"),
						tabPanel(title="PCA Plot", value="PCA Plot",
							#actionButton(ns("plot_PCA"), "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
							actionButton(ns("pcaplot"), "Save to output"),
							plotOutput(ns("pcaplot"),height = 800),
						),
						tabPanel(title="Eigenvalues", value="Eigenvalues",
							plotOutput(ns("Eigenvalues"), height = 800)
						),
						tabPanel(title="PCA Loading Table", value="PCA Loading Table",
							downloadButton(ns("loadingPCA_download_button"), "Download as CSV", icon("download",lib = "font-awesome")), br(),
							DT::dataTableOutput(ns("loadingPCATable"))
						),
						tabPanel(title="PCA 3D Plot", value="PCA 3D Plot",
							plotOutput(ns("pca_legend"), height = 100),
							rglwidgetOutput(ns("plot3d"),  width = 800, height = 800)
						),
						tabPanel(title="PCA 3D Interactive", value="PCA 3D Interactive",
							plotlyOutput(ns("plotly3d"),  width = 1000, height = 800)
						)
					)
				),
				tabPanel(title="Covariates", value="Covariates",
					tabsetPanel(id=ns("covartiate_tabset"),
						tabPanel(title="Summary", value="Summary",
							actionButton(ns("compute_PC"), "Compute/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
							textOutput(ns("N_pairs")),
							tags$br(),
							dataTableOutput(ns("covar_table")
							)
						),
						tabPanel(title="Categorical Covariates", value="Categorical Covariates",
							h5("After changing parameters, please click Refresh button in the Summary panel to generate new plots."),
							actionButton(ns("covar_cat"), "Save to output"),
							sliderInput(ns("covar_cat_height"), "Plot Height:", min = 200, max = 3000, step = 50, value = 800),
							tags$br(),
							uiOutput(ns("plot.PC_covariatesC")
							)
						),
						tabPanel(title="Numerical Covariates", value="Numerical Covariates",
							h5("After changing parameters, please click Refresh button in the Summary panel to generate new plots."),
							actionButton(ns("covar_num"), "Save to output"),
							sliderInput(ns("covar_num_height"), "Plot Height:", min = 200, max = 3000, step = 50, value = 800),
							tags$br(),
							uiOutput(ns("plot.PC_covariatesN")
							)
						)
					)
				),
				tabPanel(title="AlignQC", value="AlignQC",
					tabsetPanel(id=ns("AlignQC_tabset"),
						tabPanel(title="Top Gene Ratio", value="Top Gene Ratio",
							actionButton(ns("alignQC_TGR"), "Save to output"),
							uiOutput(ns("alignQC_TGR_plot")
							)
						),
						tabPanel(title="Top Gene List", value="Top Gene List",
							actionButton(ns("alignQC_TGL"), "Save to output"),
							uiOutput(ns("alignQC_TG_plot")
							)
						),
						tabPanel(title="Mapped Read Allocation", value="Mapped Read Allocation",
							actionButton(ns("alignQC_RA"), "Save to output"),
							uiOutput(ns("alignQC_RA_plot")
							)
						),
						tabPanel(title="Plot Other Variables", value="Plot Other Variables",
							actionButton(ns("alignQC_OV"), "Save to output"),
							uiOutput(ns("alignQC_OV_plot")
							)
						)
					)
				),
				tabPanel(title="Sample-sample Distance", value="Sample-sample Distance",
					actionButton(ns("SampleDistance"), "Save to output"),
					plotOutput(ns("pheatmap"),height = 800)
				),
				tabPanel(title="Dendrograms", value="Dendrograms",
					actionButton(ns("Dendrograms"), "Save to output"),
					plotOutput(ns("Dendrograms"),height = 800)
				),
				tabPanel(title="Box Plot", value="Box Plot",
					actionButton(ns("QCboxplot"), "Save to output"),
					plotOutput(ns("QCboxplot"),height = 800)
				),
				tabPanel(title="CV Distribution", value="CV Distribution",
					actionButton(ns("histplot"), "Save to output"),
					plotOutput(ns("histplot"),height = 800)
				),
				tabPanel(title="Help", value="Help",
					htmlOutput("help_QC")
				)
			)
		)
	)
}

qcplot_server <- function(id) {
	shiny::moduleServer(id,
		function(input, output, session) {
			source("PC_Covariates.R")
			ns <- shiny::NS(id)

			output$loadedprojects <- renderUI({
				req(length(working_project()) > 0)
				radioButtons(ns("current_dataset"), label = "Change Working Dataset", choices=names(DataInSets), inline = F, selected=working_project())
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
				req(DataInSets[[working_project()]]$MetaData)
				MetaData=DataInSets[[working_project()]]$MetaData
				ProteinGeneNameHeader = DataInSets[[working_project()]]$ProteinGeneNameHeader
				exp_unit <- DataInSets[[working_project()]]$exp_unit
				updateSelectizeInput(session,'colpalette', choices=rownames(brewer.pal.info), selected="Dark2")
				attributes=sort(setdiff(colnames(MetaData), c("sampleid", "Order", "ComparePairs") ))
				updateSelectInput(session, "PCAcolorby", choices=attributes, selected="group")
				updateSelectInput(session, "PCAshapeby", choices=c("none", attributes), selected="none")
				#updateSelectInput(session, "PCAsizeby", choices=c("none", attributes), selected="none")
				sampleIDs=sort(setdiff(colnames(MetaData), c("Order", "ComparePairs") ))
				updateRadioButtons(session,'PCA_label', inline = TRUE, choices=sampleIDs, selected="sampleid")
				updateSelectInput(session, "covar_variates", choices=attributes, selected=attributes)
				updateRadioButtons(session,'genelabel', inline = TRUE, choices=ProteinGeneNameHeader, selected="Gene.Name")
				updateTextInput(session, "Ylab", value=exp_unit)
				if (!is.null(MetaData)) {
					if (nrow(MetaData)>100) {updateRadioButtons(session,'PCA_subsample', selected="None")} #when there are too many samples, don't show  labels
				}
			})

			DataQCReactive <- reactive({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$MetaData)
				req(DataInSets[[working_project()]]$data_long)
				req(DataInSets[[working_project()]]$data_wide)
				req(DataInSets[[working_project()]]$group_order)
				req(DataInSets[[working_project()]]$sample_order)

				MetaData = DataInSets[[working_project()]]$MetaData
				data_long = DataInSets[[working_project()]]$data_long
				data_wide = DataInSets[[working_project()]]$data_wide
				input_groups = DataInSets[[working_project()]]$group_order
				input_samples = DataInSets[[working_project()]]$sample_order
				genelabel <- input$genelabel

				input_keep = which(MetaData$sampleid %in% input_samples)
				data_wide  <- data_wide[apply(data_wide, 1, function(x) sum(length(which(x==0 | is.na(x)))) < 3),]
				tmp_data_wide = data_wide[,input_keep] %>% as.matrix()

				tmp_data_long = dplyr::filter(data_long, (group %in% input_groups) & (sampleid %in% input_samples))
				tmp_group = MetaData$group[input_keep]
				tmp_sampleid = MetaData$sampleid[input_keep]

				return(list('tmp_data_wide'=tmp_data_wide,'tmp_data_long'=tmp_data_long,'tmp_group' = tmp_group, 'tmp_sampleid'=tmp_sampleid, "MetaData"=MetaData[input_keep, ] ))
			})

			DataPCAReactive <- reactive({
				DataQC <-  DataQCReactive()
				tmp_sampleid <- DataQC$tmp_sampleid
				validate(need(length(tmp_sampleid)>1, message = "Please select at least two samples (please note samples are filtered by group selection as well)."))
				tmp_data_wide <- DataQC$tmp_data_wide
				tmp_group = DataQC$tmp_group

				MaxPCANum <- as.numeric(input$MaxPCANum)
				tmp_data_wide[is.na(tmp_data_wide)] <- 0
				pca <- 	stats::prcomp(t(tmp_data_wide), rank. = MaxPCANum, scale. = FALSE)
				percentVar <- 	round((pca$sdev)^2/sum(pca$sdev^2), 3) * 100
				scores <- as.data.frame(pca$x)
				loading <- as.data.frame(pca$rotation)

				rownames(scores) <- tmp_sampleid
				scores$group <- factor(tmp_group, levels = DataInSets[[working_project()]]$group_order)
				attributes=setdiff(colnames(DataQC$MetaData), c("Order", "ComparePairs", "group") )
				MetaData=DataQC$MetaData
				colsel=match(attributes, colnames(MetaData) )
				scores=cbind(scores, MetaData[, colsel, drop=F])

				return(list('scores'=scores,'loading'=loading, 'percentVar'=percentVar, 'PCA'=pca))
			})

			#Eigenvalue bar chart
			Eigenvalues_plot <- reactive({
				req(DataPCAReactive())
				PCAlist <- DataPCAReactive()
				pca <- PCAlist$PCA
				percentVar <- PCAlist$percentVar
				maxPercent <- as.numeric(round(percentVar[[1]]))
				MaxPCANum <- as.numeric(input$MaxPCANum)
				p <- fviz_eig(pca, ncp = MaxPCANum, addlabels=TRUE,  linecolor ="red") +
				#ylim(0, maxPercent + 5) +
				labs(x = "Principal Components", y = "Percentage of variances") +
				theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18))
				#print(get_eig(pca))
				return(p)
			})

			output$Eigenvalues <- renderPlot({
				Eigenvalues_plot()
			})

			######## PCA
			plot_pca_control <- reactiveVal(0)

			observeEvent(input$plot_PCA, {
				plot_pca_control(plot_pca_control()+1)
			})

			#pcaplot_out <- eventReactive (plot_pca_control(), {
			pcaplot_out <- reactive({
				req(DataPCAReactive())
				pcnum <- as.numeric(input$pcnum)
				validate(need(length(pcnum)==2, message = "Select 2 Prinical Components."))

				PCAlist <- DataPCAReactive()
				scores <- PCAlist$scores
				percentVar <- PCAlist$percentVar
				samples=scores$sampleid

				xlabel <- paste("PC",pcnum[1],"(",round(percentVar[pcnum[1]]),"%)",sep="")
				ylabel <- paste("PC",pcnum[2],"(",round(percentVar[pcnum[2]]),"%)",sep="")

				PC1 <- paste("PC",pcnum[1],sep="")
				PC2 <- paste("PC",pcnum[2],sep="")

				#n <- length(unique(as.character(unlist(scores[, colnames(scores)==input$PCAcolorby]))))
				#colorpal = colorRampPalette(brewer.pal(8, input$colpalette))(n)

				colorpal <- UserColorPlalette(colpalette = input$colpalette, items = as.character(unlist(scores[, colnames(scores)==input$PCAcolorby])))

				if (input$PCA_subsample=="None" ) {
					labels=NULL
				} else {
					label_sel=match(input$PCA_label, names(scores))
					labels=unlist(scores[, label_sel])
					if (input$PCA_subsample=="Subset") {
						PCA_list=str_split(input$PCA_list, "\n")[[1]]
						N_sel=match(PCA_list, samples)
						N_sel=N_sel[!is.na(N_sel)]
						validate(need(length(N_sel)>0, message = "Enter at least one valid sampleid to label"))
						keep_s=rep(FALSE, length(labels))
						keep_s[N_sel]=TRUE
						labels[!keep_s]=""
					}
				}

				if (input$PCAshapeby=="none") {
					shape_by=19
				} else {
					shape_by=input$PCAshapeby
				}

				#if (input$PCAsizeby=="none") {
				size_by=as.numeric(input$PCAdotsize)
				#} else {
				#	size_by= input$PCAsizeby
				#}

				colorby = input$PCAcolorby

				if (is.numeric(scores[[input$PCAcolorby]])) {  #when colorby is numeric, don't use color palette
					p <- ggpubr::ggscatter(scores, x =PC1, y=PC2, color = colorby, shape=shape_by, size = size_by, ellipse = input$ellipsoid, mean.point = input$mean_point, rug = input$rug, label =labels, font.label = input$PCAfontsize, repel = TRUE,  ggtheme = theme_bw(base_size = 20) )
				} else {
					p <- ggpubr::ggscatter(scores, x =PC1, y=PC2, color = colorby, shape=shape_by, size = size_by, palette= colorpal, ellipse = input$ellipsoid, mean.point = input$mean_point, rug = input$rug,	label =labels, font.label = input$PCAfontsize, repel = TRUE,  ggtheme = theme_bw(base_size = 20) )
				}

				p <- ggpubr::ggpar(p, xlab = xlabel, ylab = ylabel)

				if (input$convex == "Yes") {
					p <- p + ggpubr::stat_chull(aes_string(color = colorby, fill = colorby), alpha = 0.1, geom = "polygon")
				}

				if (input$loading == "Yes") {
					ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
					num_factors <- input$contrib_factors
					genelabel <- input$genelabel
					x <- paste("PC",pcnum[1], sep="")
					y <- paste("PC",pcnum[2], sep="")

					UniqueID_selected <- data.frame(UniqueID=rownames(PCAlist$loading), abs(PCAlist$loading)) %>%
					dplyr::mutate(Max = pmax(get(x),get(y))) %>%
					dplyr::arrange(desc(Max)) %>%
					dplyr:: slice_head(n = num_factors) %>%
					dplyr::pull(UniqueID)

					datapc <- data.frame(UniqueID=rownames(PCAlist$loading), PCAlist$loading) %>%
					dplyr::filter(UniqueID %in% UniqueID_selected) %>%
					dplyr::inner_join(ProteinGeneName, by = "UniqueID")
					datapc$labelgeneid = datapc[,match(genelabel,colnames(datapc))]

					mult <- min(
						(max(scores[,y]) - min(scores[,y])/(max(datapc[,y])-min(datapc[,y]))),
						(max(scores[,x]) - min(scores[,x])/(max(datapc[,x])-min(datapc[,x])))
					)

					datapc <- transform(datapc,
						v1 = .7 * mult * (get(x)),
						v2 = .7 * mult * (get(y))
					)

					p <- p + geom_text_repel(data=datapc, aes(x=v1, y=v2, label=labelgeneid), size = 5, vjust=1, color="red") +
					geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="black")
				}
				#	p <- ggpubr::ggpar(p, legend.title ="", xlab = xlabel, ylab = ylabel, legend = "bottom") #works only when use color by.
				p <- p + guides(color = guide_legend(override.aes = list(label="")))

				return(p)
			})

			output$pcaplot <- renderPlot({
				withProgress(message = 'Drawing PCA Plot...', value = 0, {
					print(pcaplot_out())
				})
			})

			observeEvent(input$pcaplot, {
				saved_plots$pcaplot <- pcaplot_out()
			})

			output$pca_legend <- renderPlot({
				PCAlist <- DataPCAReactive()
				scores <- PCAlist$scores
				PCAcolorby <- input$PCAcolorby
				#colpalette <- input$colpalette

				tmp_group = as.character(unlist(scores[, colnames(scores)==PCAcolorby]))
				PCAcolorby = sym(PCAcolorby)

				#n <- length(unique(tmp_group))
				#colorpal = colorRampPalette(brewer.pal(8, colpalette))(n)

				colorpal <- UserColorPlalette(colpalette = input$colpalette, items = tmp_group)

				tmp_plot <- ggplot(scores, aes(x=PC1, y=PC2, color = !!PCAcolorby)) +
				geom_point() +
				scale_color_manual(values=colorpal) +
				theme_cowplot(12)

				legend_only <- get_legend(tmp_plot + theme(legend.position = "bottom",  legend.title = element_text(size = 16),
				legend.text = element_text(size = 14)) +
					guides(color = guide_legend(override.aes = list(size=8)))
				)

				cowplot::plot_grid(legend_only)
			})

			######## PCA 3D
			output$plot3d <- rgl::renderRglwidget({
				PCAlist <- DataPCAReactive()
				scores <- PCAlist$scores
				percentVar <- PCAlist$percentVar
				PCAcolorby <- input$PCAcolorby
				#colpalette <- input$colpalette

				xlabel <- paste("PC1(",round(percentVar[1]),"%)",sep="")
				ylabel <- paste("PC2(",round(percentVar[2]),"%)",sep="")
				zlabel <- paste("PC3(",round(percentVar[3]),"%)",sep="")

				sampleid <- rownames(scores)

				tmp_group = as.character(unlist(scores[, colnames(scores)==PCAcolorby]))
				#n <- length(unique(tmp_group))
				#colorpal = colorRampPalette(brewer.pal(8, colpalette))(n)

				colorpal <- UserColorPlalette(colpalette = input$colpalette, items = tmp_group)

				scores$tmp_group <- scores[, colnames(scores)==PCAcolorby]

				#rgl.open(useNULL=T)
				options(rgl.useNULL=TRUE)
				if (input$ellipsoid3d == "Yes") {
					ellipsoid3d = TRUE
				} else {
					ellipsoid3d = FALSE
				}

				if (any(table(tmp_group) <= 3))
				ellipsoid3d = FALSE

				if (input$dotlabel == "Yes") {
					dotlabel=TRUE
				} else {
					dotlabel=FALSE
				}

				car::scatter3d(PC3 ~ PC1 + PC2 | tmp_group, data= scores,
					axis.col= c("black", "black", "black"),
					xlab=xlabel, ylab=ylabel,  zlab=zlabel, labels = as.factor(sampleid), id=dotlabel, id.n=length(sampleid),
					axis.scales=FALSE,
					axis.ticks=FALSE,
					ellipsoid = ellipsoid3d,
					surface=FALSE,
					grid = FALSE,
					cex.lab=3,
					surface.alpha=1,
					surface.col = colorpal,
					point.col = colorpal
				)
				#legend3d("right", legend = levels(scores$tmp_group),  col = colorpal, pch = 12)
				#text3d(x=1.1, y=c(.9,1,1.1), z=1.1, levels(scores$tmp_group), col = colorpal)
				#points3d(x=1.2, y=c(.9,1,1.1), z=1.1, col = colorpal, size=5)
				rglwidget(width = 800, height = 800)
			})

			output$plotly3d <- renderPlotly({
				PCAlist <- DataPCAReactive()
				scores <- PCAlist$scores
				scores <- scores%>%mutate_if(is_character, as.factor)
				percentVar <- PCAlist$percentVar
				symbol_list=rep(c('circle', 'square',  'diamond',  'circle-open','square-open','diamond-open'), 2) #symbols which work with plotly scatter3d
				plot_symbols=symbol_list[unique(as.numeric(unlist(scores[, colnames(scores)==input$PCAshapeby])))]

				xlabel <- paste("PC1(",round(percentVar[1]),"%)",sep="")
				ylabel <- paste("PC2(",round(percentVar[2]),"%)",sep="")
				zlabel <- paste("PC3(",round(percentVar[3]),"%)",sep="")

				sampleid <- str_c(scores$sampleid, "\n", scores$group)
				#n <- length(unique(as.character(unlist(scores[, colnames(scores)==input$PCAcolorby]))))
				#colorpal = colorRampPalette(brewer.pal(8, input$colpalette))(n)

				colorpal <- UserColorPlalette(colpalette = input$colpalette, items = unique(as.character(unlist(scores[, colnames(scores)==input$PCAcolorby]))))

				if (input$PCAshapeby=="none"){
					p <- plot_ly(scores, x = ~PC1, y = ~PC2, z = ~PC3, color = as.formula(paste0("~", input$PCAcolorby)),
					colors = colorpal,text = sampleid) %>%
					add_markers() %>%
					layout(scene = list(xaxis = list(title = xlabel), yaxis = list(title = ylabel),  zaxis = list(title = zlabel)))

				} else{
					p <- plot_ly(scores, x = ~PC1, y = ~PC2, z = ~PC3, color = as.formula(paste0("~", input$PCAcolorby)),
						symbol=as.formula(paste0("~", input$PCAshapeby)),symbols=plot_symbols,
					colors = colorpal,text = sampleid) %>%
					add_markers() %>%
					layout(scene = list(xaxis = list(title = xlabel), yaxis = list(title = ylabel),  zaxis = list(title = zlabel)))
				}
				p$elementId <- NULL
				p
			})

			data_return <- reactive({
				PCAlist <- DataPCAReactive()
				loading <- as.data.frame(PCAlist$loading)
				loading <- loading %>% mutate_if(is.numeric, round, 6)

				if (input$abs == "Yes") {
					loading <- abs(loading)
				}

				return(loading)
			})

			output$loadingPCATable <- DT::renderDataTable({
				ttab <- DT::datatable(
					data_return(),
					rownames = TRUE,
					escape = FALSE,
					options = list(pageLength = 15),
					extensions = 'Buttons'
				)
			})

			# PCA Loading Download handler, conducts download function
			output$loadingPCA_download_button <- downloadHandler(
				filename = "PCA_Loading.csv",
				content = function(file) {S
					write.csv(data_return(), file)
				}
			)

			########## boxplot
			QCboxplot_out <- reactive({
				#QCboxplot_out <- eventReactive(input$tabset=="Box Plot", {

				#withProgress(message = 'Making box plot', value = 0, {
				DataQC <-  DataQCReactive()
				tmp_sampleid <- DataQC$tmp_sampleid
				tmp_data_long <- DataQC$tmp_data_long %>% dplyr::filter(expr !=0) %>% sample_n(1000)
				tmp_data_long$sampleid <- factor(tmp_data_long$sampleid, levels=DataInSets[[working_project()]]$sample_order)

				tmp_group = DataQC$tmp_group
				#colpalette <- input$colpalette

				#n <- brewer.pal.info[input$colpalette,'maxcolors']
				#if (length(unique(tmp_data_long$sampleid)) <=  n) {
				#	user_color <- RColorBrewer::brewer.pal(length(unique(tmp_data_long$sampleid)), input$colpalette)
				#} else {
				#	user_color=colorRampPalette(RColorBrewer::brewer.pal(n, input$colpalette))(length(unique(tmp_data_long$sampleid)))
				#}

				colorpal <- UserColorPlalette(colpalette = input$colpalette, items = unique(tmp_data_long$sampleid))

				p <- ggplot(tmp_data_long, aes(x=sampleid, y=expr)) +
				geom_boxplot(aes(color = factor(sampleid)), outlier.colour = NA) +
				scale_fill_manual(values = colorpal) +
				scale_colour_manual(values = colorpal) +
				coord_cartesian(ylim = range(boxplot(tmp_data_long$expr, plot=FALSE)$stats)*c(.9, 1.2)) +
				ylab(input$Ylab) +
				xlab(input$Xlab) +
				theme_bw(base_size = 20) +
				theme(legend.position = "none",
					legend.title=element_blank(),
					plot.margin=unit(c(1,1,1,1),"mm"),
					text = element_text(size=input$axisfontsize),
					axis.text.x = element_text(angle = input$Xangle, hjust=0.5, vjust=0.5)
				) +
				guides(col = guide_legend(ncol = 8))
				return(p)
				#})
			})

			output$QCboxplot <- renderPlot({
				QCboxplot_out()
			})

			observeEvent(input$QCboxplot, {
				saved_plots$QCboxplot <- QCboxplot_out()
			})

			############heatmap
			pheatmap_out <- reactive({
				#pheatmap_out <- eventReactive(input$tabset=="Sample-sample Distance", {
				DataQC <-  DataQCReactive()
				tmp_sampleid <- DataQC$tmp_sampleid
				tmp_data_wide <- DataQC$tmp_data_wide
				tmp_group = DataQC$tmp_group
				MetaData = DataQC$MetaData
				selCol = input$PCAcolorby

				annotation <- MetaData %>% dplyr::select(c(one_of(c("sampleid","group", selCol)))) %>%
				tibble::remove_rownames() %>%
				tibble::column_to_rownames(var = "sampleid") %>%
				dplyr::mutate(group = factor(group, levels=DataInSets[[working_project()]]$group_order))
				sample_annot = HeatmapAnnotation(df = annotation)

				#row_annot = rowAnnotation(df = annotation)
				sampleDistMatrix <- as.matrix(dist(t(tmp_data_wide)))
				rownames(sampleDistMatrix) <- tmp_sampleid
				colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(32)
				p <-  ComplexHeatmap::Heatmap(sampleDistMatrix, top_annotation = sample_annot, col=colors)
				return(p)
			})

			output$pheatmap <- renderPlot({
				ht = pheatmap_out()
				draw(ht, merge_legend = T,  auto_adjust = FALSE)
			})

			observeEvent(input$SampleDistance, {
				saved_plots$SampleDistance <- pheatmap_out()
			})

			############Dendrograms
			Dendrograms_out <- reactive({
				#Dendrograms_out <- eventReactive(input$tabset=="Dendrograms", {
				DataQC <-  DataQCReactive()
				tmp_data_wide <- DataQC$tmp_data_wide
				hc <- hclust(dist(t(tmp_data_wide)))

				if (input$dendroformat=="tree") {
					p <- factoextra::fviz_dend(hc, k = input$DendroCut, cex = input$DendroFont, k_colors =input$colpalette,	color_labels_by_k = TRUE, rect = TRUE, rect_border = "jco",	rect_fill = TRUE)
				} else 	if (input$dendroformat=="horiz") {
					p <- factoextra::fviz_dend(hc, k = input$DendroCut, cex = input$DendroFont, k_colors = input$colpalette,	 horiz = TRUE, color_labels_by_k = TRUE, rect = TRUE, rect_border = "jco", rect_fill = TRUE)
				} else if (input$dendroformat=="circular") {
					p <- factoextra::fviz_dend(hc, k = input$DendroCut, cex = input$DendroFont, k_colors = input$colpalette, type = "circular")
				}
				return(p)
			})

			output$Dendrograms <- renderPlot({
				Dendrograms_out()
			})

			observeEvent(input$Dendrograms, {
				saved_plots$Dendrograms <- Dendrograms_out()
			})

			############histplot
			histplot_out <- reactive({
				#histplot_out <- eventReactive(input$tabset=="CV Distribution", {
				withProgress(message = 'Calculating.',  detail = 'This may take a while...', value = 0, {
					DataQC <-  DataQCReactive()
					tmp_sampleid <- DataQC$tmp_sampleid
					tmp_data_long <- DataQC$tmp_data_long
					tmp_group = DataQC$tmp_group
					CV.df <- tmp_data_long %>%
					dplyr::group_by(.,  group, id) %>%
					dplyr::summarise( mean=mean(expr, na.rm = TRUE), sd=sd(expr, na.rm = TRUE), .groups = 'drop') %>%
					dplyr::mutate(CV=100*(sd/mean)) %>%
					dplyr::ungroup()

					mu <- group_by(CV.df, group) %>%
					dplyr::summarise(median = round(median(CV, na.rm = TRUE),1), .groups = 'drop') %>%
					dplyr::ungroup()

					interval <- seq.int(0, 100, 5)
					xlimmin <- interval[cut(min(mu$median), interval, include.lowest = TRUE, labels = FALSE)]
					xlimmax <- interval[cut(max(mu$median), interval, include.lowest = TRUE, labels = FALSE) +1]

					#n <- brewer.pal.info[input$colpalette,'maxcolors']
					#if (length(unique(tmp_group)) <=  n) {
					#	user_color <- RColorBrewer::brewer.pal(length(unique(tmp_group)), input$colpalette)
					#} else {
					#	user_color=colorRampPalette(RColorBrewer::brewer.pal(n, input$colpalette))(length(unique(tmp_group)))
					#}
					#
					colorpal <- UserColorPlalette(colpalette = input$colpalette, items = unique(tmp_group))


					p <- ggplot(CV.df, aes(x=CV, color = group)) +
					scale_fill_manual(values = colorpal) +
					scale_colour_manual(values = colorpal) +
					geom_freqpoly (position="dodge", na.rm = TRUE, bins = 10) +
					geom_vline(data=mu, aes(xintercept=median, color=group), linetype="dashed") +
					geom_text(data=mu, mapping=aes(x=median, y=0, label=paste(median,"(",group,")", sep="")), size=4, angle=90, vjust=-0.4, hjust=0) +
					scale_x_continuous(breaks = seq(xlimmin, xlimmax, by=5), limits=c(xlimmin,xlimmax)) +
					theme_bw(base_size = 20) +
					theme(legend.position = "bottom")
					return(p)
				})
			})

			output$histplot <- renderPlot({
				histplot_out()
			})

			observeEvent(input$histplot, {
				saved_plots$histplot <- histplot_out()
			})

			############PC_covariates QC Plots
			PC_covariates_out <-  eventReactive(input$compute_PC,{
				DataQC <-  DataQCReactive()
				tmp_data_wide <- DataQC$tmp_data_wide
				MetaData=DataQC$MetaData
				meta=MetaData[, !(colnames(MetaData) %in% c("sampleid", "Order", "ComparePairs")), drop=FALSE]
				meta=meta[, (colnames(meta) %in% input$covar_variates), drop=FALSE]
				rownames(meta)=MetaData$sampleid
				res<-Covariate_PC_Analysis(tmp_data_wide, meta, out_prefix=NULL, PC_cutoff=input$covar_PC_cutoff,
				FDR_cutoff=input$covar_FDR_cutoff, N_col=input$covar_ncol)
				return(res)
			})

			output$covar_table <- DT::renderDT(server=FALSE,{
				results <- PC_covariates_out()$selVar_All
				if (!is.null(results)) {
					results["P-value"]=as.numeric(formatC(unlist(results["P-value"]), format="e", digits=2))
					results["FDR"]=as.numeric(formatC(unlist(results["FDR"]), format="e", digits=2))
				}

				DT::datatable(results,  extensions = 'Buttons',
					options = list(
						dom = 'lBfrtip', pageLength = 25,
						buttons = list(
							list(extend = "csv", text = "Download Page", filename = "Page_results",	exportOptions = list(modifier = list(page = "current"))),
							list(extend = "csv", text = "Download All", filename = "All_Results",	exportOptions = list(modifier = list(page = "all")))
						)
					), rownames= T
				)
			})

			output$PC_covariatesC <- renderPlot({
				data=PC_covariates_out()$sel_dataC
				if (!is.null(data)) {
					data$plot
				}
			})

			output$plot.PC_covariatesC=renderUI({
				tagList(
					textOutput(ns("N_pairs_C")),
					plotOutput(ns("PC_covariatesC"), height = input$covar_cat_height)
				)
			})

			output$PC_covariatesN <- renderPlot({
				data=PC_covariates_out()$sel_dataN
				if (!is.null(data)) {
					data$plot
				}
			})

			output$plot.PC_covariatesN=renderUI({
				tagList(
					textOutput(ns("N_pairs_N")),
					plotOutput(ns("PC_covariatesN"), height = input$covar_num_height)
				)
			})

			Npairs_cov<-reactive({
				res<-PC_covariates_out()
				C=res$sel_dataC$selVar
				if (is.null(C)) {N1=0} else {N1=nrow(C)}
				N=res$sel_dataN$selVar
				if (is.null(N)) {N2=0} else {N2=nrow(N)}
				return(c(N1, N2))
			})

			observe({
				H_C=ceiling(Npairs_cov()[1]/PC_covariates_out()$ncol)*400
				if (H_C>0)  { updateSliderInput(session, "covar_cat_height", value = H_C)}
				H_N=ceiling(Npairs_cov()[2]/PC_covariates_out()$ncol)*400
				if (H_N>0)  { updateSliderInput(session, "covar_num_height", value = H_N)}
			})

			output$N_pairs_C<-renderText({str_c("There are ", Npairs_cov()[1], " significant categorical covariate-PC pairs.")})
			output$N_pairs_N<-renderText({str_c("There are ", Npairs_cov()[2], " significant numeric covariate-PC pairs.")})
			output$N_pairs<-renderText({str_c("There are ", Npairs_cov()[1]+Npairs_cov()[2], " significant covariate-PC pairs.")})


			observeEvent(input$covar_cat, {
				data=PC_covariates_out()$sel_dataC
				saved_plots$covar_cat <- data$plot
			})

			observeEvent(input$covar_num, {
				data=PC_covariates_out()$sel_dataN
				saved_plots$covar_num<- data$plot
			})

		}
	)
}
