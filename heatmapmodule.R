###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: heatmap.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 05/23/2023
##@version 3.0
###########################################################################################################
#pkgs:"heatmaply","dendextend","ComplexHeatmap","circlize", "gplots", "dplyr", "stringr"

library(heatmaply)
library(dendextend)
library(ComplexHeatmap)
library(circlize)
library(gplots) #heatmap.2

heatmap_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
		column(3,
			wellPanel(
				uiOutput(ns('loadedprojects')),
				uiOutput(ns("selectGroupSample")),
				radioButtons(ns("subset"), label="Genes used for heatmap", choices=c("All","Subset","Upload Genes", "Geneset"), inline = TRUE, selected="All"),
				conditionalPanel(ns = ns, "input.subset=='Upload Genes'",
					radioButtons(ns("upload_type"), label="Select upload type", inline = TRUE, choices = c("Gene List","Annotated Gene File"), selected = "Gene List"),
					conditionalPanel(ns = ns, "input.upload_type=='Gene List'",
						textAreaInput(ns("list"), "Enter Gene List", "", cols = 5, rows=6)
					),
					conditionalPanel(ns = ns, "input.upload_type=='Annotated Gene File'",
						uiOutput(ns("gene_annot_file"))
					)
				),
				conditionalPanel(ns = ns, "input.subset=='Geneset'",
					selectizeInput(ns("sel_geneset"), label="Available GeneSet", choices = NULL, multiple = FALSE),
					textAreaInput(ns("geneset_genes"), "Genes in Geneset", "", cols = 5, rows=6)
				),
				conditionalPanel(ns = ns, "input.subset=='All'",
					radioButtons(ns("submethod"), label= "Plot Random Genes or Variable Genes", choices= c("Random"="Random","Variable"="Variable"), inline = TRUE),
					numericInput(ns("maxgenes"), label="Choose Gene Number", min=1, max= 5000, value=100, step=1)
				),
				conditionalPanel(ns = ns ,"input.subset=='Subset'",
					selectInput(ns("test"), label="Select Genes from Test:", choices=NULL),
					fluidRow(
						column(width=6, numericInput(ns("fccut"), label= "Fold Change Cutoff", value = 1.2, min=1, step=0.1)),
						column(width=6, numericInput(ns("pvalcut"), label= "P Value Cutoff", value=0.01, min=0, step=0.001))
					),
					radioButtons(ns("psel"), label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"),inline = TRUE),
					span(textOutput(ns("filteredgene")), style = "color:red; font-size:15px; font-family:arial; font-style:italic"),
					uiOutput(ns("Test_to_sample")),
					tags$hr(style="border-color: black;")
				),
				conditionalPanel(ns = ns, "input.tabset=='Static Heatmap Layout 1'",
					selectizeInput(ns("annot"), label="Annotate Samples", choices=NULL, multiple = TRUE)
				),
				fluidRow(
					column(width=6, selectInput(ns("dendrogram"), "Apply Clustering:", c("both" ,"none", "row", "column"), selected="row")),
					column(width=6, selectInput(ns("scale"), "Apply Scaling:", c("none","row", "column"), selected="row"))
				),
				conditionalPanel(ns = ns, "input.tabset=='Static Heatmap Layout 2'",
					fluidRow(
						column(width=6, selectInput(ns("key"), "Color Key:", c("TRUE", "FALSE"))),
						column(width=6, selectInput(ns("srtCol"), "angle of label", c("45", "60","90")))
					),
					fluidRow(
						column(width=6, sliderInput(ns("hxfontsize"), "Column Font Size:", min = 0, max = 3, step = 0.5, value = 1)),
						column(width=6, sliderInput(ns("hyfontsize"), "Row Font Size:", min = 0, max = 3, step = 0.5, value = 1))
					),
					fluidRow(
						column(width=6, sliderInput(ns("right"), "Set Margin Width", min = 0, max = 20, value = 5)),
						column(width=6, sliderInput(ns("bottom"), "Set Margin Height", min = 0, max = 20, value = 5))
					)
				),
				conditionalPanel(ns = ns, "input.tabset=='Interactive Heatmap'",
					fluidRow(
						column(width=6, selectInput(ns("key"), "Color Key:", c("TRUE", "FALSE"))),
						column(width=6, selectInput(ns("srtCol"), "angle of label", c("45", "60","90")))
					),
					fluidRow(
						column(width=6, sliderInput(ns("hxfontsizei"), "Column Font Size:", min = 0, max = 3, step = 0.5, value = 1)),
						column(width=6, sliderInput(ns("hyfontsizei"), "Row Font Size:", min = 0, max = 3, step = 0.5, value = 1))
					),
					fluidRow(
						column(width=6, sliderInput(ns("l"), "Set Margin Width", min = 0, max = 200, value = 120)),
						column(width=6,	sliderInput(ns("b"), "Set Margin Height", min = 0, max = 200, value = 120))
					)
				),
				conditionalPanel(ns = ns, "input.tabset=='Static Heatmap Layout 1'",
					fluidRow(
						column(width=6, sliderInput(ns("hxfontsizep"), "Column Font Size:", min = 0, max = 20, step = 1, value = 10)),
						column(width=6, sliderInput(ns("hyfontsizep"), "Row Font Size:", min = 0, max = 20, step = 1, value = 7))
					),
					radioButtons(ns("label"),label="Gene Label",inline = TRUE, choices=""),
					sliderInput(ns("N_genes"), "Max Number of Genes to Label:", min = 0, max = 500, step = 10, value = 100),
					radioButtons(ns("highlight"), label="Highlight Subset of Genes:", inline = TRUE, choices = c("Yes","No"), selected = "No"),
					conditionalPanel(ns = ns, "input.highlight=='Yes'",
						uiOutput(ns("gene_highlight_file")),
						sliderInput(ns("hl_font_size"), "Font Size:", min = 0, max = 20, step = 1, value = 9)
					),
					sliderInput(ns("height"), "Heatmap Height:", min = 200, max = 3000, step = 50, value = 800),
					fluidRow(
						column(width=6, radioButtons(ns("row_dend"), label="Show: Row Dendrogram", inline = TRUE, choices = c("No" = FALSE,"Yes" = TRUE), selected = TRUE)),
						column(width=6, radioButtons(ns("col_dend"), label="Column Dendrogram", inline = TRUE, choices = c("No" = FALSE,"Yes" = TRUE), selected = TRUE))
					),
					fluidRow(
						column(width=4, colourpicker::colourInput(ns("lowColor"), "Low", "blue")),
						column(width=4, colourpicker::colourInput(ns("midColor"), "Mid", "white")),
						column(width=4, colourpicker::colourInput(ns("highColor"), "High", "red"))
					),
					fluidRow(
						column(width=6, selectInput(ns("distanceMethod"), "Distance Metric:", c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"))),
						column(width=6, selectInput(ns("agglomerationMethod"), "Linkage Algorithm:", c("complete", "single", "average", "centroid", "median", "mcquitty", "ward.D", "ward.D2")))
					),
					fluidRow(
						column(width=6, sliderInput(ns("cutreerows"), "cutree_rows:", min = 0, max = 8, step = 1, value = 0)),
						column(width=6, sliderInput(ns("cutreecols"), "cutree_cols:", min = 0, max = 8, step = 1, value = 0))
					),
					radioButtons(ns("custom_color"), label="Upload Colors for Annotations", inline = TRUE, choices = c("Yes","No"), selected = "No"),
					conditionalPanel(ns = ns, "input.custom_color=='Yes'",
						uiOutput(ns("annot_color_file"))
					),
					fluidRow(
						column(width=6, textInput(ns("row_title"), "Row Title", width = "100%")),
						column(width=6, sliderInput(ns("row_title_font_size"), "Row Title Font Size:", min = 0, max = 30, step = 1, value = 16))
					),
					fluidRow(
						column(width=6, textInput(ns("column_title"), "Column Title", width = "100%")),
						column(width=6, sliderInput(ns("column_title_font_size"), "Column Title Font Size:", min = 0, max = 30, step = 1, value = 16))
					),
					h5("After changing parameters, please click Plot/Refresh button in the plot panel to generate heatmap.")
				)
			)
		),
		column(9,
			tabsetPanel(id=ns("tabset"),
				tabPanel(title="Static Heatmap Layout 1", value="Static Heatmap Layout 1",
					actionButton(ns("plot_heatmap"), "Plot/Refresh"),
					actionButton(ns("pheatmap2"), "Save to output"),
					uiOutput(ns("plot.heatmap"))
				),
				tabPanel(title="Static Heatmap Layout 2", value="Static Heatmap Layout 2",
					actionButton(ns("staticheatmap"), "Save to output"),
					plotOutput(ns("staticheatmap"), height = 800)
				),
				tabPanel(title="Interactive Heatmap", value="Interactive Heatmap",
					actionButton(ns("action_heatmaps"), "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
					p(),
					plotlyOutput(ns("interactiveheatmap"), height = 800)
				),
				tabPanel(title="Help", value="Help",
					htmlOutput("help_heatmap")
				)
			)
		)
	)
}

heatmap_server <- function(id) {
	shiny::moduleServer(id,
		function(input, output, session){
			ns <- session$ns
			output$loadedprojects <- renderUI({
				req(length(working_project()) > 0)
				radioButtons(ns("current_dataset"), label = "Change Working Dataset", choices=names(DataInSets), inline = F, selected=working_project())
			})

			observeEvent(input$current_dataset, {
				working_project(input$current_dataset)
			})

			output$selectGroupSample <-	renderUI({
				req(length(working_project()) > 0)
				sample_info <- paste("Selected ", length(DataInSets[[working_project()]]$group_order), " out of ", length(DataInSets[[working_project()]]$groups), " Groups, ",
					length(DataInSets[[working_project()]]$sample_order), " out of ", length(DataInSets[[working_project()]]$samples),
				" Samples. (Update Selection at: Top Menu -> Groups and Samples.)", sep="")
				tagList(
					tags$p(sample_info),
					tags$hr(style="border-color: black;")
				)
			})

			#observeEvent(DataInSets[[working_project()]]$MetaData, {
			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$MetaData)#cat("load file UI for", ProjectInfo$ProjectID, "\n")
				updateRadioButtons(session, "subset",  selected="All")
				output$gene_highlight_file=renderUI({
					tagList(fileInput("file_gene_highlight", "Highlight Genes (csv with headers like Genes, Pathways, Color)"))
				})
				updateRadioButtons(session, "highlight",  selected="No")
				output$gene_annot_file=renderUI({
					tagList(fileInput("file_gene_annot", "Choose gene annotation file (csv with headers like Genes, Pathways, Color)"))
				})
				updateRadioButtons(session, "custom_color",  selected="No")
				output$annot_color_file=renderUI({
					tagList(fileInput("annot_color_file", "Upload annotation Colors (csv with 3 headers: Attribute, Value and Color)"))
				})
				updateTabsetPanel(session, "Tables", selected = "Project Overview")
			})

			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$MetaData)
				MetaData = DataInSets[[working_project()]]$MetaData
				req(DataInSets[[working_project()]]$tests_order)
				tests = DataInSets[[working_project()]]$tests_order
				ProteinGeneNameHeader = DataInSets[[working_project()]]$ProteinGeneNameHeader
				updateSelectizeInput(session,'test',choices=tests, selected=tests[1])
				updateRadioButtons(session,'label', inline = TRUE, choices=ProteinGeneNameHeader, selected="Gene.Name")
				attributes = sort(setdiff(colnames(MetaData), c("sampleid", "Order", "ComparePairs") ))
				updateSelectInput(session, "annot", choices=attributes, selected="group")
			})

			output$plot.heatmap=renderUI({
				plotOutput(ns("pheatmap2"), height = input$height)
			})

			filteredGene = reactive({
				test = input$test
				fccut =log2(as.numeric(input$fccut))
				pvalcut = as.numeric(input$pvalcut)

				results_long = DataInSets[[working_project()]]$results_long

				if (input$psel == "Padj") {
					filteredGene = results_long %>% filter(test %in% test & abs(logFC) > fccut & Adj.P.Value < pvalcut) %>%
					dplyr::pull(UniqueID)
				} else {
					filteredGene = results_long %>% filter(test %in% test & abs(logFC) > fccut & P.Value < pvalcut) %>%
					dplyr::pull(UniqueID)
				}
				#cat("Selected Genes:",length(filteredGene), "\n") #debug
				return(filteredGene)
			})

			output$filteredgene <- renderText({ paste("Selected Genes:",length(filteredGene()),sep="")})

			#################
			observeEvent(input$subset , {
				req(length(working_project()) > 0)
				req(input$subset == "Geneset")
				genesetnames <- GetGeneSetNames()
				updateSelectizeInput(session, "sel_geneset", choices =  c('Type to Search' = '', genesetnames), server = TRUE)
			})

			observeEvent(input$sel_geneset, {
				req(length(working_project()) > 0)
				req(input$subset == "Geneset")
				req(input$sel_geneset!="")
				sel_geneset <- input$sel_geneset
				geneset_genenames <- GetGenesFromGeneSet(sel_geneset)
				updateTextAreaInput(session, "geneset_genes", value=paste(geneset_genenames, collapse=","))
			})
			###############

			DataHeatMapReactive <- reactive({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$results_long)
				results_long = DataInSets[[working_project()]]$results_long
				ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
				MetaData = DataInSets[[working_project()]]$MetaData
				tmpgroups = DataInSets[[working_project()]]$group_order
				tmpsamples = DataInSets[[working_project()]]$sample_order

				#tmpkeep = which((MetaData$group %in% tmpgroups)&(MetaData$sampleid %in% tmpsamples))
				tmpkeep = which(MetaData$sampleid %in% tmpsamples)
				gene_annot_info=NULL
				tmp_group = MetaData$group[tmpkeep]
				tmp_sampleid = MetaData$sampleid[tmpkeep]
				annotation = data.frame("group" = tmp_group, sampleid=tmp_sampleid)
				rownames(annotation) <- tmp_sampleid
				annotation <- annotation %>% left_join(MetaData, by = join_by(group, sampleid))
				annotation$group = factor(tmp_group, levels=DataInSets[[working_project()]]$group_order)

				if(length(tmpkeep)>0) {
					y <- DataInSets[[working_project()]]$group_order
					x= MetaData$group[tmpkeep]
					z = MetaData$sampleid[tmpkeep]
					new_order <- as.character(z[order(match(x, y))])
					tmpdat  <- DataInSets[[working_project()]]$data_wide %>% dplyr::select(all_of(new_order))
					tmpdat[is.na(tmpdat)] <- 0
					rownames(tmpdat) <-  rownames(DataInSets[[working_project()]]$data_wide)
				}

				if (input$subset == "Subset") {
					if(length(filteredGene())>0) {
						gene_list=intersect(rownames(tmpdat), filteredGene()) #user only genes that are in expression.
						tmpdat  <-  tmpdat[gene_list,]
					}
				}

				if (input$subset == "All") {
					if (nrow(tmpdat)>input$maxgenes) {
						if (input$submethod=="Random") {
							tmpdat=tmpdat[sample(1:nrow(tmpdat), input$maxgenes),] #this will keep rownames
							#tmpdat <- tmpdat %>% sample_n(input$maxgenes) #this will remove rownames
						} else {
							dataSD=apply(tmpdat, 1, function(x) sd(x,na.rm=T))
							dataM=rowMeans(tmpdat)
							diff=dataSD/(dataM+median(dataM)) #SD/mean, added median value to penalized lowly expressed genes
							tmpdat=tmpdat[order(diff, decreasing=TRUE)[1:input$maxgenes], ]
						}
					}
				}

				if (input$subset == "Upload Genes" | input$subset == "Geneset") {
					if (input$subset == "Upload Genes"){
						if (input$upload_type=='Gene List') {
							gene_list <- input$list
						} #else {
						#req(input$file_gene_annot)
						#annot_genes=read_csv(input$file_gene_annot$datapath)
						#gene_list <- unlist(annot_genes[, 1])
						#}
					}

					if (input$subset == "Geneset") {
						req(input$geneset_genes)
						gene_list  <- input$geneset_genes
					}

					gene_list <- ProcessUploadGeneList(gene_list)

					validate(need(length(gene_list )>2, message = "Please input at least 2 valid genes."))

					uploadlist <- dplyr::filter(ProteinGeneName, (UniqueID %in% gene_list ) | (Protein.ID %in% gene_list ) | (toupper(Gene.Name) %in% toupper(gene_list )))  %>%
					dplyr::pull(UniqueID)
					validate(need(length(uploadlist)>2, message = "Please input at least 2 valid genes."))

					#restore order of the input list
					sel1 = match(uploadlist, ProteinGeneName$UniqueID)
					ID_order <- ProteinGeneName[sel1, ] %>% mutate(N1=match(UniqueID, gene_list), N2=match(Protein.ID, gene_list), N3=match(toupper(Gene.Name), toupper(gene_list)), N=pmin(N1, N2, N3, na.rm=T))%>%arrange(N)
					tmpdat  <-  tmpdat[ID_order$UniqueID,]
					sel_rows1=rowSums(is.na(tmpdat))<ncol(tmpdat) #remove data rows with all NAs
					sel_rows2=rownames(tmpdat) %in% rownames(DataInSets[[working_project()]]$data_wide) #remove duplicate rows caused by matching (e.g."ALDH7A1_P49419"   "ALDH7A1_P49419-2")
					tmpdat  <-  tmpdat[sel_rows1 & sel_rows2, ]

					#if (input$upload_type=='Annotated Gene File') {
					#	gene_annot_info=data.frame(UniqueID=ID_order$UniqueID, annot_genes[ID_order$N, ])
					#	gene_annot_info=gene_annot_info[sel_rows1 & sel_rows2, ]
					#}
				}

				if (nrow(tmpdat)>5000 ) {tmpdat=tmpdat[sample(1:nrow(tmpdat), 5000),]; cat("Reduce data pionts to 5K\n")} #Use at most 5000 genes so the App won't crash

				df <- data.matrix(tmpdat)
				#use selected gene label
				sel = match(rownames(df), ProteinGeneName$UniqueID)
				selCol = match(input$label, names(ProteinGeneName))

				if (sum(is.na(sel))==0 & sum(is.na(selCol)==0)) {
					rownames(df) = unlist(ProteinGeneName[sel, selCol])
				} else {
					cat("gene lables not updated",sum(is.na(sel)), sum(is.na(selCol)), "\n")
				}
				#match sampleid order
				new_order = match(colnames(df), annotation$sampleid)
				annotation = annotation[new_order, ]
				return(list("df"=df, "annotation"=annotation, "gene_annot_info"=gene_annot_info))
			})

			pheatmap2_out <- eventReactive(input$plot_heatmap, {
				withProgress(message = 'Making static heatmap 1:', value = 0, {
					DataHeatMap <- DataHeatMapReactive()
					data.in <- DataHeatMap$df
					annotation <- DataHeatMap$annotation
					gene_annot_info <- DataHeatMap$gene_annot_info
					sample_annot=NULL #column annotation
					if (!is.null(input$annot)) {
						sel_col=match(input$annot, names(annotation))
						df_annot=annotation[, sel_col, drop=FALSE]
						sample_annot=HeatmapAnnotation(df = df_annot)
						if (input$custom_color=="Yes") {
							req(input$annot_color_file)
							annot_color=read_csv(input$annot_color_file$datapath)
							annot_color<-annot_color%>%dplyr::filter(Attribute %in% names(df_annot))
							#validate(need(nrow(annot_color)>0, message = "Please input valid annotate attributes."))
							if (nrow(annot_color)>0) {
								attr_list=unique(annot_color$Attribute)
								color_list=NULL
								#browser() #debug
								for (attr in attr_list) {
									subdata<-annot_color%>%filter(Attribute==attr)
									colorV=subdata$Color; names(colorV)=subdata$Value
									color_list[[attr]]=colorV
								}
								sample_annot=HeatmapAnnotation(df = df_annot, col=color_list)
							} else {
								cat("Annotation Color File Attributes not matching MetaData!\n")
							}
						}
					}
					cluster_rows = FALSE;
					cluster_cols=FALSE
					if (input$dendrogram == "both" | input$dendrogram == "row")
					cluster_rows = TRUE
					if (input$dendrogram == "both" | input$dendrogram == "column")
					cluster_cols = TRUE

					cexRow = as.numeric(as.character(input$hyfontsizep))
					cexCol = as.numeric(as.character(input$hxfontsizep))

					labCol = TRUE
					labRow = TRUE
					# cat("pheatmap ", dim(data.in), date(), "\n") #debug
					if (cexRow  == 0 | nrow(data.in) > input$N_genes) {
						labRow = FALSE
						cexRow = 5
					}
					if (cexCol == 0) {
						labCol = FALSE
						cexCol  = 5
					}

					cutree_rows = input$cutreerows
					cutree_cols = input$cutreecols

					#clean up SD=0 rows and columns
					if (input$scale=="row") {
						row_SD=apply(data.in, 1, function(x) sd(x,na.rm=T))
						data.in=data.in[row_SD!=0, ]
					}
					if (input$scale=="column") {
						col_SD=apply(data.in, 2, function(x) sd(x,na.rm=T))
						data.in=data.in[, col_SD!=0]
					}

					#now reproduce in Heatmap
					if (input$scale=="none") {
						data_range=quantile(unlist(data.in), probs=c(0.01, 0.5, 0.99), na.rm=T)
						col_fun=colorRamp2(data_range, c(input$lowColor,input$midColor, input$highColor) )
						legend_text="Value"
					} else {
						if (input$scale=="row") {
							data.in=t(scale(t(data.in)) )
						}	else {
							data.in=scale(data.in)
						}
						data_range=quantile(unlist(abs(data.in)), probs=c(0.01, 0.5, 0.99), na.rm=T)
						max_s=data_range[3]
						col_fun=circlize::colorRamp2(c(0-max_s, 0, max_s),  c(input$lowColor,input$midColor, input$highColor) )
						legend_text=str_c("Scaled Value")
					}
					if (cluster_cols==F) {cutree_cols=0}
					if (input$highlight=="No") {row_label_side="right"} else (row_label_side="left")


					p <- ComplexHeatmap::Heatmap(data.in, col=col_fun, cluster_rows = cluster_rows, cluster_columns = cluster_cols,
						clustering_distance_rows=input$distanceMethod, clustering_distance_columns=input$distanceMethod,
						clustering_method_rows=input$agglomerationMethod, clustering_method_columns=input$agglomerationMethod,
						row_km=cutree_rows, column_km=cutree_cols, row_km_repeats = 100, column_km_repeats = 100,
						show_row_names = labRow, show_column_names = labCol, row_names_side=row_label_side,
						show_row_dend=as.logical(input$row_dend), show_column_dend = as.logical(input$col_dend),
						top_annotation = sample_annot,	row_names_gp = gpar(fontsize = cexRow),
						column_names_gp = gpar(fontsize = cexCol), heatmap_legend_param = list(title = legend_text, color_bar = "continuous")
					)

					#highlight genes
					if (!is.null(gene_annot_info)) {
						df=gene_annot_info[, 3:ncol(gene_annot_info), drop=F]
						sel_col_path=match(c("Color", "Pathways"), names(df))
						if (sum(is.na(sel_col_path))==0) {
							df_color<-df%>%filter(!duplicated(Color))
							pathway_color=df_color$Color
							names(pathway_color)=df_color$Pathways
							rowAnnot=rowAnnotation(Pathways=gene_annot_info$Pathways, col=list(Pathways=pathway_color) )
						} else {
							rowAnnot=rowAnnotation(df=df)
						}
						p <- p+rowAnnot
					}


					if (input$highlight=="Yes"){
						req(input$file_gene_highlight)
						annot_genes=read_csv(input$file_gene_highlight$datapath)
						ccl <- which(toupper(rownames(data.in)) %in% toupper(annot_genes$gene_name) )
						validate(need(length(ccl)>0, message = "Please input at least one valid gene to highlight."))

						sel_col=match(toupper(rownames(data.in)[ccl]), toupper(annot_genes$gene_name) )
						ccl_color <- as.character(annot_genes$Color[sel_col])
						nameZoom = rowAnnotation(link = anno_mark(at = ccl, labels = rownames(data.in)[ccl], labels_gp = gpar(fontface = "bold",col = ccl_color,fontsize = input$hl_font_size), padding = 0.2))
						p <- p + nameZoom
						#Add pathway legend if no gene annotation
						hasP <- match("Pathways", names(annot_genes))

						if (is.null(gene_annot_info) & hasP) {
							Pathways=rep("", nrow(data.in));
							Pathways[ccl]=annot_genes$Pathways[sel_col]
							Pathways=str_wrap(Pathways,width=16)
							logjs(Pathways)
							legend_height = (max(str_count(Pathways,"\n"))+1) * 0.36
							Colors=rep("", nrow(data.in));
							Colors[ccl]=annot_genes$Color[sel_col]
							p_colors = structure(unique(as.character(Colors)), names=unique(as.character(Pathways)))
							p_colors = p_colors[-which(names(p_colors)=="")]

							pathway_legend <-  ComplexHeatmap::Heatmap( data.frame(Pathways), name = "Pathways",  rect_gp = gpar(type = "none"), show_column_names= FALSE, width = unit(0, "mm"), col = p_colors, heatmap_legend_param = list(title_position="topleft", labels_gp = gpar(lineheight=0.8), grid_height = unit(legend_height, "cm")))
							p <- p + pathway_legend
						}
					}
					return(p)
				})
			})

			output$pheatmap2 <- renderPlot({
				ht = pheatmap2_out()
				draw(ht, merge_legend = T, auto_adjust = FALSE)
			})

			observeEvent(input$pheatmap2, {
				saved_plots$pheatmap2 <- pheatmap2_out()
			})

			staticout <- reactive({
				withProgress(message = 'Making static heatmap 2:', value = 0, {
					DataHeatMap <- DataHeatMapReactive()

					data.in <- DataHeatMap$df
					annotation <- DataHeatMap$annotation

					cutree_rows = input$cutreerows
					cutree_cols = input$cutreecols
					if (cutree_rows == 0)
					cutree_rows = NULL
					if (cutree_cols == 0)
					cutree_cols = NULL

					if (input$dendrogram == "both" | input$dendrogram == "row")
					dend_r <- data.in %>% dist(method = input$distanceMethod) %>% hclust(method = input$agglomerationMethod) %>% as.dendrogram
					#%>% ladderize %>%  color_branches (k=cutree_rows)
					if (input$dendrogram == "both" | input$dendrogram == "column")
					dend_c <- t(data.in) %>% dist(method = input$distanceMethod) %>% hclust(method = input$agglomerationMethod) %>% as.dendrogram
					#%>% ladderize %>% color_branches (k=cutree_cols)

					#  cat(date(), dim(data.in), "layout 2\n") #debug
					cexRow = as.numeric(as.character(input$hyfontsize))
					cexCol = as.numeric(as.character(input$hxfontsize))

					labCol = labRow = NULL

					if (cexRow  == 0 | nrow(data.in) > 50) {
						labRow = FALSE
						cexRow = 0.2
					}

					if (cexCol == 0) {
						labCol = FALSE
						cexCol  = 0.2
					}

					p <-	heatmap.2(
						data.in,
						trace = "none",
						scale = input$scale,
						dendrogram = input$dendrogram,
						key = input$key,
						labRow = labRow,
						labCol = labCol,
						cexRow = cexRow,
						cexCol = cexCol,
						Rowv = if (input$dendrogram == "both" | input$dendrogram == "row") dend_r else FALSE,
						Colv = if (input$dendrogram == "both" | input$dendrogram == "column") dend_c else FALSE,
						col = colorpanel (32, low = input$lowColor,mid = input$midColor, high = input$highColor),
						srtCol = as.numeric(as.character(input$srtCol)),
						margins = c(input$bottom,input$right)
					)
					obj = recordPlot()
					return(obj)
				})
			})

			output$staticheatmap <- renderPlot({
				replayPlot(staticout())
			})

			observeEvent(input$staticheatmap, {
				saved_plots$staticheatmap <- staticout()
			})

			interactiveHeatmap <- eventReactive(input$action_heatmaps, {

				DataHeatMap <- DataHeatMapReactive()

				data.in <- DataHeatMap$df %>% as.data.frame()
				validate(need(nrow(data.in) <= 100, message = "Limited 100 genes."))

				annotation <- DataHeatMap$annotation
				cutree_rows = input$cutreerows
				cutree_cols = input$cutreecols
				if (cutree_rows == 0)
				cutree_rows = NULL
				if (cutree_cols == 0)
				cutree_cols = NULL


				if (input$dendrogram == "both" | input$dendrogram == "row")
				dend_r <- data.in %>% dist(method = input$distanceMethod) %>% hclust(method = input$agglomerationMethod) %>% as.dendrogram %>% ladderize %>%  dendextend::color_branches(k=cutree_rows)
				if (input$dendrogram == "both" | input$dendrogram == "column")
				dend_c <- t(data.in) %>% dist(method = input$distanceMethod) %>% hclust(method = input$agglomerationMethod) %>% as.dendrogram %>% ladderize %>% dendextend::color_branches(k=cutree_cols)

				cexRow = as.numeric(as.character(input$hyfontsizei))
				cexCol = as.numeric(as.character(input$hxfontsizei))

				labCol = colnames(data.in)
				labRow = rownames(data.in)


				if (cexRow  == 0 | nrow(data.in) > 50) {
					labRow = NA
					cexRow = 0.2
				}

				if (cexCol == 0) {
					labCol = NA
					cexCol  = 0.2
				}

				hide_colorbar=FALSE
				if (input$key == "FALSE")
				hide_colorbar=TRUE

				heatmaply::heatmaply(data.in#,
					#dendrogram = input$dendrogram,
					#colors = colorpanel (32, low = input$lowColor,mid = input$midColor, high = input$highColor),
					#Rowv = if (input$dendrogram == "both" | input$dendrogram == "row") dend_r else FALSE,
					#Colv = if (input$dendrogram == "both" | input$dendrogram == "column") dend_c else FALSE,
					#labRow = labRow,
					#labCol = labCol,
					#cexRow = cexRow,
					#cexCol = cexCol,
					#srtCol = as.numeric(as.character(input$srtCol)),
					#hide_colorbar = hide_colorbar
				) %>%
				layout(margin = list(l = input$l, b = input$b))
			})

			output$interactiveheatmap <- renderPlotly({
				withProgress(message = 'Making interactive heatmap:', value = 0, {
					interactiveHeatmap()
				})
			})
		}
	)
}
