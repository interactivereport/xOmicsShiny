###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: patternmodule.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 05/23/2023
##@version 3.0

##########################################################################################################
## Pattern Clustering
##########################################################################################################
#pkgs:  "Mfuzz","factoextra", "cluster", "shiny","DT", "dplyr","stringr", "tidyr"

library(Mfuzz)
#library(NbClust)
library(factoextra)
library(cluster)

pattern_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
		column(3,
			wellPanel(
				uiOutput(ns('loadedprojects')),
				uiOutput(ns("selectGroupSample")),
				radioButtons(ns("subset"),label="Use subset genes or upload your own subset?", choices=c("subset","upload genes","Geneset"), inline = TRUE, selected="subset"),
				conditionalPanel(ns = ns, "input.subset=='subset'",
					fluidRow(
						column(width=6,numericInput(ns("fccut"), label= "Choose Fold Change Threshold",value = 1.2, min=1, step=0.1)),
						column(width=6,numericInput(ns("pvalcut"), label= "Choose P-value Threshold", value=0.01, min=0, step=0.001))
					),
					radioButtons(ns("psel"), label= "P value or P.adj Value?", choices= c("Pval"="Pval","Padj"="Padj"), inline = TRUE),
					span(textOutput(ns("filteredgene")), style = "color:red; font-size:15px; font-family:arial; font-style:italic")
				),
				conditionalPanel(ns = ns, "input.subset=='upload genes'",
					textAreaInput(ns("list"), "list", "", cols = 5, rows=6)
				),
				conditionalPanel(ns = ns, "input.subset=='Geneset'",
					selectizeInput(ns("sel_geneset"), label="Available GeneSet", choices = NULL, multiple = FALSE),
					textAreaInput(ns("geneset_genes"), "Genes in Geneset", "", cols = 5, rows=6)
				),
				conditionalPanel("input.Pattern_tabset=='Clustering of Centroid Profiles'",
					radioButtons(ns("ClusterMehtod"), label="Cluster Method", inline = FALSE, choices = c("Soft Clustering" = "mfuzz", "K-means" = "kmeans", "Partitioning Around Medoids (take longer time)" = "pam")),
					fluidRow(
						column(width=6, sliderInput(ns("k"), "Cluster Number:", min = 2, max = 12, step = 1, value = 6)),
						column(width=6, sliderInput(ns("ncol"), label= "Column Number", min = 1, max = 6, step = 1, value = 3))
					),
					conditionalPanel(ns = ns,"input.ClusterMehtod=='kmeans'",
						fluidRow(
							column(width=6, sliderInput(ns("font"), "Font Size:", min = 12, max = 24, step = 1, value = 14)),
							column(width=6, sliderInput(ns("Xangle"), label= "X Angle", min = 0, max = 90, step = 15, value = 45))
						),
						radioButtons(ns("plot_Y_scale"), label="Y Axis Scale", inline = TRUE, choices = c("Auto","Manual"), selected = "Auto"),
						conditionalPanel(ns = ns, "input.plot_Y_scale=='Manual'",
							fluidRow(
								column(width=6, numericInput(ns("plot_Ymin"), label= "Y Min",  value = 0, step=0.1)),
								column(width=6, numericInput(ns("plot_Ymax"), label= "Y Max",  value=5, step=0.1))
							)
						)
					)
				),
				conditionalPanel("input.Pattern_tabset=='Optimal Number of Clusters'",
					radioButtons(ns("nbclustMehtod"), label="Method", inline = FALSE, choices = c("silhouette" = "silhouette", "wss" = "wss", "gap_stat (take longer time)" = "gap_stat")),
					sliderInput(ns("kmax"), "Max Cluster Number:", min = 2, max = 12, step = 1, value = 8),
				),
				conditionalPanel("input.Pattern_tabset=='Data Table'",
					radioButtons(ns("DataFormat"), label="Data Output Format:", inline = TRUE, choices = c("Wide Format" = "wide", "Long Format" = "long"))
				)
			)
		),
		column(9,
			tabsetPanel(id="Pattern_tabset",
				tabPanel(title="Clustering of Centroid Profiles",
					actionButton(ns("plot_pattern"), "Plot/Refresh"),
					actionButton(ns("pattern"), "Save to output"),
					plotOutput(ns("pattern"), height=800)
				),
				tabPanel(title="Optimal Number of Clusters",
					plotOutput(ns("nbclust"), height=800)
				),
				#tabPanel(title="Data Table",actionButton("Pattern_data", "Save to output"), DT::dataTableOutput("dat_pattern")),
				tabPanel(title="Data Table",
					DT::dataTableOutput(ns("dat_pattern"))
				),
				tabPanel(title="Help", htmlOutput("help_pattern"))
			)
		)
	)
}

pattern_server <- function(id) {
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
				req(DataInSets[[working_project()]]$tests_order)
				tests = c("ALL", DataInSets[[working_project()]]$tests_order)
				allgroups = DataInSets[[working_project()]]$groups
				groups = DataInSets[[working_project()]]$group_order
				updateSelectizeInput(session,'group', choices=allgroups, selected=groups)
				updateSelectizeInput(session,'sel_test', choices=tests, selected=tests[1])
			})

			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$ProteinGeneName)
				req(DataInSets[[working_project()]]$results_long)
				ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
				results_long = DataInSets[[working_project()]]$results_long
				fccut = log2(as.numeric(input$fccut))
				pvalcut = as.numeric(input$pvalcut)

				req(input$psel)
				if (input$psel == "Padj") {
					filteredgene1 = results_long %>%
					dplyr::filter(abs(logFC) > fccut & Adj.P.Value < pvalcut) %>%
					dplyr::pull(UniqueID)
				} else {
					filteredgene1 = results_long %>%
					dplyr::filter(abs(logFC) > fccut & P.Value < pvalcut) %>%
					dplyr::pull(UniqueID)
				}
				output$filteredgene <- renderText({paste("Selected Genes:",length(filteredgene1),sep="")})
			})

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

			DatapatternReactive <- reactive({
				req(length(working_project()) > 0)
				ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
				sample_group <- DataInSets[[working_project()]]$sample_group
				fccut = log2(as.numeric(input$fccut))
				pvalcut = as.numeric(input$pvalcut)
				sel_group = DataInSets[[working_project()]]$group_order
				results_long <- DataInSets[[working_project()]]$results_long
				data_long <- DataInSets[[working_project()]]$data_long

				if (input$subset == "subset") {
					if (input$psel == "Padj") {
						filteredgene = results_long %>%
						dplyr::filter(abs(logFC) > fccut & Adj.P.Value < pvalcut) %>%
						dplyr::pull(UniqueID)
					} else {
						filteredgene = results_long %>%
						dplyr::filter(abs(logFC) > fccut & P.Value < pvalcut) %>%
						dplyr::pull(UniqueID)
					}
				}

				if (input$subset=="Upload" | input$subset=="Geneset") {
					if (input$subset=="Upload") {
						req(input$gene_list)
						gene_list <- input$gene_list
					} else {
						req(input$geneset_genes)
						gene_list <- input$geneset_genes
					}

					gene_list <- ProcessUploadGeneList(gene_list)

					validate(need(length(gene_list)>2, message = "input gene list"))
					filteredgene <- dplyr::filter(ProteinGeneName, (UniqueID %in% gene_list) | (Protein.ID %in% gene_list) | (toupper(Gene.Name) %in% toupper(gene_list)))  %>%
					dplyr::pull(UniqueID)
				}

				subdatlong <- dplyr::filter(data_long, (group %in% sel_group) & (UniqueID %in% filteredgene)) %>%
				group_by(., group, UniqueID) %>%
				dplyr::summarise(mean=mean(expr, na.rm = TRUE), .groups = 'drop')

				subdatwide <- subdatlong  %>%
				tidyr::spread(.,group, mean, fill = 0) %>% as.data.frame() %>%
				remove_rownames(.) %>%
				column_to_rownames(.,var="UniqueID") %>%
				dplyr::select(all_of(sel_group))

				return(list("subdatlong"= subdatlong,"subdatwide"= subdatwide, "filteredgene" = filteredgene))
			})

			out <- eventReactive(input$plot_pattern, {withProgress(message = 'Processing...', value = 0, {

				Datapattern <- DatapatternReactive()
				subdatwide <- Datapattern$subdatwide
				subdatlong <- Datapattern$subdatlong
				sel_group <- DataInSets[[working_project()]]$group_order

				k=input$k
				set.seed(123)
				if (input$ClusterMehtod == "kmeans") {

					cl <- kmeans(subdatwide, k)
					cluster<-cl$cluster

					cluster.df <- data.frame(UniqueID=names(cluster), cluster=paste('Cluster',cluster,sep=' '), row.names=NULL)

					subdatlong <- subdatlong  %>%
					left_join(., cluster.df, by="UniqueID")

					subdatlong$group = factor(subdatlong$group,levels = sel_group)

					p <- ggplot(subdatlong, aes(x=group, y=mean)) +
					facet_wrap(~ cluster,scales = "free", ncol = input$ncol) +
					geom_line(aes(group=UniqueID, color="UniqueID")) +
					stat_summary(aes(color="red", group=1), fun=mean, geom="line", size=1.2, group=1)

					if (input$plot_Y_scale=="Manual") {
						p <- p + ylim(input$plot_Ymin, input$plot_Ymax)
					}

					p <- p +	theme_bw(base_size = input$font) + ylab("expr") + xlab(" ") +
					theme (plot.margin = unit(c(1,1,1,1), "cm"), axis.text.x = element_text(angle = input$Xangle),legend.position="none")
					p
				} else if (input$ClusterMehtod == "pam") {
					clpam <- cluster::pam(subdatwide, k)
					cluster <- clpam$clustering

					cluster.df <- data.frame(UniqueID=names(cluster), cluster=cluster, row.names=NULL)
					subdatlong <- subdatlong  %>%
					left_join(., cluster.df, by="UniqueID")
					subdatlong$group = factor(subdatlong$group,levels = sel_group)

					p <- ggplot(subdatlong, aes(x=group, y=mean)) +
					facet_wrap(~ cluster,scales = "free", ncol = 3) +
					geom_line(aes(group=UniqueID, color="UniqueID")) +
					stat_summary(aes(color="red", group=1), fun=mean, geom="line", size=1.2, group=1)
					p <- p + theme_bw(base_size = 14) + ylab("expr") + xlab(" ") +
					theme (plot.margin = unit(c(1,1,1,1), "cm"), axis.text.x = element_text(angle = 45),legend.position="none")
					return(p)

				} else if (input$ClusterMehtod == "mfuzz") {
					tmp_expr <- new('ExpressionSet', exprs = as.matrix(subdatwide))
					m1 <- mestimate(tmp_expr)
					cl <- mfuzz(tmp_expr, c = k, m = m1, iter.max = 200)
					nrow=ceiling(k/input$ncol)
					mfuzz.plot(tmp_expr, cl = cl, mfrow = c(nrow, input$ncol), min.mem=0.4, time.labels=colnames(subdatwide),new.window = FALSE)
					p = recordPlot()
					return(p)
				}
			})
		})

		output$pattern<- renderPlot({
			if (input$ClusterMehtod == "kmeans") {
				out()
			} else if (input$ClusterMehtod == "pam") {
				out()
			} else if (input$ClusterMehtod == "mfuzz") {
				replayPlot(out())
			}
		})

		observeEvent(input$pattern, {
			if (input$ClusterMehtod == "kmeans") {
				saved_plots$patternkmeans <- out()
			} else if (input$ClusterMehtod == "pam") {
				saved_plots$patternpam <- out()
			} else if (input$ClusterMehtod == "mfuzz") {
				saved_plots$patternmfuzz <- out()
			}
		})


		output$dat_pattern <- DT::renderDataTable({withProgress(message = 'Processing...', value = 0, {
			set.seed(123)
			Datapattern <- DatapatternReactive ()
			subdatwide <- Datapattern$subdatwide
			subdatlong <- Datapattern$subdatlong
			sel_group <- DataInSets[[working_project()]]$group_order
			k=input$k
			if (input$ClusterMehtod == "kmeans") {
				cl <- kmeans(subdatwide, k)
				cluster <- cl$cluster
				cluster.df <- data.frame(UniqueID=names(cluster), cluster=cluster, row.names=NULL)
			} else if (input$ClusterMehtod == "pam") {
				clpam <- pam(subdatwide, k)
				cluster <- clpam$clustering
				cluster.df <- data.frame(UniqueID=names(cluster), cluster=cluster, row.names=NULL)
			} else if (input$ClusterMehtod == "mfuzz") {
				tmp_expr <- new('ExpressionSet', exprs = as.matrix(subdatwide))
				m1 <- mestimate(tmp_expr)
				cl <- mfuzz(tmp_expr, c = k, m = m1, iter.max = 200)
				cluster <- cl$cluster
				cluster.df <- data.frame(UniqueID=names(cluster), cluster=cluster, row.names=NULL)
			}

			if (input$DataFormat == "long") {
				subdatlong[,sapply(subdatlong,is.numeric)] <- signif(subdatlong[,sapply(subdatlong,is.numeric)],3)
				subdatlong <- subdatlong  %>%
				left_join(., cluster.df, by="UniqueID")
				DT::datatable(subdatlong)
			} else if (input$DataFormat == "wide") {
				subdatwide[,sapply(subdatwide,is.numeric)] <- signif(subdatwide[,sapply(subdatwide,is.numeric)],3)
				subdatwide  <- subdatwide %>%
				rownames_to_column(.,var="UniqueID") %>%
				left_join(., cluster.df, by="UniqueID")%>%
				separate(UniqueID, c("Gene", "ID"), sep = "_")

				DT::datatable(subdatwide,	extensions = 'Buttons',
					options = list(dom = "Blfrtip",	buttons = list("copy", list(extend = "collection",buttons = c("csv", "excel", "pdf"),	text = "Download")),
					lengthMenu = list( c(10, 20, -1), c(10, 20, "All")), pageLength = 10),
					filter = 'top')
				}
			})
		})

		output$nbclust <- renderPlot({withProgress(message = 'Processing...', value = 0, {
			Datapattern <- DatapatternReactive()
			subdatwide <- Datapattern$subdatwide
			factoextra::fviz_nbclust(subdatwide, kmeans, method = input$nbclustMehtod, k.max = input$kmax, nboot=20) + theme_bw(base_size = 14)
		})
	})
}
)
}

