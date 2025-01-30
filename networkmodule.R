###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: network.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 01/27/2025
##@version 3.0

##########################################################################################################
## Correlation Network
##########################################################################################################
#pkgs: "visNetwork" "networkD3", "DT",  "Hmisc","dplyr"

library(visNetwork)
library(networkD3)

network_ui <- function(id) {
	ns <- shiny::NS(id)
	fluidRow(
		column(3,
			wellPanel(
				uiOutput(ns('loadedprojects')),
				radioButtons(ns("label"),label="Select Gene by", inline = TRUE, choices=c("UniqueID", "Gene.Name"), selected="Gene.Name"),
				selectizeInput(ns("sel_gene"),	label="Gene Name (Select 1 or more)",	choices = NULL,	multiple=TRUE, options = list(placeholder =	'Type to search')),
				fluidRow(
					column(width=6, numericInput(ns("rcut"), label= "r Cutoff",  min = 0.7, max = 1, value = 0.8, step=0.05)),
					column(width=6, numericInput(ns("pvalcut"), label= "P Value Cutoff", value=0.01, min=0, step=0.001))
				),
				radioButtons(ns("psel"), label= "P value or P.adj Value?", inline = TRUE, choices= c("P.Value" = "p", "Adj.P.Value" = "Padj"), selected="Padj"),
				textOutput(ns("networkstat")),
				uiOutput(ns("myTabUI")),
				tags$hr(style="border-color: black;"),
				conditionalPanel(ns = ns, "input.tabset=='visNetwork'",
					tags$h4("Node Options:"),
					tags$h5("Label inside: ellipse, circle"),
					tags$h5("Label outside: diamond, star, triangle, triangleDown, hexagon, square."),
					fluidRow(
						column(width=6, colourpicker::colourInput(ns("nodecolor"), "Node Color", "#1E90FF", palette = "limited")),
						column(width=6, selectizeInput(ns("nodeshape"), label= "Node Shape", choices = c("ellipse"="ellipse","circle"="circle","diamond"="diamond","star"="star","triangle"="triangle","triangleDown"="triangleDown","hexagon"="hexagon","square"="square"),selected = "ellipse"))
					),
					sliderInput(ns("nodesize"), "Node Size", min = 10, max = 40, step = 5, value = 25),
					tags$h4("Edge Options:"),
					fluidRow(
						column(width=6, colourpicker::colourInput(ns("poscolor"), "Positive corr", "#FF0000", palette = "limited")),
						column(width=6, colourpicker::colourInput(ns("negcolor"), "Negative corr", "#1E90FF", palette = "limited"))
					),
					sliderInput(ns("edgewidth"), "Edge Width (relative to abs(r)):", min = 10, max = 50, step = 10, value = 20),
					tags$h4("Label Options:"),
					fluidRow(
						column(width=4, colourpicker::colourInput(ns("labelcolor"), "Label Color", "#FF0000", palette = "limited")),
						column(width=8, sliderInput(ns("fontsize"), "Label Font Size", min = 10, max = 30, step = 2, value = 16))
					),
					tags$h4("Format Options:"),
					radioButtons(ns("layout"), label= "Layout", choices= c("Random"="Random", "Hierarchical"="Hierarchical"), inline = TRUE, selected = "Random"),
					radioButtons(ns("highlighted"), label= "Highlighted Input", choices = c("no" = "no", "yes" = "yes"), inline = TRUE, selected = "no"),
					radioButtons(ns("freeze"), label= "Plot Freeze", choices = c("no" = "no", "yes" = "yes"), inline = TRUE, selected = "no")
			))
		),
		column(9,
			tabsetPanel(id=ns("tabset"),
				tabPanel(title="visNetwork", value ="visNetwork",
					actionButton(ns("plotvisnetwork"), "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
					visNetworkOutput(ns("visnetwork"), height="800px"), style = "background-color: #eeeeee;"),
				tabPanel(title="networkD3", value ="networkD3",
					actionButton(ns("plotnetworkD3"), "Plot/Refresh", style="color: #0961E3; background-color: #F6E98C ; border-color: #2e6da4"),
					forceNetworkOutput(ns("networkD3"), height="800px"), style = "background-color: #eeeeee;"),
				tabPanel(title="Selected Ids Data Table", DT::dataTableOutput(ns("dat_network"))),
				tabPanel(title="Correlation Summary", DT::dataTableOutput(ns("summary_network"))),
				tabPanel(title="Help", htmlOutput("help_network"))
			)
		)
	)
}

network_server <- function(id) {
	shiny::moduleServer(id,
		function(input, output, session) {
			ns <- session$ns
			output$loadedprojects <- renderUI({
					req(length(working_project()) > 0)
					radioButtons(ns("current_dataset"), label = "Change Working Dataset", choices=DS_names(), inline = F, selected=working_project())
			})

			observeEvent(input$current_dataset, {
					working_project(input$current_dataset)
			})

			DataNetworkReactive <- reactive({
					req(length(working_project()) > 0)
					req(DataInSets[[working_project()]]$ProteinGeneName)
					ProteinGeneName <- DataInSets[[working_project()]]$ProteinGeneName
					ProjectID <- DataInSets[[working_project()]]$ProjectID
		
					req(input$psel)
					req(input$pvalcut)
					
					run_network=FALSE
					CorResFile <- DataInSets[[working_project()]]$file2

					if (is.null(CorResFile))  {
						run_network=TRUE
						} else if (file.exists(CorResFile)) {
						load(CorResFile)
					} else {
						run_network=TRUE
					}

					if (run_network) {withProgress(message = 'Compute correlation network data.',	detail = 'This may take a few minutes...', value = 0,	{
								data_wide <- DataInSets[[working_project()]]$data_wide
								#if data_wide has many genes, trim down to 10K
								if (nrow(data_wide)>10000 ) {
									dataSD=apply(data_wide, 1, function(x) sd(x,na.rm=T))
									dataM=rowMeans(data_wide)
									diff=dataSD/(dataM+median(dataM))
									data_wide=data_wide[order(diff, decreasing=TRUE)[1:10000], ]
									cat("reduce gene size to 10K for project ", ProjectID, "\n")
								}
								cor_res <- Hmisc::rcorr(as.matrix(t(data_wide)))
								cormat <- cor_res$r
								pmat <- cor_res$P

								ut <- upper.tri(cormat)
								adj.p <- p.adjust(pmat[ut], method = "BH")

								network <- tibble (
									from = rownames(cormat)[row(cormat)[ut]],
									to = rownames(cormat)[col(cormat)[ut]],
									cor  = signif(cormat[ut], 2),
									p = signif(pmat[ut], 4),
									adj.p = adj.p,
									direction = as.integer(sign(cormat[ut]))
								)
								network <- network %>% mutate_if(is.factor, as.character) %>%
								dplyr::filter(!is.na(cor) & abs(cor) > 0.7 & p < 0.05)

								if (nrow(network)>2e6) {
									network <- network %>% mutate_if(is.factor, as.character) %>%
									dplyr::filter(!is.na(cor) & abs(cor) > 0.8 & p < 0.005)
								}
								if (nrow(network)>2e6) {
									network <- network %>% mutate_if(is.factor, as.character) %>%
									dplyr::filter(!is.na(cor) & abs(cor) > 0.85 & p < 0.005)
								}
								save(network,
									file =  paste("networkdata/", DataInSets[[working_project()]]$ProjectID, ".RData", sep = ""))
								DataInSets[[working_project()]]$file2=paste("networkdata/", DataInSets[[working_project()]]$ProjectID, ".RData", sep = "")
						})
					}

					if (!("Padj" %in% colnames(network))) { #for the previous preprocessed data
						network$Padj <- p.adjust(network$p, method = "BH") 
						network <- network %>%
						dplyr::relocate(Padj, .after = p)
						
					}

					rcutoff <- as.numeric(input$rcut)
					pvalcut <- as.numeric(input$pvalcut)
					pval_sel <- input$psel
					psel = sym(pval_sel)
					
					sel_gene = input$sel_gene
					tmpids = ProteinGeneName[unique(na.omit(c(apply(ProteinGeneName, 2, function(k) match(sel_gene, k))))), ]
					edges.sel <- network %>% dplyr::filter((from %in% tmpids$UniqueID) | (to %in% tmpids$UniqueID))

					edges <- dplyr::filter(edges.sel, abs(cor) > rcutoff & !!psel < pvalcut)
					networks_ids <-	unique(c(as.character(edges$from), as.character(edges$to)))
					nodes <- ProteinGeneName %>%
					dplyr::filter(UniqueID %in% networks_ids) %>%
					dplyr::select(UniqueID, Gene.Name) %>%
					dplyr::rename(id = UniqueID, label = Gene.Name)

					net <- list("network"=network, "nodes" = nodes, "edges" = edges)
					return(net)
			})

			observe({
					req(length(working_project()) > 0)
					req(DataInSets[[working_project()]]$ProteinGeneName)
					ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
					req(input$label)
					if (input$label=="UniqueID") {
						DataIngenes <- ProteinGeneName %>% dplyr::pull(UniqueID)
					} else {
						DataIngenes <- ProteinGeneName %>% dplyr::pull(Gene.Name)
					}
					updateSelectizeInput(session,'sel_gene', choices= DataIngenes, server=TRUE)
			})

			observe({
					net <-	DataNetworkReactive()
					output$networkstat <- renderText({
							sprintf("\nNodes:%d  Edges:%d",	nrow(net$nodes), nrow(net$edges))
					})

					if (nrow(net$nodes) == 0) {
						output$myTabUI <- renderUI({
								tags$h4("Zero node. Try lower cutoffs or select other genes.", style="color: red;")
						})
						} else if (nrow(net$nodes) > 500) {
						output$myTabUI <- renderUI({
								tags$h4("Too many nodes (limit to 500). Try higher cutoffs or select fewer genes.", style="color: red;")
						})
					} else {
						output$myTabUI <- renderUI({
								tags$h4("Generate plot by clicking plot/refersh button", style="color: green;")
						})
					}
			})

			observeEvent(input$plotvisnetwork,{
					output$visnetwork <- renderVisNetwork({
							withProgress(message = 'Making Network:', value = 0, {
									isolate({
											net <-	DataNetworkReactive()
											nodes <- net$nodes
											edges <- net$edges
											req(nrow(nodes)>0)

											inputgenes = input$sel_gene
											nodes <- nodes %>%
											dplyr::mutate(sel = ifelse(label %in% inputgenes, "input gene", ""))

											rcutoff <- as.numeric(input$rcut)
											nodesize <- as.numeric(input$nodesize)
											edgewidth <- as.numeric(input$edgewidth)
											edges$width <- (abs(edges$cor) - rcutoff)*edgewidth
											edges <- edges %>%
											dplyr::mutate(color = ifelse(direction == 1,input$poscolor, input$negcolor))

											pvis <- visNetwork(nodes, edges, height = "100%", width = "100%") %>%
											visNodes(color = list(background = input$nodecolor, border = input$nodecolor, highlight = "yellow"),
												shape = input$nodeshape, size =  nodesize,
												font = list(color = input$labelcolor, size = input$fontsize)) %>% visInteraction(navigationButtons = TRUE)

											if (input$layout == "Random"){
												pvis <- pvis %>% visLayout(randomSeed = 123)
											} else {
												pvis <- pvis %>% visHierarchicalLayout()
											}

											if (input$highlighted == "yes") {
												pvis <- pvis %>% visOptions(selectedBy = list(variable = "sel", selected = "input gene"))
											} else {
												pvis
											}

											if (input$freeze == "yes") {
												pvis <- pvis %>% visInteraction(dragNodes = FALSE, dragView = FALSE, zoomView = FALSE)
											} else {
												pvis
											}

									})
							})
					})
			})

			observeEvent(input$plotnetworkD3,{
					output$networkD3 <- renderForceNetwork({
							withProgress(message = 'Making Network:', value = 0, {
									isolate({
											net <-	DataNetworkReactive()
											nodes <- net$nodes
											edges <- net$edges
											req(nrow(nodes)>0)
											nodes$group = 1
											nodes$size = 10
											edgelist <- as.data.frame(edges)
											nodes <- as.data.frame(nodes)
											sources <- edgelist$from
											targets <- edgelist$to
											node_names <- factor(sort(unique(c(as.character(sources),  as.character(targets)))))
											links <- data.frame(source = match(sources, node_names) - 1,target = match(targets, node_names) - 1, value = edgelist$cor)
											nodes <- nodes[match(node_names, nodes$id),]
											forceNetwork(Links = links, Nodes = nodes, Source = "source",
												Target = "target", Value = "value", NodeID = "label",fontSize=7,
											Group = "group", opacity = 0.9,zoom = TRUE, opacityNoHover = 1)
									})
							})
					})
			})

			output$dat_network <- DT::renderDT(server=FALSE,{
					net <-	DataNetworkReactive()
					results <- as.data.frame(net$edges)
					results[,sapply(results,is.numeric)] <- signif(results[,sapply(results,is.numeric)],3)
					DT::datatable(results, extensions = 'Buttons', options = list(dom = 'lBfrtip', pageLength = 15,
							buttons = list(
								list(extend = "csv", text = "Download Page", filename = "Page_results",	exportOptions = list(modifier = list(page = "current"))),
								list(extend = "csv", text = "Download All", filename = "All_Results",	exportOptions = list(modifier = list(page = "all")))
							)
					))
			})

			output$summary_network <- DT::renderDT(server=FALSE, {
					req(length(working_project()) > 0)
					req(DataInSets[[working_project()]]$ProteinGeneName)
					ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
					req(input$label)
					
					net <-	DataNetworkReactive()
					network <- net$network
					
					network_switched <- network %>%
					dplyr::relocate(to, .before = from) %>%
					dplyr::rename(to = from,from = to)
					
					network <- network %>%
					dplyr::bind_rows(.,network_switched) %>%
					dplyr::distinct()
					rm(network_switched)
					
					from_count <- network %>%
					dplyr::select(from) %>%
					dplyr::group_by(from) %>%
					count() %>%
					dplyr::rename(UniqueID = from)
					
					to_count <- network %>%
					dplyr::select(to) %>%
					dplyr::group_by(to) %>%
					count() %>%
					dplyr::rename(UniqueID = to)
					
					correlation_count <- bind_rows(from_count, to_count) %>%
					dplyr::group_by(UniqueID) %>%
					dplyr::summarise(count=sum(n), .groups = 'drop') %>%
					dplyr::rename(Before_Filter = count)
					rm(from_count, to_count)
					
					pval_sel <- input$psel
					psel = sym(pval_sel)
					pvalcut <- as.numeric(input$pvalcut)
					rcutoff <- as.numeric(input$rcut)
					
					network_filtered <- dplyr::filter(network, abs(cor) > rcutoff & !!psel < pvalcut)
					
					from_count_filtered <- network_filtered %>%
					dplyr::select(from) %>%
					dplyr::group_by(from) %>%
					count() %>%
					dplyr::rename(UniqueID = from)
					
					to_count_filtered <- network_filtered %>%
					dplyr::select(to) %>%
					dplyr::group_by(to) %>%
					count() %>%
					dplyr::rename(UniqueID = to)
					
					correlation_count_filtered <- bind_rows(from_count_filtered, to_count_filtered) %>%
					dplyr::group_by(UniqueID) %>%
					dplyr::summarise(count=sum(n), .groups = 'drop') %>%
					dplyr::rename(After_Filter = count)
					rm(network_filtered, from_count_filtered, to_count_filtered)
					
					correlation_count <- correlation_count %>%
					dplyr::left_join(.,correlation_count_filtered, by="UniqueID") %>%
					dplyr::arrange(desc(After_Filter))
					rm(correlation_count_filtered)
					
					correlation_count <- dplyr::inner_join(ProteinGeneName,correlation_count, by = "UniqueID") %>%
					  dplyr::select(-any_of(c('id')))

					DT::datatable(correlation_count, extensions = 'Buttons', options = list(dom = 'lBfrtip', pageLength = 15,
							buttons = list(
								list(extend = "csv", text = "Download Page", filename = "Page_results",	exportOptions = list(modifier = list(page = "current"))),
								list(extend = "csv", text = "Download All", filename = "All_Results",	exportOptions = list(modifier = list(page = "all")))
							)
					))
			})

		}
	)
}
