###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: network.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 05/23/2023
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
					radioButtons(ns("label"),label="Select Gene Label",inline = TRUE, choices=c("UniqueID", "Gene.Name"), selected="Gene.Name"),
					selectizeInput(ns("sel_net_gene"),	label="Gene Name (Select 1 or more)",	choices = NULL,	multiple=TRUE, options = list(placeholder =	'Type to search')),
					sliderInput(ns("rcut"), label= "Choose r Cutoff",  min = 0.5, max = 1, value = 0.9, step=0.02),
					selectInput(ns("pcut"), label= "Choose P Value Cutoff", choices= c("0.0001"=0.0001,"0.001"=0.001,"0.01"=0.01,"0.05"=0.05),selected=0.01),
					textOutput(ns("networkstat")),
					uiOutput(ns("myTabUI"))
				)
			),
			column(9,
				tabsetPanel(id="Network_tabset",
					tabPanel("visNetwork", visNetworkOutput(ns("visnetwork"), height="800px"), style = "background-color: #eeeeee;"),
					tabPanel("networkD3", forceNetworkOutput(ns("networkD3"), height="800px"), style = "background-color: #eeeeee;"),
					tabPanel(title="Data Table",	DT::dataTableOutput(ns("dat_network"))),
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
					if (nrow(data_wide)>3000 ) {
						dataSD=apply(data_wide, 1, function(x) sd(x,na.rm=T))
						dataM=rowMeans(data_wide)
						diff=dataSD/(dataM+median(dataM))
						data_wide=data_wide[order(diff, decreasing=TRUE)[1:10000], ]
						cat("reduce gene size to 10K for project ", ProjectID, "\n")
					}
					cor_res <- Hmisc::rcorr(as.matrix(t(data_wide)))
					cormat <- cor_res$r
					pmat <- cor_res$P
				#	cormat <- coop::pcor(as.matrix(t(data_wide)))

					ut <- upper.tri(cormat)
					network <- tibble (
						from = rownames(cormat)[row(cormat)[ut]],
						to = rownames(cormat)[col(cormat)[ut]],
						cor  = signif(cormat[ut], 2),
						p = signif(pmat[ut], 2),
						direction = as.integer(sign(cormat[ut]))
					)
					network <- network %>% mutate_if(is.factor, as.character) %>%
				dplyr::filter(!is.na(cor) & abs(cor) > 0.7 & p < 0.05)
				#	dplyr::filter(!is.na(cor) & abs(cor) > 0.7)

					if (nrow(network)>2e6) {
						network <- network %>% mutate_if(is.factor, as.character) %>%
						dplyr::filter(!is.na(cor) & abs(cor) > 0.8 & p < 0.005)
						#dplyr::filter(!is.na(cor) & abs(cor) > 0.8)
					}
					if (nrow(network)>2e6) {
						network <- network %>% mutate_if(is.factor, as.character) %>%
						dplyr::filter(!is.na(cor) & abs(cor) > 0.85 & p < 0.005)
						#dplyr::filter(!is.na(cor) & abs(cor) > 0.85)
					}
					save(network,
					file =  paste("networkdata/", DataInSets[[working_project()]]$ProjectID, ".RData", sep = ""))
						DataInSets[[working_project()]]$file2=paste("networkdata/", DataInSets[[working_project()]]$ProjectID, ".RData", sep = "")
					})
				}

				sel_gene = input$sel_net_gene
				tmpids = ProteinGeneName[unique(na.omit(c(apply(ProteinGeneName, 2, function(k) match(sel_gene, k))))), ]

				edges.sel <- network %>% dplyr::filter((from %in% tmpids$UniqueID) |	(to %in% tmpids$UniqueID))
				rcutoff <- as.numeric(input$rcut)
				pvalcutoff <- as.numeric(as.character(input$pcut))
				#edges <- dplyr::filter(edges.sel, abs(cor) > rcutoff & p < pvalcutoff)
				edges <- dplyr::filter(edges.sel, abs(cor) > rcutoff)
				networks_ids <-	unique(c(as.character(edges$from), as.character(edges$to)))
				nodes <- ProteinGeneName %>%
				dplyr::filter(UniqueID %in% networks_ids) %>%
				dplyr::select(UniqueID, Gene.Name) %>%
				dplyr::rename(id = UniqueID, label = Gene.Name)

				net <- list("nodes" = nodes, "edges" = edges)
				return(net)
			})

			observe({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$ProteinGeneName)
				ProteinGeneName = DataInSets[[working_project()]]$ProteinGeneName
				req(input$label)
				if (input$label=="UniqueID") {
					#DataIngenes <- ProteinGeneName %>% dplyr::select(UniqueID) %>% collect %>% .[["UniqueID"]] %>%	as.character()
					DataIngenes <- ProteinGeneName %>% dplyr::pull(UniqueID) 
				} else {
					#DataIngenes <- ProteinGeneName %>% dplyr::select(Gene.Name) %>% collect %>% .[["Gene.Name"]] %>%	as.character()
					DataIngenes <- ProteinGeneName %>% dplyr::pull(Gene.Name) 
				}
				updateSelectizeInput(session,'sel_net_gene', choices= DataIngenes, server=TRUE)
			})


			observe({
				net <-	DataNetworkReactive()
				output$networkstat <- renderText({
					sprintf("\nNodes:%d  Edges:%d",	nrow(net$nodes), nrow(net$edges))
				})

				if (nrow(net$nodes) > 0 & nrow(net$nodes) < 200){
					output$myTabUI <- renderUI({
						actionButton(ns("gennet"),"Generate")
					})
				} else if (nrow(net$nodes) == 0) {
					output$myTabUI <- renderUI({
						"Zero node. Try lower cutoffs or select other genes."
					})
				} else {
					output$myTabUI <- renderUI({
						"Too many nodes. Try higher cutoffs or select fewer genes."
					})
				}
			})

			observeEvent(input$gennet,{
				output$visnetwork <- renderVisNetwork({
					withProgress(message = 'Making Network:', value = 0, {
						isolate({
							net <-	DataNetworkReactive()
							visNetwork(net$nodes,net$edges,  height = "800px", width = "100%") %>%
							visLayout(randomSeed = 123) %>%
							visInteraction(navigationButtons = TRUE)

						})
					})
				})
			})

			observeEvent(input$gennet,{
				output$networkD3 <- renderForceNetwork({
					withProgress(message = 'Making Network:', value = 0, {
						isolate({
							net <-	DataNetworkReactive()
							net$nodes$group = 1
							net$nodes$size = 10
							edgelist <- as.data.frame(net$edges)
							nodes <- as.data.frame(net$nodes)
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
				DataIn <- DataReactive()
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
		}
	)
}
