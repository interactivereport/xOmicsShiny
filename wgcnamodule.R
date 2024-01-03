###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: network.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com); Kyra Griffin-Mitchell (kyra.griffinmitchell@Biogen.com)
##@Date : 02/23/2022
##@version 1.0
###########################################################################################################

##########################################################################################################
## WGCNA
##########################################################################################################
library(WGCNA)

wgcna_ui <- function(id) {
	ns <- shiny::NS(id)

	fluidRow(
		column(3,
			wellPanel(
				radioButtons(ns("WGCNAgenelable"),label="Select Gene Label",inline = TRUE, choices=c("id", "UniqueID", "Protein.ID", "Gene.Name", "Intensity"), selected="UniqueID"),
				sliderInput(ns("wgcna_rcut"), label= "Choose r Cutoff",  min = 0.7, max = 1, value = 0.9, step=0.02),
				# selectInput("wgcna_pcut", label= "Choose P Value Cutoff", choices= c("0.0001"=0.0001,"0.001"=0.001,"0.01"=0.01,"0.05"=0.05),selected=0.01),
				numericInput(ns("WGCNAtopNum"), label= "Top Number of Genes:",  value=250, min=250, step=25),
				numericInput(ns("minModuleSize"), label= "Mininum Module Size:",  value=30, min= 1, max = 1000),
				numericInput(ns("maxBlockSize"), label= "Max Block Size:",  value=4000, min = 100, max = 30000),
				actionButton(ns("plotwgcna"),"Generate")

			)
		),
		column(9,
			tabsetPanel(id="WGCNA_tabset",
				tabPanel(title="Dendrogram", value="Dendrogram",
					plotOutput(ns("Dendrogram"), height=800)
				),
				#tabPanel(title="Heatmap", value="Heatmap", uiOutput(ns("Heatmap"), style = "background-color: #eeeeee;")), #height="800px"
				#tabPanel(title="Adjacency Matrix", value="Adjacency Matrix", 	DT::dataTableOutput(ns("adj_WGCNA"))),
				tabPanel(title="Help", htmlOutput('help_WGCNA'))
			)
		)
	)
}

wgcna_server <- function(id) {
	shiny::moduleServer(id,
		function(input, output, session) {
			ns <- session$ns

			WGCNAReactive <- reactive({
				req(length(working_project()) > 0)
				req(DataInSets[[working_project()]]$data_wide)
				req(DataInSets[[working_project()]]$ProjectID)
				ProjectID <- DataInSets[[working_project()]]$ProjectID

				data_wide = DataInSets[[working_project()]]$data_wide
				ProteinGeneName  = DataInSets[[working_project()]]$ProteinGeneName

				wgcnafile <- paste("wgcnadata/", ProjectID, ".RData", sep = "")
				if (file.exists(wgcnafile)) {
					load(wgcnafile)
				} else {

					# Top number of genes
					topNum <- as.numeric(input$WGCNAtopNum)
					# Gene Label
					gene_label <- input$WGCNAgenelable

					dataExpr <- data_wide %>%
					na.omit()
					gene.names=rownames(dataExpr)
					dataExpr <- rownames_to_column(dataExpr, var = "UniqueID")
					dataExpr <- dplyr::left_join(dataExpr, ProteinGeneName)
					dataExpr$labelgeneid = dataExpr[,match(gene_label,colnames(dataExpr))]
					dataExpr <- dataExpr %>%
					dplyr::select(-c(id, UniqueID, Gene.Name, Protein.ID))

					dataExpr <- column_to_rownames(dataExpr, var = "labelgeneid")
					dataExpr[is.na(dataExpr) | dataExpr=="Inf"] = NA

					SubGeneNames=gene.names[1:topNum]
					
					# Ensure all columns are numeric before transposing; otherwise cell values may
					# become character, causing problems in WGCNA::blockwiseModules, as happened to
					# the Mouse_microglia_RNA-Seq data
					dataExpr <-  dataExpr %>%
					  dplyr::select(tidyselect::where(is.numeric))

					dataExpr = as.data.frame(t(dataExpr))
					dataExpr= dataExpr[,1:topNum]

					WGCNA::allowWGCNAThreads()
					ALLOW_WGCNA_THREADS=8
					enableWGCNAThreads()

					# Choose a set of soft-thresholding powers
					powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
					r2_cutoff <- input$wgcna_rcut

					cor <- WGCNA::cor
					sft <- WGCNA::pickSoftThreshold(dataExpr, dataIsExpr = TRUE, powerVector = powers,	corFnc = cor, corOptions = list(use = 'p'),	networkType = "unsigned")

					# Generating adjacency and TOM similarity matrices based on the selected softpower
					picked_power <- softPower <- sft$powerEstimate

					##calclute the adjacency matrix
					#adj= WGCNA::adjacency(dataExpr,type = "unsigned", power = softPower)
					#
					##turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
					#TOM=WGCNA::TOMsimilarityFromExpr(dataExpr,networkType = "unsigned", TOMType = "unsigned", power = softPower)
					#
					#colnames(TOM) = rownames(TOM) = SubGeneNames
					#dissTOM = 1 - TOM
					#
					##Module Detection
					##hierarchical clustering of the genes based on the TOM dissimilarity measure
					#geneTree = flashClust::flashClust(as.dist(dissTOM),method="average")
					#
					## #plot the resulting clustering tree (dendrogram)
					## plot(geneTree, xlab="", sub="",cex=0.35)
					##
					#
					## Set the minimum module size
					#minModuleSize = 20;
					#
					## Module identification using dynamic tree cut
					#dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
					#
					## #the following command gives the module labels and the size of each module.
					## #Lable 0 is reserved for unassigned genes
					## table(dynamicMods)
					#
					#dynamicColors = WGCNA::labels2colors(dynamicMods)
					## table(dynamicColors)
					##
					## plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
					##                     dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
					##                     guideHang = 0.05, main = "Gene dendrogram and module colors")
					#
					##discard the unassigned genes, and focus on the rest
					#restGenes= (dynamicColors != "grey")
					#diss1=1-WGCNA::TOMsimilarityFromExpr(dataExpr[,restGenes], power = softPower)
					#
					#colnames(diss1) =rownames(diss1) =SubGeneNames[restGenes]
					## hier1=flashClust(as.dist(diss1), method="average" )
					## plotDendroAndColors(hier1, dynamicColors[restGenes], "Dynamic Tree Cut",
					##                     dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
					##                     guideHang = 0.05, main = "Gene dendrogram and module colors")
					#
					##set the diagonal of the dissimilarity to NA
					#diag(diss1) = NA;
					##
					## #Visualize the Tom plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
					## sizeGrWindow(7,7)
					## TOMplot(diss1, hier1) #, as.character(dynamicColors[restGenes]))
					#
					## plot heatmap using plotly
					#plotly_heatmap <- plot_ly(z = diss1, type = "heatmap", colors = "YlOrRd")
					## heatmap <- d3heatmap(nba_players, scale = "column", color = "YlOrRd")
					## heat <- heatmap(diss1)
					#

					temp_cor <- cor
					cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
					netwk <- blockwiseModules(dataExpr,                # <= input here

						# == Adjacency Function ==
						power = picked_power,                # <= power here
						networkType = "signed",

						# == Tree and Block Options ==
						deepSplit = 2,
						pamRespectsDendro = F,
						# detectCutHeight = 0.75,
						minModuleSize = input$minModuleSize, #30,
						maxBlockSize = input$maxBlockSize,#4000,

						# == Module Adjustments ==
						reassignThreshold = 0,
						mergeCutHeight = 0.25,

						# == TOM == Archive the run results in TOM file (saves time)
						saveTOMs = F,
						saveTOMFileBase = "ER",

						# == Output Options
						numericLabels = T,
					verbose = 3)
					cor <- temp_cor
				}
				return(netwk)

			})

			observeEvent(input$plotwgcna,{
				output$Dendrogram <- renderPlot({
					netwk <-	WGCNAReactive()
					mergedColors = labels2colors(netwk$colors)

					plotDendroAndColors(
						netwk$dendrograms[[1]],
						mergedColors[netwk$blockGenes[[1]]],
						"Module colors",
						dendroLabels = FALSE,
						hang = 0.03,
						addGuide = TRUE,
					guideHang = 0.05 )
				})
			})
		}
	)
}
