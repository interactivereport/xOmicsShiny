###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: network.R
##@Developer : Lin Tinchi(tinchi.lin@biogen.com); Benbo Gao (benbo.gao@Biogen.com); Kyra Griffin-Mitchell (kyra.griffinmitchell@Biogen.com)
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
	  rclipboard::rclipboardSetup(),
		column(3,
			wellPanel(
			  uiOutput(ns('loadedprojects')),
				radioButtons(ns("WGCNAgenelable"),label="Select Gene Label",inline = TRUE, choices=c("Gene.Name","UniqueID"), selected="Gene.Name"),
				#sliderInput(ns("wgcna_rcut"), label= "R-Squared Cutoff for Picking Soft-threshold Power",  min = 0.7, max = 1, value = 0.9, step=0.02),
				# selectInput("wgcna_pcut", label= "Choose P Value Cutoff", choices= c("0.0001"=0.0001,"0.001"=0.001,"0.01"=0.01,"0.05"=0.05),selected=0.01),
				numericInput(ns("WGCNAtopNum"), label= "Select Top N Genes, where N is :",  value=250L, min=250L, step=25L, max = 10000L),
				numericInput(ns("mergeCutHeight"), label= "Dendrogram Cut Height for Merging:",  value=0.25, min= 0, max = 1.0, step = 0.01),
				#numericInput(ns("minModuleSize"), label= "Mininum Module Size:",  value=30L, min= 1L, max = 1000L),
				#numericInput(ns("maxBlockSize"), label= "Max Block Size:",  value=4000, min = 100, max = 30000),
				actionButton(ns("plotwgcna"),"Re-run"),
				br(),
				span("1. If the data is one of the saved projects in the CSV file,the app will load precomputed results based on up to 10,000 genes with default parameter values." ,style="color:red", inline = TRUE),
				br(),
				span("2. If you wish to re-run WGCNA on the saved project with a different parameter or number of genes, please click 'Re-run' but refrain from clicking it repeatedly.",style="color:red", inline = TRUE)

			)
		),
		column(9,
		  
			tabsetPanel(id="WGCNA_tabset",
				tabPanel(title="Dendrogram", value="Dendrogram",
					plotOutput(ns("Dendrogram"), height=800)
					#uiOutput(NS(id, 'dendro_container_ui'))
				),
				#tabPanel(title="Heatmap", value="Heatmap", uiOutput(ns("Heatmap"), style = "background-color: #eeeeee;")), #height="800px"
				#tabPanel(title="Adjacency Matrix", value="Adjacency Matrix", 	DT::dataTableOutput(ns("adj_WGCNA"))),
				tabPanel(title="Gene Clusters", DT::dataTableOutput(ns("gene_cluster"))),
				tabPanel(title="Help", htmlOutput('help_WGCNA'))
			)
		)
	)
}

wgcna_server <- function(id) {
	shiny::moduleServer(id,
		function(input, output, session) {
			ns <- session$ns
      
			output$loadedprojects <- renderUI({
			  req(length(working_project()) > 0)
			  radioButtons(ns("current_dataset"), label = "Change Working Dataset", choices=DS_names(), inline = F, selected=working_project())
			})
      
			toListen <- reactive({
			  req(input$current_dataset) 
			  req(input$WGCNAgenelable)
			})
			
			observeEvent(toListen(), {

			  req(length(working_project()) > 0)
			  req(DataInSets[[working_project()]]$data_wide)
			  req(DataInSets[[working_project()]]$ProjectID)
			  req(DataInSets[[working_project()]]$ProteinGeneName)
			  
			  working_project(input$current_dataset)
			  
			  data_wide <- DataInSets[[working_project()]]$data_wide
			  
			  if (nrow(data_wide)>10000 ) {
			    data_wide <- data_wide %>% na.omit()
			    dataSD=apply(data_wide, 1, function(x) sd(x,na.rm=T))
			    dataM=rowMeans(data_wide)
			    diff=dataSD/(dataM+median(dataM))
			    data_wide=data_wide[order(diff, decreasing=TRUE)[1:10000], ] 
			    dataExpr <- data_wide
			  } else {
			    dataExpr <- data_wide %>%
			      na.omit()
			  }
			  
			  default_n_gene <- min(10000, nrow(dataExpr))
			  
			  updateNumericInput(session, "WGCNAtopNum", 
			                     label= "Select Top N Genes, where N is :",  value=default_n_gene, min=250L, step=25L, max = default_n_gene)
			  
			  ProjectID <- DataInSets[[working_project()]]$ProjectID
			  
			  wgcnafile <- paste("data/wgcna_data/wgcna_", ProjectID, ".RData", sep = "")
			  
			  load(wgcnafile)
			  wgcna <- netwk
			  
			  mergedColors = labels2colors(wgcna$colors)
			  
			  output$Dendrogram <- renderPlot({
			    withProgress(message = "Creating plot using pre-calculated data", value = 0, {
  			    plotDendroAndColors(
  			      wgcna$dendrograms[[1]],
  			      mergedColors[wgcna$blockGenes[[1]]],
  			      "Module colors",
  			      dendroLabels = FALSE,
  			      hang = 0.03,
  			      addGuide = TRUE,
  			      guideHang = 0.05 )
			    })
			  })
			  
			  ProteinGeneName  <- DataInSets[[working_project()]]$ProteinGeneName
			  gene_label <- input$WGCNAgenelable
			  
			  # t0: merge WGCNA output with ProteinGeneName so that genes can be shown as UniqueID or Gene.Name
			  t0 <- tibble::tibble(UniqueID = names(wgcna$colors), color = labels2colors(wgcna$colors)) %>%
			    dplyr::left_join(ProteinGeneName[, c("UniqueID","Gene.Name")], by = "UniqueID") %>%
			    dplyr::select(color,all_of(gene_label)) %>%
			    dplyr::rename(gene = gene_label)
			  
			  # t1: collapse all genes in a cluster into a cell
			  t1 <- t0 %>%
			    dplyr::group_by(color) %>% 
			    dplyr::summarize(n_gene = n(),
			                     gene_group = paste0(gene, collapse = ",")) %>%
			    dplyr::ungroup()
			  
			  # t2: add the copy button
			  t2 <- t1
			  t2$copy <- vapply(1L:nrow(t1), function(i){
			    as.character(
			      rclipButton(
			        paste0("clipbtn_", i), 
			        label = "Copy all genes in cluster", 
			        clipText = t1[i, "gene_group"], 
			        #icon = icon("clipboard"),
			        icon = icon("copy", lib = "glyphicon"),
			        class = "btn-primary btn-sm"
			      )
			    )
			  }, character(1L))
			  
			  # rearrange columns
			  t2 <- t2 %>% dplyr::select(color, n_gene, copy, gene_group)
			  
			  output$gene_cluster <- DT::renderDT({
			    DT::datatable(
			      t2,
			      escape = FALSE,
			      selection = "none",
			      colnames=c("Color of cluster", "Number of genes", "Action","Genes in cluster")
			    )
			  })
			})
			
			# use eventReactive to control reactivity of WGCNAReactive;
			# otherwise, whenever an input change, WGCNAReactive will be re-calculated
			# and its re-calculation could take a long time.
			WGCNAReactive <- eventReactive(input$plotwgcna, {
			  withProgress(message = "Running WGCNA", detail = 'This may take a while...', value = 0.2, {
			    
			    # what if the user-imported data doesn't have $data_wide, $ProjectID..etc?
  			  req(length(working_project()) > 0)
  			  req(DataInSets[[working_project()]]$data_wide)
  			  req(DataInSets[[working_project()]]$ProjectID)
  			  req(DataInSets[[working_project()]]$ProteinGeneName)
  			  ProjectID <- DataInSets[[working_project()]]$ProjectID
  			  
  			  data_wide = DataInSets[[working_project()]]$data_wide
  			  ProteinGeneName  = DataInSets[[working_project()]]$ProteinGeneName
  			  
  			  if (nrow(data_wide)>10000 ) {
  			    data_wide <- data_wide %>% na.omit()
  			    dataSD=apply(data_wide, 1, function(x) sd(x,na.rm=T))
  			    dataM=rowMeans(data_wide)
  			    diff=dataSD/(dataM+median(dataM))
  			    data_wide=data_wide[order(diff, decreasing=TRUE)[1:10000], ] 
  			    dataExpr <- data_wide
  			    cat("reduce gene size to 10K for project ", ProjectID, "\n")
  			  } else {
  			    dataExpr <- data_wide %>%
  			      na.omit()
  			  }
  			  print(paste0("**** dim of dataExpr after-preprocssin is ****", dim(dataExpr)))
  			  
  			  # Note: if launching app from the server, the path for `load_`  files should be
  			  # paste0("/mnt/depts/dept04/compbio/projects/xOmicsShiny/data/wgcna_data/TOM
  			  load_wgcna_file <- paste("data/wgcna_data/load_", ProjectID, ".RData", sep = "")
  			  
  			  default_n_gene <- min(10000, nrow(dataExpr))
  			  
  			  if (file.exists(load_wgcna_file) & default_n_gene==input$WGCNAtopNum){
  			  
  			    # If file exist and the number of genes selected rename the same, load 
  			    # pre-computed result and TOM file (blockwiseModules(loadTom = T)) to 
  			    # reduce running time
  			    
  			    # The load_*.RData contains two objects, dataExpr and picked_power, so that
  			    # the app doesn't need to recalculate either from scratch
  			     load(load_wgcna_file)
  			    
  			    print(paste0("**** scenario 1 ****"))
 
  			    WGCNA::allowWGCNAThreads()
  			    ALLOW_WGCNA_THREADS=8L
  			    enableWGCNAThreads() # this causes much longer time if app launch from local machine, but not so from server
  			    cor <- WGCNA::cor
  			    
  			    temp_cor <- cor
  			    cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
  			    netwk <- blockwiseModules(dataExpr,                # <= input here
  			                              
  			                              # == Adjacency Function ==
  			                              power = picked_power,                # <= power here
  			                              networkType = "signed",
  			                              
  			                              # == Tree and Block Options ==
  			                              deepSplit = 2L,
  			                              pamRespectsDendro = F,
  			                              # detectCutHeight = 0.75,
  			                              minModuleSize = min(20, ncol(dataExpr/2)), # al# 30, #input$minModuleSize, #30,
  			                              # set block size to be number of genes, so that all
  			                              # genes will be analyzed in a single block
  			                              maxBlockSize = input$WGCNAtopNum,
  			                              
  			                              # == Module Adjustments ==
  			                              reassignThreshold = 0,
  			                              mergeCutHeight = input$mergeCutHeight,#,0.25,
  			                              
  			                              # == TOM == Archive the run results in TOM file (saves time)
  			                              saveTOMs = F,
  			                              loadTOM = TRUE,
  			                              
  			                              # Note: When launching from server, the path for TOM should be
  			                              # paste0("/mnt/depts/dept04/compbio/projects/xOmicsShiny/data/wgcna_data/TOM_",x)
  			                              saveTOMFileBase = paste0("./data/wgcna_data/TOM_", ProjectID),
  			                              
  			                              # == Output Options
  			                              numericLabels = T,
  			                              verbose = 3L)
  			    cor <- temp_cor
  			    
  			  } else if (file.exists(load_wgcna_file) & (default_n_gene - input$WGCNAtopNum)/default_n_gene < 0.1) {
  			    
  			    # If file exist and the number of genes selected is within 10% of 
  			    # the default number of genes, load pre-computed result
  			    # but do not load TOM file (blockwiseModules(loadTom = F))
  			    
  			    load(load_wgcna_file)
  			    
  			    print(paste0("**** scenario 2 ****"))
  			    
  			    dataExpr= dataExpr[,1L:input$WGCNAtopNum]
  			    
  			    WGCNA::allowWGCNAThreads()
  			    ALLOW_WGCNA_THREADS=8L
  			    enableWGCNAThreads() # this causes much longer time if app launch from local machine, but not so from server
  			    cor <- WGCNA::cor
  			    
  			    temp_cor <- cor
  			    cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
  			    netwk <- blockwiseModules(dataExpr,                # <= input here
  			                              
  			                              # == Adjacency Function ==
  			                              power = picked_power,                # <= power here
  			                              networkType = "signed",
  			                              
  			                              # == Tree and Block Options ==
  			                              deepSplit = 2L,
  			                              pamRespectsDendro = F,
  			                              # detectCutHeight = 0.75,
  			                              minModuleSize = min(20, ncol(dataExpr/2)), # al# 30, #input$minModuleSize, #30,
  			                              # set block size to be number of genes, so that all
  			                              # genes will be analyzed in a single block
  			                              maxBlockSize = input$WGCNAtopNum,
  			                              
  			                              # == Module Adjustments ==
  			                              reassignThreshold = 0,
  			                              mergeCutHeight = input$mergeCutHeight,#,0.25,
  			                              
  			                              # == TOM == Archive the run results in TOM file (saves time)
  			                              saveTOMs = F,
  			                              loadTOM = FALSE,

  			                              # Note: When launching from server, the path for TOM should be
  			                              # paste0("/mnt/depts/dept04/compbio/projects/xOmicsShiny/data/wgcna_data/TOM_",x)
  			                              saveTOMFileBase = paste0("./data/wgcna_data/TOM_", ProjectID),
  			                              
  			                              # == Output Options
  			                              numericLabels = T,
  			                              verbose = 3L)
  			    cor <- temp_cor

  			  } else {
  			    
  			    print(paste0("**** compute everything from scratch ****"))
  			    
  			    ## Top number of genes
  			    topNum <- as.numeric(input$WGCNAtopNum)
  			    # Gene Label
  			    gene_label <- input$WGCNAgenelable
  			    
  			    dataExpr <- data_wide %>%
  			      na.omit()
  			    gene.names=rownames(dataExpr)
  			    dataExpr <- rownames_to_column(dataExpr, var = "UniqueID")
  			    dataExpr <- dplyr::left_join(dataExpr, ProteinGeneName)
  			    # row names of dataExpr should be unique; allowing for other non-unique
  			    # vector would cause problems
  			    #dataExpr$labelgeneid = dataExpr[,match(gene_label,colnames(dataExpr))]
  			    dataExpr <- dataExpr %>%
  			      dplyr::select(-c(id, Gene.Name, Protein.ID))
  			    
  			    # row names of dataExpr should be unique; allowing for other non-unique
  			    dataExpr <- column_to_rownames(dataExpr, var = "UniqueID")
  			    dataExpr[is.na(dataExpr) | dataExpr=="Inf"] = NA
  			    
  			    SubGeneNames=gene.names[1L:topNum]
  			    
  			    # Ensure all columns are numeric before transposing; otherwise cell values may
  			    # become character, causing problems in WGCNA::blockwiseModules, as happened to
  			    # the Mouse_microglia_RNA-Seq data
  			    dataExpr <-  dataExpr %>%
  			      dplyr::select(tidyselect::where(is.numeric))
  			    
  			    dataExpr = as.data.frame(t(dataExpr))
  			    dataExpr= dataExpr[,1L:topNum]
  			    
  			    WGCNA::allowWGCNAThreads()
  			    ALLOW_WGCNA_THREADS=8L
  			    enableWGCNAThreads() # this causes much longer time if app launch from local machine, but not so from server
  			    
  			    # Choose a set of soft-thresholding powers
  			    powers <- c(c(1L:10L), seq(from = 12L, to = 20L, by = 2L))
  			    #r2_cutoff <- input$wgcna_rcut
  			    
  			    cor <- WGCNA::cor
  			    sft <- WGCNA::pickSoftThreshold(dataExpr, dataIsExpr = TRUE, powerVector = powers,	corFnc = cor, corOptions = list(use = 'p'),	networkType = "signed")
  			    
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
  			                              deepSplit = 2L,
  			                              pamRespectsDendro = F,
  			                              # detectCutHeight = 0.75,
  			                              minModuleSize = min(20, ncol(dataExpr/2)), # al# 30, #input$minModuleSize, #30,
  			                              # set block size to be number of genes, so that all
  			                              # genes will be analyzed in a single block
  			                              maxBlockSize = input$WGCNAtopNum,#4000,
  			                              
  			                              # == Module Adjustments ==
  			                              reassignThreshold = 0,
  			                              mergeCutHeight = input$mergeCutHeight,#,0.25,
  			                              
  			                              # == TOM == Archive the run results in TOM file (saves time)
  			                              saveTOMs = F,
  			                              saveTOMFileBase = "ER",
  			                              
  			                              # == Output Options
  			                              numericLabels = T,
  			                              verbose = 3L)
  			    cor <- temp_cor
  			  }
  			  netwk
			  })
			})
			
			#### generate dendrogram and gene cluster table #####
			# use input$WGCNAReactive() as event handler to ensure observeEvent() depends on it only
			# and does not directly depends on input$, which ensure WGCNAReactive() will be calculated first.
			observeEvent(WGCNAReactive(),{

			  wgcna <- WGCNAReactive()
			  mergedColors = labels2colors(wgcna$colors)
			  
			    output$Dendrogram <- renderPlot({
	
			      plotDendroAndColors(
			        wgcna$dendrograms[[1]],
			        mergedColors[wgcna$blockGenes[[1]]],
			        "Module colors",
			        dendroLabels = FALSE,
			        hang = 0.03,
			        addGuide = TRUE,
			        guideHang = 0.05 )
			      
			    })
				
				# generate table showing clustered genes #
				ProteinGeneName  <- DataInSets[[working_project()]]$ProteinGeneName
				gene_label <- input$WGCNAgenelable
				
				# t0: merge WGCNA output with ProteinGeneName so that genes can be shown as UniqueID or Gene.Name
				t0 <- tibble::tibble(UniqueID = names(wgcna$colors), color = labels2colors(wgcna$colors)) %>%
				  dplyr::left_join(ProteinGeneName[, c("UniqueID","Gene.Name")], by = "UniqueID") %>%
				  dplyr::select(color,all_of(gene_label)) %>%
				  dplyr::rename(gene = gene_label)
				
				# t1: collapse all genes in a cluster into a cell
				t1 <- t0 %>%
				  dplyr::group_by(color) %>% 
				  dplyr::summarize(n_gene = n(),
				                   gene_group = paste0(gene, collapse = ",")) %>%
				  dplyr::ungroup()
				
				# t2: add the copy button
				t2 <- t1
				t2$copy <- vapply(1L:nrow(t1), function(i){
				  as.character(
				    rclipButton(
				      paste0("clipbtn_", i), 
				      label = "Copy all genes in cluster", 
				      clipText = t1[i, "gene_group"], 
				      #icon = icon("clipboard"),
				      icon = icon("copy", lib = "glyphicon"),
				      class = "btn-primary btn-sm"
				    )
				  )
				}, character(1L))
				
				# rearrange columns
				t2 <- t2 %>% dplyr::select(color, n_gene, copy, gene_group)
				
				output$gene_cluster <- DT::renderDT({
				  DT::datatable(
				    t2,
				    escape = FALSE,
				    selection = "none",
				    colnames=c("Color of cluster", "Number of genes", "Action","Genes in cluster")
				    )
				})
				
			})
			
		}
	)
}
