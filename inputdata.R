###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: inputdata.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 05/23/2023
##@version 3.0
###########################################################################################################
#pkgs: "stringr", "dplyr", "DBI", "pool", "RMariaDB", "RMySQL",  "tidyr", "tibble", "DT"

working_project <- reactiveVal()
DataInSets <- reactiveValues()
DS_names<- reactiveVal() #track loaded projects (not NULL) in DataInSets
saved_plots <- reactiveValues()
saved_table <- reactiveValues()
##################

observe({
	req(input$select_dataset %in% c('Saved Projects in Database', 'Saved Projects in CSV file', 'Public Data(DiseaseLand)'))
	projects = NULL
	if (input$select_dataset=='Saved Projects in Database') {
		library(DBI)
		library(pool)
		pool <- dbPool(RMariaDB::MariaDB(),	dbname='rshiny', host='10.9.68.75', port=3306, user='bgao',	password='sqladmin')
		statement = "SELECT `Project.Title` FROM `projects` ORDER BY `Project.ID` DESC"
		query <- sqlInterpolate(pool, statement)
		projects <- dbGetQuery(pool, query) %>% dplyr::pull(Project.Title)
		con <- poolCheckout(pool)
		poolReturn(con)
		poolClose(pool)
		updateSelectInput(session, "sel_project", choices=c("", projects), selected="")
	}

	if (input$select_dataset=='Saved Projects in CSV file') {
		saved_projects = read.csv("data/saved_projectsNEW.csv")
		projects = saved_projects$ProjectID
		updateSelectInput(session, "sel_project", choices=c("", projects), selected="")
	}

	if (input$select_dataset=='Public Data(DiseaseLand)') {
		library(DBI)
		library(pool)
		pool <- dbPool(RMariaDB::MariaDB(),	dbname='diseaseland',	host='10.9.68.75',	port=3306,	user='bgao', password='sqladmin')
		statement = "SELECT ProjectID from projects"
		query <- sqlInterpolate(pool, statement)
		projects <- dbGetQuery(pool, query) %>% dplyr::pull(ProjectID)
		con <- poolCheckout(pool)
		poolReturn(con)
		poolClose(pool)

		updateSelectizeInput(session, "sel_project", choices = projects, server = TRUE)
	}
})

SavedProjectReactive <- reactive({

	if (input$select_dataset=='Saved Projects in Database') {
		library(DBI)
		library(pool)
		pool <- dbPool(RMariaDB::MariaDB(),	dbname='rshiny', host='10.9.68.75',	port=3306,	user='bgao', password='sqladmin', idleTimeout = 600000)

		saved_projects = dbReadTable(pool, "projects")
		saved_projects <- saved_projects %>%
		dplyr::select(Project.Title, Project.Description, Species) %>%
		dplyr::mutate(ShortNames = Project.Title) %>%
		dplyr::rename(names = Project.Description) %>%
		dplyr::rename(ProjectID = Project.Title)
		con <- poolCheckout(pool)
		poolReturn(con)
		poolClose(pool)
	} else if (input$select_dataset=='Saved Projects in CSV file') {
		saved_projects = read.csv("data/saved_projectsNEW.csv")
	}

	clientData <- session$clientData
	query <- parseQueryString(clientData$url_search)
	if (!is.null(query[['project']])) {
		library(DBI)
		library(pool)
		pool <- dbPool(RMariaDB::MariaDB(),	dbname='rshiny', host='10.9.68.75',	port=3306,	user='bgao', password='sqladmin', idleTimeout = 600000)

		saved_projects = dbReadTable(pool, "projects")
		saved_projects_db <- saved_projects %>%
		dplyr::select(Project.Title, Project.Description, Species) %>%
		dplyr::mutate(ShortNames = Project.Title) %>%
		dplyr::rename(Names = Project.Description) %>%
		dplyr::rename(ProjectID = Project.Title)

		saved_projects_file = read.csv("data/saved_projectsNEW.csv") %>%
		dplyr::select(ProjectID, Names, Species, ShortNames)

		saved_projects <- dplyr::bind_rows(saved_projects_db, saved_projects_file) %>%
		dplyr::distinct()
		con <- poolCheckout(pool)
		poolReturn(con)
		poolClose(pool)
	}
	return(saved_projects)
})

returnlist <- vector(mode='list', length=28)
names(returnlist) <- c("ProjectID", "Name", "Species", "ShortName", "Path", "file1", "file2",
	"MetaData", "MetaData_long","groups", "group_order", "samples", "sample_order", "tests", "tests_order",
	"results_long", "data_long", "ProteinGeneName", "ProteinGeneNameHeader","data_wide", "data_results",
	"comp_info", "sel_comp", "exp_unit",
"statresult", "results_drc","results_omics","results_lin")

#####

DataReactiveRData <- reactive({
	withProgress(message = 'Fetching data.',  detail = 'This may take a while...', value = 0, {
		#browser()
		query <- parseQueryString(isolate(session$clientData$url_search))
		req((!is.null(query[['project']]) & !(query[["project"]] %in% names(DataInSets)))|| input$sel_project!="" || (input$select_dataset=='Upload RData File' & !is.null(input$file1))) #by bgao 0212204

		ProjectID=NULL; ProjectName=NULL;  ShortName=NULL; file1=NULL; file2=NULL; ProjectPath=NULL
		Species = "human"
		exp_unit = "Expression Level"

		if (!is.null(query[['project']])) {
			if (!(query[["project"]] %in% names(DataInSets))) {
				saved_projects <- SavedProjectReactive()
				ProjectID = query[['project']]
				validate(need(ProjectID %in% saved_projects$ProjectID , message = "Please pass a valid ProjectID from URL."))
				ProjectName=saved_projects$Name[saved_projects$ProjectID==ProjectID]
				Species=saved_projects$Species[saved_projects$ProjectID==ProjectID]
				ShortName=saved_projects$ShortNames[saved_projects$ProjectID==ProjectID]
				file1= paste("data/",  ProjectID, ".RData", sep = "")  #data file
				file2= paste("networkdata/", ProjectID, ".RData", sep = "") #Correlation results
				ProjectPath="data/"
				if (is.null(file1) || !file.exists(file1)){
					shinyalert("Oops!", "File does NOT exit.", showConfirmButton = FALSE, showCancelButton = TRUE, type = "error") #by bgao 0212204
				}
				validate(need(file.exists(file1), message = "File does NOT exit.")) #by bgao 0212204
				load(file1)
			}
		}

		if (input$select_dataset %in% c('Saved Projects in Database', 'Saved Projects in CSV file') & input$sel_project!="") {
			saved_projects <- SavedProjectReactive()
			ProjectID=input$sel_project
			ProjectName=saved_projects$Name[saved_projects$ProjectID==ProjectID]
			Species=saved_projects$Species[saved_projects$ProjectID==ProjectID]
			ShortName=saved_projects$ShortNames[saved_projects$ProjectID==ProjectID]
			file1= paste("data/",  ProjectID, ".RData", sep = "")  #data file
			file2= paste("networkdata/", ProjectID, ".RData", sep = "") #Correlation results
			ProjectPath="data/"
			if (is.null(file1) || !file.exists(file1)){
				shinyalert("Oops!", "File does NOT exit.", showConfirmButton = FALSE, showCancelButton = TRUE, type = "error")
			}
			validate(need(file.exists(file1), message = "File does NOT exit."))
			load(file1)
		}

		if (!is.null(query[['serverfile']])) {
			ProjectID = query[['serverfile']]
			if (!is.null(server_dir)) {
				validate(need(file.exists(stringr::str_c(server_dir, "/",  ProjectID, ".csv")),
				message = "Please pass a valid ProjectID from URL. Files must be located in server file folder" ))
				unlisted_project=read.csv(stringr::str_c(server_dir, "/",  ProjectID, ".csv"))
				ProjectName=unlisted_project$Name
				Species=unlisted_project$Species
				ShortName=unlisted_project$ShortName
				file1= paste(server_dir, "/",   ProjectID, ".RData", sep = "")  #data file
				file2= paste(server_dir, "/",  ProjectID, "_network.RData", sep = "") #Correlation results
				ProjectPath=server_dir
				if ("ExpressionUnit" %in% names(unlisted_project)) {
					exp_unit= unlisted_project$ExpressionUnit[1]
				}
			}
			if (is.null(file1) || !file.exists(file1)){
				shinyalert("Oops!", "File does NOT exit.", showConfirmButton = FALSE, showCancelButton = TRUE, type = "error")
			}
			validate(need(file.exists(file1), message = "File does NOT exit."))
			load(file1)

		}

		if (!is.null(query[['unlisted']])) {
			ProjectID = query[['unlisted']]
			validate(need(file.exists(str_c("unlisted/",  ProjectID, ".csv")),
			message = "Please pass a valid ProjectID from URL. Files must be located in unlisted folder" ))
			unlisted_project=read.csv(str_c("unlisted/", ProjectID, ".csv"))
			ProjectID=ProjectID
			ProjectName=unlisted_project$Name
			Species=unlisted_project$Species
			ShortName=unlisted_project$ShortName
			file1 = paste("unlisted/",  ProjectID, ".RData", sep = "")  #data file
			file2 = paste("unlisted/", ProjectID, "_network.RData", sep = "") #Correlation results
			ProjectPath="unlisted/"
			if ("ExpressionUnit" %in% names(unlisted_project)) {
				exp_unit= unlisted_project$ExpressionUnit[1]
			}
			if (is.null(file1) || !file.exists(file1)){
				shinyalert("Oops!", "File does NOT exit.", showConfirmButton = FALSE, showCancelButton = TRUE, type = "error")
			}
			validate(need(file.exists(file1), message = "File does NOT exit."))
			load(file1)
		}

		if (input$select_dataset=='Upload RData File' & !is.null(input$file1)) {
			ProjectID = str_replace(input$file1$name, regex(".RData", ignore_case = TRUE), "")
			ProjectName = ShortName =  input$project_name
			Species = input$species
			file1 = input$file1$datapath
			file2 = input$file2$datapath
			ProjectPath = NULL
			if (is.null(file1) || !file.exists(file1)){
				shinyalert("Oops!", "File does NOT exit.", showConfirmButton = FALSE, showCancelButton = TRUE, type = "error")
			}
			validate(need(file.exists(file1), message = "File does NOT exit."))
			load(file1)
		}

		returnlist[["ProjectID"]] = ProjectID
		returnlist[["Name"]] = ProjectName
		returnlist[["Species"]] = Species
		returnlist[["ShortName"]] = ShortName
		returnlist[["Path"]] = ProjectPath
		returnlist[["file1"]] = file1
		returnlist[["file2"]] = file2
		returnlist[["exp_unit"]] <- exp_unit

		#rename list
		lookup <- c("UniqueID" = "uniqueID", "group" = "Group")

		if (exists("data_wide")) {
			if (!is.data.frame(data_wide)) {
				data_wide <- data.frame(data_wide, check.names = FALSE)
			}  #change data_wide to data frame from numeric matrix if needed
			returnlist[["data_wide"]] = data_wide
		}

		if (exists("ProteinGeneName")) {
			if (!"Protein.ID" %in% names(ProteinGeneName)) {
				ProteinGeneName$Protein.ID=NA
			} #Add Protein.ID column as it is required for certain tools.
			ProteinGeneNameHeader = colnames(ProteinGeneName)
			ProteinGeneNameHeader <- intersect(ProteinGeneNameHeader, c("UniqueID","Gene.Name","Protein.ID"))

			# modify by bgao 08/02/2023, keep ProteinGeneName ids only in result table
			ProteinGeneName <- ProteinGeneName %>%
			dplyr::rename(any_of(lookup)) %>%
			dplyr::mutate_if(is.factor, as.character)  %>%
			dplyr::filter(UniqueID %in% (results_long %>% dplyr::pull(UniqueID) %>% unique()))

			returnlist[["ProteinGeneName"]] = ProteinGeneName
			returnlist[["ProteinGeneNameHeader"]] = ProteinGeneNameHeader
		}

		if (exists("data_results")) {
			returnlist[["data_results"]] = data_results %>% dplyr::rename(any_of(lookup))
		}

		if (exists("results_long")) {
			results_long <- results_long %>%
			dplyr::rename(any_of(lookup)) %>%
			dplyr::mutate_if(is.factor, as.character)  %>% dplyr::inner_join(ProteinGeneName, ., by = "UniqueID")
			tests <- tests_order <- unique(as.character(results_long$test))
			returnlist[["tests"]] = tests
			returnlist[["tests_order"]] = tests_order
			returnlist[["results_long"]] = results_long
		}

		if (exists("MetaData")) {
			MetaData <- MetaData %>%
			dplyr::rename(any_of(lookup))

			if ("Order" %in% names(MetaData))
			groups <- group_order <- as.character(MetaData$Order[MetaData$Order != ""])
			else
			groups <- group_order <- MetaData %>% dplyr::pull(group) %>% unique()

			samples <- sample_order <- as.character(MetaData$sampleid[order(match(MetaData$group,groups))])

			### meta data to long form
			MetaData_long <- MetaData %>%
			dplyr::select(-any_of(c("Order", "ComparePairs"))) %>%
			#dplyr::mutate_if(is.numeric, as.character) %>% #this will fail when there are columns in Time format
			dplyr::mutate_all(as.character) %>%
			tidyr::pivot_longer(cols = -sampleid,  names_to = "type",values_to = "group")

			returnlist[["MetaData"]] = MetaData
			returnlist[["MetaData_long"]] = MetaData_long
			returnlist[["groups"]] = groups
			returnlist[["group_order"]] = group_order
			returnlist[["samples"]] = samples
			returnlist[["sample_order"]] = sample_order
		}

		if (exists("data_long")) {
			data_long <- data_long %>%
			dplyr::rename(any_of(lookup))

			if (exists("ProteinGeneName")) {
				data_long <- data_long %>% dplyr::mutate_if(is.factor, as.character)  %>% dplyr::inner_join(ProteinGeneName, ., by = "UniqueID")
				if (length(groups) == 0) {
					groups <- group_order <- as.character(unique(data_long$group))
				}
			} else {
				groups <- group_order <- as.character(unique(data_long$group))
			}
			returnlist[["data_long"]] = data_long
			returnlist[["groups"]] = groups
			returnlist[["group_order"]] = group_order
		}

		if (exists("comp_info")){
			sel_comp<-data.frame(Comparison=rownames(comp_info), comp_info)%>%dplyr::filter(Group_name!="", !is.na(Group_name)); dim(sel_comp)
			#sel_comp<-data.frame(Comparison=rownames(comp_info), comp_info)%>%dplyr::filter(str_detect(Subsetting_group, ":")); dim(sel_comp)
			if (nrow(sel_comp)>0) {
				sel_comp<-sel_comp%>%dplyr::mutate(N_samples=0, sample_list=NA, subset_list=NA)
				for (i in 1:nrow(sel_comp)) {
					sel_samples<-rep(TRUE, nrow(MetaData))
					if  (str_detect(sel_comp$Subsetting_group[i], ":")){
						sg1<-str_split(sel_comp$Subsetting_group[i], ";")[[1]]
						for (j in 1:length(sg1)){
							sub_values=str_split(sg1[j], ":")[[1]]
							sel_j=MetaData[[sub_values[1]]]==sub_values[2]
							sel_samples=sel_samples & sel_j
						}
						#sel_comp$N_samples[i]=sum(sel_samples)
						sel_comp$subset_list[i]=paste(MetaData$sampleid[sel_samples], collapse = ",")
					}
					#now further filter for Group_test and Group_ctrl
					sel_DEG_samples<- (MetaData[[sel_comp$Group_name[i]]] %in% c(sel_comp$Group_test[i], sel_comp$Group_ctrl[i]) )
					sel_samples = sel_samples & sel_DEG_samples
					sel_comp$N_samples[i]=sum(sel_samples)
					sel_comp$sample_list[i]=paste(MetaData$sampleid[sel_samples], collapse = ",")
				}
				# cat(i, sg1, sum(sel_samples), paste(MetaData$sampleid[sel_samples], collapse = ","), "\n\n")
			}
			sel_comp<-sel_comp%>%dplyr::filter(N_samples>0)
			if (nrow(sel_comp)==0) {sel_comp=NULL}

			returnlist[["sel_comp"]]=sel_comp
		}

		if (exists("results_stat")) {
			returnlist[["statresult"]] = results_stat %>% dplyr::rename(any_of(lookup))
			pcutoffdf  <- rbind(data.frame("Type" = "Pvalue", "LessThan0.05" = sum(results_stat$pvalue < 0.05,na.rm = TRUE), "LessThan0.01" = sum(results_stat$pvalue < 0.01,na.rm = TRUE), "LessThan0.001" = sum(results_stat$pvalue < 0.001,na.rm = TRUE)),
				data.frame("Type" = "Padjust", "LessThan0.05" = sum(results_stat$padjust < 0.05,na.rm = TRUE), "LessThan0.01" = sum(results_stat$padjust < 0.01,na.rm = TRUE), "LessThan0.001" = sum(results_stat$padjust < 0.001,na.rm = TRUE))
			)
			returnlist[["pcutoffdf"]] = pcutoffdf
		}

		if (exists("results_drc")) {

			returnlist[["results_drc"]] = results_drc %>% dplyr::rename(any_of(lookup))
		}

		if (exists("results_omics")) {
			results_omics <- results_omics %>% dplyr::filter(model!="NA") %>% dplyr::rename(any_of(lookup))
			returnlist[["results_omics"]] = results_omics
		}


		if (exists("results_lin")) {
			returnlist[["results_lin"]] = results_lin %>% dplyr::rename(any_of(lookup))
		}
		return(returnlist)
	})
})

DataReactiveDB <- reactive({
	withProgress(message = 'Fetching data.',  detail = 'This may take a while...', value = 0, {

		#validate(
		#	need(input$table_cell_clicked!="", message = "Please select a project")
		#)

		#info = input$table_cell_clicked
		#ProjectID = info$value

		#ProjectID = "GSE51799"
		#ProjectID = "GSE51684"
		#ProjectID = "GSE11227"
		#ProjectID = "GSE1145"
		#ProjectID = "GSE18956"

		ProjectID <- input$sel_project
		library(DBI)
		library(pool)
		pool <- dbPool(RMariaDB::MariaDB(),	dbname='diseaseland',	host='10.9.68.75',	port=3306,	user='bgao', password='sqladmin', idleTimeout = 600000)

		#geneannotation
		statement = "SELECT GeneIndex, GeneID, GeneName FROM `geneannotation`"
		query <- sqlInterpolate(pool, statement)
		geneannotation <- dbGetQuery(pool, query) %>%
		dplyr::mutate(UniqueID = GeneID)

		ProteinGeneName <- geneannotation %>% dplyr::select(one_of(c("UniqueID","GeneName","GeneID")))  %>%
		`colnames<-`(c("UniqueID", "Gene.Name", "Protein.ID"))

		geneannotation <- geneannotation %>% dplyr::select(one_of(c("GeneIndex","UniqueID")))

		#project info
		statement = paste("SELECT * from projects WHERE ProjectID = '", ProjectID, "'", sep="")
		query <- sqlInterpolate(pool, statement)
		projectinfo <- dbGetQuery(pool, query)

		#samples info
		statement = paste("SELECT * FROM `samples` WHERE 	ProjectName = '", ProjectID, "'", sep="")
		query <- sqlInterpolate(pool, statement)
		samples <- dbGetQuery(pool, query)

		#comparison info for the project
		statement = paste("SELECT ComparisonIndex, variable, value FROM `comparisonslong` WHERE ComparisonIndex IN (SELECT ComparisonIndex FROM `comparisonslong` WHERE variable = 'ProjectName' AND value = '",ProjectID, "')", sep="")
		query <- sqlInterpolate(pool, statement)
		comparisons <- dbGetQuery(pool, query)
		ComparisonContrast <- dplyr::filter(comparisons, variable == "ComparisonContrast" ) %>%
		dplyr::select(value) %>%
		tidyr::separate(., value, c("type","contrast"),sep = " => ")

		types = c()
		for (i in unique(ComparisonContrast$type)){
			types = c(types, unlist(strsplit(i, "[:]")))
		}
		types <- unique(types)

		variable.need <- c(paste("Case", types,sep="."), "Case.SampleIDs")
		Case <- comparisons %>%
		dplyr::filter(grepl("^Case", variable)) %>%
		dplyr::filter(variable %in% variable.need) %>%
		spread(variable,value)
		colnames(Case) <- gsub("Case.","", colnames(Case))

		variable.need <- c(paste("Control", types,sep="."), "Control.SampleIDs")
		Control <- comparisons %>%
		dplyr::filter(grepl("^Control", variable)) %>%
		dplyr::filter(variable %in% variable.need) %>%
		spread(variable,value)
		colnames(Control) <- gsub("Control.","", colnames(Control))

		Case.Control <- bind_rows(Case,Control) %>%
		dplyr::select(-"ComparisonIndex")

		groups = data.frame()
		for (i in 1:nrow(Case.Control)) {
			SampleID <- Case.Control[i,'SampleIDs']
			SampleID <- unlist(strsplit(SampleID, ";"))
			SampleID.df <- as.data.frame(SampleID) %>%
			mutate_if(is.factor, as.character)
			types <- dplyr::select(Case.Control[i,], -one_of('SampleIDs'))
			#groups <- bind_rows(groups, bind_cols(SampleID.df,types[rep(seq_len(nrow(types)), each= nrow(SampleID.df)),]))
			groups <- bind_rows(groups, cbind(SampleID.df,types,row.names = NULL))
		}

		groups.1 <- data.frame()
		for (i in unique(groups$SampleID)){
			oneID <- filter(groups, SampleID==i)
			groups.1 <- bind_rows(groups.1, unlist(apply(oneID, 2, function(x)unique(x[!is.na(x)]))))
		}

		if (ncol(groups.1) ==2) {
			MetaData <- groups.1 %>%
			`colnames<-`(c("sampleid", "group"))
		} else {
			MetaData <- groups.1 %>%
			dplyr::rename(sampleid = SampleID) %>%
			tidyr::unite("group", -sampleid, sep=".", remove = FALSE)
		}

		### meta data to long form
		MetaData_long <- MetaData %>%
		tidyr::pivot_longer(cols = -sampleid,  names_to = "type",values_to = "group")

		#results_long
		comparisonslist <- comparisons %>% dplyr::pull(ComparisonIndex) %>% unique() %>% paste(.,collapse =",")
		statement = paste("SELECT ComparisonIndex, GeneIndex, Log2FoldChange, PValue, AdjustedPValue
		FROM `comparisondata`	WHERE PValue IS NOT NULL AND ComparisonIndex IN (", comparisonslist , ")", sep="")
		query <- sqlInterpolate(pool, statement)
		results_long <- dbGetQuery(pool, query)

		results_long <- results_long %>%
		dplyr::filter(!is.na(PValue))  %>%
		#dplyr::mutate_if(is.factor, as.character)  %>%
		dplyr::group_by(ComparisonIndex, GeneIndex) %>%
		dplyr::slice_min(PValue, with_ties = FALSE) %>%
		dplyr::right_join(geneannotation,., by = "GeneIndex")	%>%
		dplyr::left_join(.,(dplyr::filter(comparisons, variable == "ComparisonContrast") %>% dplyr::select(-variable)),  by = "ComparisonIndex")	%>%
		tidyr::separate(., value, c("type","contrast"),sep = " => ") %>%
		dplyr::select(one_of(c("UniqueID","contrast", "AdjustedPValue", "PValue", "Log2FoldChange"))) %>%
		`colnames<-`(c("UniqueID", "test", "Adj.P.Value","P.Value","logFC"))

		#samples in the project
		SampleIdxlist <- samples %>% dplyr::pull(SampleIndex) %>% unique() %>% paste(.,collapse =",")
		#data_long
		statement = paste("SELECT SampleIndex, GeneIndex, Count as expr
		FROM genefpkm WHERE Count != 0 AND SampleIndex IN (",SampleIdxlist,")", sep="")
		query <- sqlInterpolate(pool, statement)
		data_long <- dbGetQuery(pool, query)

		if (nrow(data_long) == 0) {
			statement = paste("SELECT SampleIndex, GeneIndex, GeneExpression as expr
			FROM genelevelexpression WHERE SampleIndex IN (",SampleIdxlist,")", sep="")
			query <- sqlInterpolate(pool, statement)
			data_long <- dbGetQuery(pool, query)
		}

		print(format(object.size(data_long), unit = 'auto'))

		data_long <- data_long %>%
		dplyr::mutate_if(is.numeric, round, 0) %>%
		dplyr::filter(expr != 0) %>%
		dplyr::mutate_if(is.factor, as.character)  %>%
		dplyr::inner_join(geneannotation,., by = "GeneIndex") %>%
		dplyr::select(-GeneIndex) %>%
		dplyr::filter(UniqueID %in% unique(results_long$UniqueID)) %>%
		dplyr::left_join(., (samples %>% dplyr::select(SampleID, SampleIndex)), by = "SampleIndex") %>%
		dplyr::rename(sampleid = SampleID) %>%
		dplyr::select(-SampleIndex) %>%
		dplyr::left_join(., MetaData, by="sampleid")

		#meanExpr <- with(data_long, ave(expr, UniqueID, FUN = mean))
		#data_long <- data_long[meanExpr >= 1, ]

		print(format(object.size(data_long), unit = 'auto'))

		data_wide <- data_long %>%
		dplyr::select(one_of(c("UniqueID","sampleid","expr"))) %>%
		spread(., sampleid, expr) %>%
		tibble::column_to_rownames(.,"UniqueID")
		print(format(object.size(data_wide), unit = 'auto'))

		group_names <- MetaData %>%
		dplyr::pull('group') %>% unique()

		sample_names <- MetaData %>%
		dplyr::pull(sampleid) %>% unique()

		tests  <- unique(as.character(ComparisonContrast$contrast))

		results_long <- results_long %>% mutate_if(is.factor, as.character)  %>%
		inner_join(ProteinGeneName,., by = "UniqueID") %>%
		dplyr::filter(UniqueID %in% rownames(data_wide))

		#data_long <- data_long %>% mutate_if(is.factor, as.character)  %>%
		#inner_join(ProteinGeneName,., by = "UniqueID") #%>%
		#dplyr::filter(UniqueID %in% rownames(data_wide)) %>%
		#dplyr::rename(expr = Count)

		poolClose(pool)

		ProteinGeneName <- ProteinGeneName %>%
		dplyr::inner_join(data.frame(UniqueID=rownames(data_wide), Intensity=apply(data_wide,1,mean, na.rm=TRUE)), by = join_by(UniqueID)) %>%
		dplyr::mutate_if(is.numeric, round, 0)

		##data_results
		#data_results <- ProteinGeneName %>%
		#dplyr::select(any_of(c("UniqueID","Gene.Name","Protein.ID"))) %>%
		#dplyr::left_join(data.frame(UniqueID=rownames(data_wide), Intensity=apply(data_wide,1,mean, na.rm = TRUE)) %>%	dplyr::filter(!duplicated(UniqueID)),by = join_by(UniqueID))
		#
		#sinfo1 <- data.frame(sampleid=names(data_wide)) %>%
		#dplyr::left_join(MetaData%>%dplyr::select(sampleid, group), by = join_by(sampleid))
		#
		#for(grp in unique(sinfo1$group)){
		#	subdata <- data.frame(UniqueID=rownames(data_wide), t(apply(data_wide[,sinfo1$group==grp, drop=FALSE],1,function(x)return(setNames(c(mean(x),sd(x)),paste(grp,c("Mean","sd"),sep="_"))))), check.names=FALSE )
		#	data_results <- data_results %>%
		#	dplyr::left_join(subdata %>% dplyr::filter(!duplicated(UniqueID)), by = join_by(UniqueID))
		#}
		#
		#for (ctr in tests) {
		#	subdata <- results_long %>%
		#	dplyr::filter(test==ctr) %>%
		#	dplyr::select(UniqueID, logFC, P.Value, Adj.P.Value)
		#	names(subdata)[2:4]=str_c(ctr, "_", names(subdata)[2:4])
		#	data_results <- data_results %>%
		#	dplyr::left_join(subdata %>% dplyr::filter(!duplicated(UniqueID)),by = join_by(UniqueID))
		#}
		#

		returnlist[["ProjectID"]] <- ProjectID
		returnlist[["Name"]] <-  ProjectID
		returnlist[["Species"]] <- unique(samples$Organism)
		returnlist[["ShortName"]] <-  ProjectID
		returnlist[["Path"]] <- NULL
		returnlist[["file1"]] <- NULL
		returnlist[["file2"]] <- NULL
		returnlist[["exp_unit"]] <- "expression level"
		returnlist[["MetaData"]] = MetaData
		returnlist[["MetaData_long"]] = MetaData_long
		returnlist[["ProteinGeneName"]] = ProteinGeneName
		returnlist[["ProteinGeneNameHeader"]] = colnames(ProteinGeneName)
		returnlist[["data_long"]] = data_long
		returnlist[["data_wide"]] = data_wide
		returnlist[["results_long"]] = results_long
		returnlist[["data_results"]] = results_long
		returnlist[["groups"]] = group_names
		returnlist[["group_order"]] = group_names
		returnlist[["samples"]] = sample_names
		returnlist[["sample_order"]] = sample_names
		returnlist[["tests"]] = tests
		returnlist[["tests_order"]] = tests

		return(returnlist)
	})
})

DataReactiveTxt <- reactive({
	withProgress(message = 'Fetching data.',  detail = 'This may take a while...', value = 0, {

		req(input$F_project_name!="" & ((!is.null(input$F_sample) & !is.null(input$F_exp)) | (!is.null(input$F_comp))))

		cleanup_empty <- function(df) {
			df.empty=(is.na(df) | df=="")
			selCol=!(colSums( df.empty)==nrow(df))
			selRow=!(rowSums( df.empty)==ncol(df))
			return(df[selRow, selCol])
		}

		#create unique project ID
		Project_name=input$F_project_name
		ProjectID=str_c("PRJ_",  make.names(Project_name) )
		if (length(ProjectID>45)) {
			ProjectID = substr(ProjectID, 1, 45)
		}
		ProjectID=str_c(ProjectID,"_", stringi::stri_rand_strings(1,6) )
		#species
		species <- input$Fspecies

		#MetaData
		if (!is.null(input$F_sample)) {
			MetaData <- read.csv(input$F_sample$datapath, header=T, check.names=F)
			MetaData <- cleanup_empty(MetaData)

			if ("Order" %in% names(MetaData))
			groups <- group_order <- as.character(MetaData$Order[MetaData$Order != ""])
			else
			groups <- group_order <- MetaData %>% dplyr::pull(group) %>% unique()

			samples <- sample_order <- as.character(MetaData$sampleid[order(match(MetaData$group,groups))])

			### meta data to long form
			MetaData_long <- MetaData %>%
			dplyr::select(-any_of(c("Order", "ComparePairs"))) %>%
			tidyr::pivot_longer(cols = -sampleid,  names_to = "type",values_to = "group")
		}

		#data_wide
		if (!is.null(input$F_exp)) {

			exp_file=input$F_exp$datapath
			if (str_detect(exp_file, "gz$") ) {
				exp_file=gzfile(exp_file, "rt")
			}

			if (str_detect(exp_file, "zip$") ) {
				fnames = as.character(unzip(exp_file, list = TRUE)$Name)
				exp_file=unz(exp_file, fnames[1])
			}

			data_wide=read.csv(exp_file, row.name=1, header=T, check.names=F)
			data_wide=cleanup_empty(data_wide)

			#data_long
			data_long <- data.table::melt(as.matrix(data_wide))
			colnames(data_long) <- c("UniqueID","sampleid","expr")
			data_long <- data_long %>% dplyr::mutate(sampleid=as.character(sampleid))
			data_long <- data_long %>% dplyr::left_join(MetaData %>% dplyr::select(sampleid, group), by = join_by(sampleid))
		}

		#results_long
		if (!is.null(input$F_comp)) {
			comp_file=input$F_comp$datapath
			if (str_detect(comp_file, "gz$")){
				comp_file=gzfile(comp_file, "rt")
			}

			if (str_detect(comp_file, "zip$")){
				fnames = as.character(unzip(comp_file, list = TRUE)$Name)
				comp_file=unz(comp_file, fnames[1])
			}

			results_long=read.csv(comp_file, header=T, check.names=F)
			results_long=cleanup_empty(results_long)


			IDs=rownames(data_wide);
			IDs2=results_long$UniqueID
			IDall=unique(c(IDs, IDs2));
			IDall=IDall[!is.na(IDall)]
			IDall=IDall[!(IDall=="")]
			tests = sort(unique(results_long$test))
		}

		#ProteinGeneName
		if (input$F_annot_auto==0) {
			ProteinGeneName = read.csv(input$F_annot$datapath)
			ProteinGeneName = cleanup_empty(ProteinGeneName)
			ProteinGeneName <- ProteinGeneName %>%
			dplyr::filter(UniqueID %in% IDall)

			if (!"Protein.ID" %in% names(ProteinGeneName)) {
				ProteinGeneName$Protein.ID=NA
			} #Add Protein.ID column as it is required for certain tools.

		} else {
			if (str_detect(input$F_ID_type, "UniProt") ) { #protein name match
				ProteinInfo <- readRDS('db/ProteinInfo.rds') %>%
				dplyr::mutate(Gene_Name=str_replace(Gene_Name, " .+", "")) #replace space, as UniProt put alias here

				if (input$F_ID_type=="UniProtKB Protein ID" ) {
					ProteinGeneName <- data.frame(id=1:length(IDall), UniqueID=IDall) %>%
					dplyr::left_join(ProteinInfo %>%
						dplyr::transmute(UniqueID=UniProtKB.AC, Gene.Name=Gene_Name, Protein.ID=UniProtKB.AC, Description=Protein_Name) %>%
						dplyr::filter(!duplicated(UniqueID)),by = join_by(UniqueID)
					)
				} else {
					ProteinGeneName <- data.frame(id=1:length(IDall), UniqueID=IDall) %>%
					dplyr::left_join(ProteinInfo%>%
						dplyr::transmute(UniqueID=UniProtKB.ID, Gene.Name=Gene_Name, Protein.ID=UniProtKB.AC, Description=Protein_Name) %>%
						dplyr::filter(!duplicated(UniqueID)), by = join_by(UniqueID)
					)
				}
				if (input$F_description==0) {
					ProteinGeneName <- ProteinGeneName %>%
					dplyr::select(-Description)
				}

				if (input$F_fillName==1) {
					ProteinGeneName <- ProteinGeneName %>%
					dplyr::mutate(Gene.Name=ifelse(is.na(Gene.Name), UniqueID, Gene.Name))
				}

			} else { #gene
				#cat("working on ",species," genes for project", ProjectID, "\n")
				if (species=="rat") {
					ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl")
				} else if (species=="mouse") {
					ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
				} else {
					ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
				}
				#setProgress(0.2, detail = "Connected to Biomart, converting IDs to gene names..."); Sys.sleep(0.1)

				if (input$F_ID_type=="Ensembl Gene ID" ) {
					filter_type="ensembl_gene_id"
					IDall_old=IDall
					IDall=str_replace(IDall, "\\.\\d+$", "")
					EID=data.frame(IDall_old, IDall)
				} else if (input$F_ID_type=="NCBI GeneID" ) {
					filter_type="entrezgene_id"
				} else if (input$F_ID_type=="Gene Symbol" ) {
					filter_type="external_gene_name"
				}

				E_attributes <- c( 'ensembl_gene_id', "external_gene_name", "gene_biotype","entrezgene_id")
				if (input$F_description==1) {
					E_attributes<-c(E_attributes, "description")
				}

				system.time( output<-getBM(attributes = E_attributes,filters = filter_type, values = IDall, mart = ensembl, useCache=FALSE) )

				if (nrow(output)==0) {
					error_message="No gene annotation extracted from Biomart. Did you select the correct Species and Unique ID Type?"
					upload_message(error_message)
				}
				validate(need(nrow(output)>0, message = "No gene annotation extracted from Biomart. Did you select the correct Species and Unique ID Type?"))

				output <- output %>% dplyr::arrange(entrezgene_id, ensembl_gene_id) #Favor IDs with smaller numbers

				if (input$F_description==0) {
					output$description=NA
				}
				F_TYPE=sym(filter_type)

				#ProteinGeneName <- data.frame(id=1:length(IDall), UniqueID=IDall) %>%
				ProteinGeneName <- data.frame(UniqueID=IDall) %>%
				dplyr::left_join(output%>%
					dplyr::transmute(UniqueID=!!F_TYPE, Gene.Name=external_gene_name, GeneType=gene_biotype, Description=description) %>%
					dplyr::filter(!duplicated(UniqueID))
				)

				if (input$F_ID_type=="Ensembl Gene ID" ) {
					ProteinGeneName <- EID %>%
					dplyr::transmute(UniqueID=IDall_old, Unique1=IDall) %>%
					left_join(ProteinGeneName%>%mutate(Unique1=UniqueID)%>%dplyr::select(-UniqueID)) %>%
					dplyr::select(-Unique1)
				}

				ProteinGeneName$Protein.ID=NA
				ProteinGeneName <- ProteinGeneName %>%
				#dplyr::select(id, UniqueID, Gene.Name, Protein.ID, GeneType, Description)
				dplyr::select(UniqueID, Gene.Name, Protein.ID, GeneType, Description)

				if (input$F_description==0) {
					ProteinGeneName <- ProteinGeneName %>%
					dplyr::select(-Description)
				}

				if (input$F_fillName==1) {
					ProteinGeneName <- ProteinGeneName %>%
					dplyr::mutate(Gene.Name=ifelse(is.na(Gene.Name), UniqueID, Gene.Name) )
				}
			}
		}


		#data_results
		if(!is.null(ProteinGeneName) & !is.null(data_wide) & !is.null(results_long)) {
			data_results <- ProteinGeneName %>%
			dplyr::select(any_of(c("UniqueID","Gene.Name","Protein.ID"))) %>%
			dplyr::left_join(data.frame(UniqueID=rownames(data_wide), Intensity=apply(data_wide,1,mean, na.rm = TRUE)) %>%	dplyr::filter(!duplicated(UniqueID)),by = join_by(UniqueID))

			sinfo1 <- data.frame(sampleid=names(data_wide)) %>%
			dplyr::left_join(MetaData%>%dplyr::select(sampleid, group), by = join_by(sampleid))

			for(grp in unique(sinfo1$group)){
				subdata <- data.frame(UniqueID=rownames(data_wide), t(apply(data_wide[,sinfo1$group==grp, drop=FALSE],1,function(x)return(setNames(c(mean(x),sd(x)),paste(grp,c("Mean","sd"),sep="_"))))), check.names=FALSE )
				data_results <- data_results %>%
				dplyr::left_join(subdata %>% dplyr::filter(!duplicated(UniqueID)), by = join_by(UniqueID))
			}

			for (ctr in tests) {
				subdata <- results_long %>%
				dplyr::filter(test==ctr) %>%
				dplyr::select(UniqueID, logFC, P.Value, Adj.P.Value)
				names(subdata)[2:4]=str_c(ctr, "_", names(subdata)[2:4])
				data_results <- data_results %>%
				dplyr::left_join(subdata %>% dplyr::filter(!duplicated(UniqueID)),by = join_by(UniqueID))
			}
		}

		file1 = stringr::str_c(input$folder_name, "/", ProjectID, ".RData")
		if (input$savetoserver == "YES"){
			save(data_long,data_results,data_wide,MetaData,ProteinGeneName,results_long,file=file1)
		}


		ProjectPath <- "unlisted/"
		exp_unit <- "Expression Level"
		#samples <- sort(unique(MetaData$sampleid))

		group_order <- groups
		sample_order <- samples
		tests_order <- tests
		ProteinGeneNameHeader = colnames(ProteinGeneName)

		results_long <- results_long %>% dplyr::mutate_if(is.factor, as.character)  %>% dplyr::inner_join(ProteinGeneName, ., by = "UniqueID")

		returnlist[["ProjectID"]] <- ProjectID
		returnlist[["Name"]] <- Project_name
		returnlist[["Species"]] <- species
		returnlist[["ShortName"]] <- Project_name
		returnlist[["Path"]] <- ProjectPath
		returnlist[["file1"]] <- file1
		# Unlike in DataReactiveRdata(), 'file2' is not defined in DataReactiveTxt()
		# Hence, set returnlist[["file2]] = NULL; otherwise error message shows up
		# when uploading csv files: "object file2 not found"
		returnlist[["file2"]] <- NULL
		returnlist[["exp_unit"]] <- exp_unit
		returnlist[["MetaData"]] = MetaData
		returnlist[["MetaData_long"]] = MetaData_long
		returnlist[["ProteinGeneName"]] = ProteinGeneName
		returnlist[["ProteinGeneNameHeader"]] = ProteinGeneNameHeader
		returnlist[["data_long"]] = data_long
		returnlist[["data_wide"]] = data_wide
		returnlist[["results_long"]] = results_long
		returnlist[["data_results"]] = data_results
		returnlist[["groups"]] = groups
		returnlist[["group_order"]] = group_order
		returnlist[["samples"]] = samples
		returnlist[["sample_order"]] = sample_order
		returnlist[["tests"]] = tests
		returnlist[["tests_order"]] = tests_order

		return(returnlist)
	})
})

##
observe({
	if (length(names(DataInSets)) == 0) {
		output$loaddata <- renderUI({
			actionButton("load", "Load Dataset")
		})
	} else {
		output$loaddata <- renderUI({
			actionButton("adddata", "Add Another Dataset")
		})
		output$removedata <- renderUI({
			actionButton("removedata", "Remove Current Working Dataset")
		})
	}
})

#load project in csv or database
observeEvent(input$load | input$adddata | input$uploadData | input$customData, {
	query <- parseQueryString(isolate(session$clientData$url_search))
	req((!is.null(query[['project']]) & !(query[['project']] %in% names(DataInSets)))|| input$sel_project!="" || (input$select_dataset=='Upload RData File' & !is.null(input$file1)))

	if (input$select_dataset == 'Public Data(DiseaseLand)') {
		DataIn <- DataReactiveDB()
	}  else
	if (input$select_dataset == 'Upload Data Files (csv)') {
		DataIn <- DataReactiveTxt()
	}else
	if (!is.null(input$load)) {
		DataIn <- DataReactiveRData()
	}else
	if (!is.null(input$adddata)) {
		DataIn <- DataReactiveRData()
	} else if (!is.null(query[['project']]) & !(query[['project']] %in% names(DataInSets))){
		DataIn <- DataReactiveRData()
	} else {
		return()
	}

	ProjectID <- DataIn$ProjectID
	working_project(ProjectID)
	DataInSets[[ProjectID]]  <-  DataIn
	DataInSets_List<-reactiveValuesToList(DataInSets)
	DataInSets_List<-DataInSets_List[!sapply(DataInSets_List, is.null)]
	DS_names(names(DataInSets_List))

	if (length(names(DataInSets)) == 1) {
		if (!(modulelist[15] %in% saved_setting$value)) {
			saved_setting$value <- c(saved_setting$value, modulelist[15])
			source(moduleFilelist[15],  local = TRUE)
			appendTab(session=session, inputId = "menu", tabPanel(modulelist[15],	groupsample_ui("15")))
			groupsample_server(id = "15")
		}
	}

	if (!is.null(DataInSets[[working_project()]]$data_wide)) {
		if (!(modulelist[1] %in% saved_setting$value)) {
			saved_setting$value <- c(saved_setting$value, modulelist[1])
			source(moduleFilelist[1],  local = TRUE)
			appendTab(session=session, inputId = "menu", tabPanel(modulelist[1],	qcplot_ui("1")))
			qcplot_server(id = "1")
		}
	}

	if (!is.null(DataInSets[[working_project()]]$results_long)) {

		if (!(modulelist[2] %in% saved_setting$value)) {
			saved_setting$value <- c(saved_setting$value, modulelist[2])
			source(moduleFilelist[2],  local = TRUE)
			appendTab(session=session, inputId = "menu", tabPanel(modulelist[2],	deg_ui("deg")))
			deg_server(id = "deg")
		}

		if (!(modulelist[4] %in% saved_setting$value)) {
			saved_setting$value <- c(saved_setting$value, modulelist[4])
			source(moduleFilelist[4],  local = TRUE)
			appendTab(session=session, inputId = "menu", tabPanel(modulelist[4],	expression_ui("4")))
			expression_server(id = "4")
		}

		if (!(modulelist[8] %in% saved_setting$value)) {
			saved_setting$value <- c(saved_setting$value, modulelist[8])
			source(moduleFilelist[8],  local = TRUE)
			appendTab(session=session, inputId = "menu", tabPanel(modulelist[8],	venn_ui("venn")))
			venn_server(id = "venn")
		}
	}

	if (!is.null(DataInSets[[working_project()]]$results_drc)) {
		if (!(modulelist[11] %in% saved_setting$value)) {
			saved_setting$value <- c(saved_setting$value, modulelist[11])
			source(moduleFilelist[11],  local = TRUE)
			appendTab(session=session, inputId = "menu", tabPanel(modulelist[11],	drc_ui("11")))
			drc_server(id = "11")
		}
	}

	if (!is.null(DataInSets[[working_project()]]$results_omics)) {
		if (!(modulelist[12] %in% saved_setting$value)) {
			saved_setting$value <- c(saved_setting$value, modulelist[12])
			source("util-basicandfitfunc.R",local=TRUE)$value
			source(moduleFilelist[12],  local = TRUE)
			appendTab(session=session, inputId = "menu", tabPanel(modulelist[12],	dromics_ui("12")))
			dromics_server(id = "12")
		}
	}

	if (!is.null(DataInSets[[working_project()]]$results_lin)) {
		if (!(modulelist[13] %in% saved_setting$value)) {
			saved_setting$value <- c(saved_setting$value, modulelist[13])
			source(moduleFilelist[13],  local = TRUE)
			appendTab(session=session, inputId = "menu",  tabPanel(modulelist[13],	monotonic_ui("13")))
			monotonic_server(id = "13")
		}
	}
})

observeEvent(input$removedata, {
	DataInSets[[working_project()]] <- NULL
	DataInSets_List<-reactiveValuesToList(DataInSets)
	DataInSets_List<-DataInSets_List[!sapply(DataInSets_List, is.null)]
	if (length(names(DataInSets_List))>0) {
		working_project(names(DataInSets_List)[1])
	} else {
		working_project(NULL)
	}

	DS_names(names(DataInSets_List))
	#make(lock_envir = FALSE)
	#rm(currentproject, envir = .subset2(DataInSets, "impl")$.values)
	#.subset2(DataInSets, "impl")$.valuesDeps$invalidate()
	#.subset2(DataInSets, "impl")$.values$remove(currentproject)
})
## save tables
observeEvent(input$results, {
	req(length(working_project()) > 0)
	req(DataInSets[[working_project()]]$data_results)
	results = DataInSets[[working_project()]]$data_results
	results[,sapply(results,is.numeric)] <- signif(results[,sapply(results,is.numeric)],3)
	saved_table$results <- results
})

observeEvent(input$sample, {
	req(length(working_project()) > 0)
	req(DataInSets[[working_project()]]$MetaData)
	saved_table$sample <- DataInSets[[working_project()]]$MetaData
})

observeEvent(input$data_wide, {
	req(length(working_project()) > 0)
	req(DataInSets[[working_project()]]$data_wide)
	data_wide <- DataInSets[[working_project()]]$data_wide %>%
	tibble::rownames_to_column(var = "ID")
	saved_table$data <- data_wide
})

observeEvent(input$ProteinGeneName, {
	req(length(working_project()) > 0)
	req(DataInSets[[working_project()]]$ProteinGeneName)
	saved_table$ProteinGeneName <- DataInSets[[working_project()]]$ProteinGeneName
})

### unit and short name
output$exp_unit <- renderUI({
	req(length(working_project()) > 0)
	textInput("exp_unit", "Expression Data Units", value="Expression Level", width="300px")
})

output$ShortName <- renderUI({
	req(length(working_project()) > 0)
	textInput("ShortName", "Rename Project Name", value="", width="300px")
})

observe({
	req(length(working_project()) > 0)
	updateTextInput(session, "exp_unit", value=DataInSets[[working_project()]]$exp_unit)
	updateTextInput(session, "ShortName", value=DataInSets[[working_project()]]$ShortName)
})

observeEvent(input$exp_unit, {
	req(length(working_project()) > 0)
	DataInSets[[working_project()]]$exp_unit <- input$exp_unit
})

observeEvent(input$ShortName, {
	req(length(working_project()) > 0)
	DataInSets[[working_project()]]$ShortName <- input$ShortName
})

###
output$loadedprojects <- renderUI({
	req(length(working_project()) > 0)
	radioButtons("current_dataset", label = "Change Working Dataset", choices=DS_names(), inline = F, selected=working_project())
})

observeEvent(input$current_dataset, {
	working_project(input$current_dataset)
})

output$project <- renderText({
	if (length(working_project()) == 0){
		""
	} else {
		paste("Project: ", working_project(), sep=" ")
	}
})

output$summary <- renderText({
	req(length(working_project()) > 0)

	summary=stringr::str_c('<style type="text/css">
		.disc {	list-style-type: disc;}
		.square { list-style-type: square; margin-left: -2em;	font-size: small}
		</style>',
		"<h2>Project ID: ", DataInSets[[working_project()]]$ProjectID, "</h2><br>",
		"<ul class='disc'>",
		"<li>Project Short Name: ", DataInSets[[working_project()]]$ShortName, "</li>",
		"<li>Description: ", DataInSets[[working_project()]]$Name, "</li>",
		"<li>Species: ", DataInSets[[working_project()]]$Species, "</li>",
		"<li>Data Path: ", DataInSets[[working_project()]]$Path, "</li>",
		"<li>Number of Samples: ", nrow(DataInSets[[working_project()]]$MetaData), "</li>",
		"<li>Number of Groups: ", length(DataInSets[[working_project()]]$groups), " (please see group table below)</li>",
		"<li>Number of Genes/Proteins: ", nrow(DataInSets[[working_project()]]$data_wide), "</li>",
		"<li>Number of Comparison Tests: ", length(DataInSets[[working_project()]]$tests), "</li>",
		'<ul class="square">', paste(stringr::str_c("<li>", DataInSets[[working_project()]]$tests, "</li>"), collapse=""), "</ul></li></ul><br><hr>",
		"<h4>Number of Samples in Each Group</h4>"
	)
})

group_info <- reactive({
	req(length(working_project()) > 0)
	req(DataInSets[[working_project()]]$MetaData)
	group_info <- DataInSets[[working_project()]]$MetaData %>% dplyr::group_by(group) %>% dplyr::count()
	return(t(group_info))
})
output$group_table <- renderTable(group_info(), colnames=F)

output$results <- DT::renderDataTable(server=TRUE,{
	req(length(working_project()) > 0)
	req(DataInSets[[working_project()]]$data_results)

	results <- DataInSets[[working_project()]]$data_results %>%
	#dplyr::select(-one_of(c("UniqueID","id")))
	dplyr::select(-any_of(c("id")))

	results[,sapply(results,is.numeric)] <- signif(results[,sapply(results,is.numeric)],3)

	DT::datatable(results,  extensions = 'Buttons',
		options = list(	dom = 'lBfrtip', pageLength = 15,
			buttons = list(
				list(extend = "csv", text = "Download Current Page", filename = "Page_Results",	exportOptions = list(modifier = list(page = "current")))
			)
		),
	rownames= T)
})

output$sample <-  DT::renderDT(server=FALSE, {
	req(length(working_project()) > 0)
	req(DataInSets[[working_project()]]$MetaData)

	meta <- DataInSets[[working_project()]]$MetaData
	DT::datatable(meta,  extensions = 'Buttons',
		options = list(dom = 'lBfrtip', pageLength = 15,
			buttons = list(
				list(extend = "csv", text = "Download Current Page", filename = "Page_Samples",	exportOptions = list(modifier = list(page = "current"))),
				list(extend = "csv", text = "Download All", filename = "All_Samples",	exportOptions = list(modifier = list(page = "all")
				)
			)
		)
	),
rownames= F)
})

output$comp_info <- renderUI ({
	req(length(working_project()) > 0)
	req(DataInSets[[working_project()]]$comparison)
	output$comparison <-  DT::renderDT(server=FALSE,{
		DT::datatable(DataInSets[[working_project()]]$comp_info,  extensions = 'Buttons',
			options = list(dom = 'lBfrtip', pageLength = 15,
				buttons = list(
					list(extend = "csv", text = "Download Current Page", filename = "Page_comp_info", exportOptions = list(modifier = list(page = "current"))),
					list(extend = "csv", text = "Download All", filename = "All_comp_info",	exportOptions = list(modifier = list(page = "all")))
				)
			)
		)
	}
)
tagList(
	#h4("Comparison Table (shown only when RData file contains comp_info)"),
	dataTableOutput('comparison')
)
})

output$data_wide <- DT::renderDataTable(server=TRUE, {
	req(length(working_project()) > 0)
	req(DataInSets[[working_project()]]$data_wide)
	data_w <- DataInSets[[working_project()]]$data_wide
	data_w[,sapply(data_w,is.numeric)] <- signif(data_w[,sapply(data_w,is.numeric)],3)

	DT::datatable(data_w, extensions = c('FixedColumns', 'Buttons'),
		options = list(
			pageLength = 15,	dom = 'lBfrtip', 	scrollX = TRUE,	fixedColumns = list(leftColumns = 1),
			buttons = list(
				list(extend = "csv", text = "Download Current Page", filename = "Page_Samples",	exportOptions = list(modifier = list(page = "current")))
			)
		)
	)
})

output$ProteinGeneName <- DT::renderDataTable(server=TRUE, {
	req(length(working_project()) > 0)
	req(DataInSets[[working_project()]]$ProteinGeneName)
	DT::datatable(DataInSets[[working_project()]]$ProteinGeneName, extensions = 'Buttons',
		options = list(dom = 'lBfrtip', pageLength = 15,
			buttons = list(
				list(extend = "csv", text = "Download Current Page", filename = "Page_GeneNames",	exportOptions = list(modifier = list(page = "current")))
			)
		),
	rownames= F)
})

## download big table
output$download_result_button <- shiny::downloadHandler(
	filename = function() {
		paste("Results-", Sys.Date(), ".csv", sep="")
	},
	content = function(file) {
		write.csv(DataInSets[[working_project()]]$data_results, file)
	}
)

output$download_data_button <- shiny::downloadHandler(
	filename = function() {
		paste("Data-", Sys.Date(), ".csv", sep="")
	},
	content = function(file) {
		write.csv(DataInSets[[working_project()]]$data_wide, file)
	}
)

output$download_ProteinGeneName_button <- shiny::downloadHandler(
	filename = function() {
		paste("ProteinGeneName-", Sys.Date(), ".csv", sep="")
	},
	content = function(file) {
		write.csv(DataInSets[[working_project()]]$ProteinGeneName, file)
	}
)

