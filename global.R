###########################################################################################################
## Proteomics Visualization R Shiny App
##
##This software belongs to Biogen Inc. All right reserved.
##
##@file: global.R
##@Developer : Benbo Gao (benbo.gao@Biogen.com)
##@Date : 05/23/2023
##@version 3.0
###########################################################################################################
options(stringsAsFactors=F)
options(ggrepel.max.overlaps = Inf)
options(shiny.maxRequestSize = 120*1024^2)  #upload files up to 100 Mb
options(rgl.useNULL = TRUE)

suppressPackageStartupMessages({
	library(shiny)
	library(shinythemes)
	library(shinyalert)
	library(shinyjqui)
        library(shinyjs)
        #library(tidyverse)
	library(dplyr)
	library(purrr)
	library(tidyr)
	library(tibble)
	library(plotly)
	library(ggpubr)
	library(ggrepel)
	library(DT)
	library(RColorBrewer)
	library(colourpicker)
	library(stringr)
	library(ggrastr)
	library(ggpmisc)
        library(rclipboard)
})

resultfilter <- function(results_long, test_sel, p_sel, direction, pvalcut, FCcut, sel_label) {
	results_long_filtered <- results_long %>%
	dplyr::mutate_if(is.factor, as.character) %>%
	dplyr::filter(if (!is.na(test_sel)) {test == test_sel} else TRUE) %>%
	dplyr::filter(if (p_sel=="Padj") {Adj.P.Value < pvalcut} else {P.Value < pvalcut}) %>%
	dplyr::filter(if (direction=="Up") {logFC >= FCcut} else if (direction=="Down") {logFC <= -FCcut} else {abs(logFC) >= FCcut}) %>%
	dplyr::arrange(desc(abs(logFC)))  %>%
	dplyr::mutate(labelid = !!sym(sel_label)) %>%
	dplyr::filter(!is.na(labelid) & labelid != "") %>%
	#dplyr::distinct(labelid, .keep_all = TRUE) %>%
	as.data.frame()

	sig_genes_FC <- results_long_filtered %>%
	dplyr::filter(!is.na(Gene.Name) & Gene.Name != "") %>%
	dplyr::group_by(Gene.Name) %>% # if duplicated gene, select one with max(abs(logFc))
	dplyr::arrange(desc(abs(logFC))) %>%
	dplyr::slice(1) %>% #don't use dplyr::top
	dplyr::ungroup() %>%
	dplyr::arrange(desc(abs(logFC))) %>%
	dplyr::pull(logFC, Gene.Name)

	sig_genes <- names(sig_genes_FC)

	all_genes <- results_long %>%
	dplyr::filter(!is.na(`Gene.Name`) & Gene.Name != "") %>%
	dplyr::pull(Gene.Name) %>%
	unique()

	return(
		list(
			"results_long_filtered" = results_long_filtered,
			"sig_genes" = sig_genes,
			"all_genes" = all_genes,
			"sig_genes_FC" = sig_genes_FC
		)
	)
}

GeneFilter <- function(results_long, test_sel, p_sel, direction, pvalcut, FCcut, sel_label) {

	results_long_filtered <- results_long %>%
	dplyr::mutate_if(is.factor, as.character) %>%
	dplyr::filter(if (!is.na(test_sel)) {test == test_sel} else TRUE) %>%
	#dplyr::filter(if (!is.na(test_sel)) {test %in% test_sel} else TRUE) %>%
	dplyr::filter(if (p_sel=="Padj") {Adj.P.Value < pvalcut} else {P.Value < pvalcut}) %>%
	dplyr::filter(if (direction=="Up") {logFC >= FCcut} else if (direction=="Down") {logFC <= -FCcut} else {abs(logFC) >= FCcut}) %>%
	dplyr::arrange(desc(abs(logFC)))  %>%
	dplyr::mutate(labelid = !!sym(sel_label)) %>%
	dplyr::filter(!is.na(labelid) & labelid != "") %>%
	#dplyr::distinct(labelid, .keep_all = TRUE) %>%
	as.data.frame()

	return(results_long_filtered)
}

datafilter <- function(data_long, sel_group, sel_sample, sel_ID="all", sel_label) {

	data_long_filtered  <- data_long %>%
	dplyr::mutate_if(is.factor, as.character) %>%
	dplyr::filter((group %in% sel_group) & (sampleid %in% sel_sample)) %>%
	dplyr::mutate(group = factor(group, levels = sel_group)) %>%
	dplyr::mutate(sampleid = factor(sampleid, levels = sel_sample)) %>%
	dplyr::filter(!is.na(expr))  %>%
	dplyr::mutate(labelid = !!sym(sel_label)) %>%
	dplyr::filter(!is.na(labelid) & labelid != "") %>%
	dplyr::filter(if (sel_ID[1] !="all")  {labelid %in% sel_ID} else TRUE)

	#labelidcount <- data_long_filtered %>%
	#group_by(sampleid, labelid) %>%
	#summarise(count = n(),.groups = "drop")

	if (sel_label == "Gene.Name") {
		data_long_filtered  <- data_long_filtered %>%
		dplyr::group_by(sampleid, labelid) %>%
		dplyr::mutate(expr = mean(expr, na.rm = TRUE)) %>%
		dplyr::slice(1) %>%
		dplyr::ungroup()
	}

	#labelidcount <- data_long_filtered %>%
	#group_by(sampleid, labelid) %>%
	#summarise(count = n(),.groups = "drop")

	data_wide_filtered <- data_long_filtered  %>%
	dplyr::select(one_of(c("labelid", "expr", "sampleid"))) %>%
	tidyr::spread(sampleid, expr, fill = NA) %>%
	tibble::column_to_rownames(var = "labelid")

	data_long_mean <- data_long_filtered %>%
	dplyr::group_by(., group, labelid) %>%
	dplyr::summarise(mean=mean(expr, na.rm = TRUE), .groups = 'drop') %>%
	dplyr::ungroup()

	data_wide_mean <-  data_long_mean  %>%
	tidyr::spread(.,group, mean, fill = NA) %>%
	remove_rownames(.) %>%
	column_to_rownames(.,var="labelid")

	return(
		list(
			"data_long_filtered" = data_long_filtered,
			"data_wide_filtered"= data_wide_filtered,
			"data_long_mean" = data_long_mean,
			"data_wide_mean"= data_wide_mean
		)
	)
}

named_group_split <- function(.tbl, ...) {
	grouped <- group_by(.tbl, ...)
	names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / ")))
	grouped %>%	group_split() %>%	rlang::set_names(names)
}


LoadedData <- reactiveValues()

GetGeneSetNames <- function() {

	if(!("gmtlist" %in% names(LoadedData))){
		load("db/gmtlist.RData")
		LoadedData[["gmtlist"]]  <-  gmtlist
	} else {
		gmtlist <- LoadedData[["gmtlist"]]
	}

	if(!("kegg.pathways" %in% names(LoadedData))){
		load("db/kegg.pathways.RData")
		LoadedData[["kegg.pathways"]]  <-  kegg.pathways
	} else{
		kegg.pathways <- LoadedData[["kegg.pathways"]]
	}

	genesetnames <- c()
	for (setname in names(gmtlist)) {
		genesetnames <- c(genesetnames, names(gmtlist[[setname]]))
	}
	genesetnames <- c(names(kegg.pathways$human$kg.sets), genesetnames)
	return(genesetnames)
}

GetGenesFromGeneSet <- function(sel_geneset) {

	if(!("hgnc" %in% names(LoadedData))){
		load("db/hgnc.RData")
		LoadedData[["hgnc"]]  <-  hgnc
	} else {
		hgnc <- LoadedData[["hgnc"]]
	}
	kegg.pathways <- LoadedData[["kegg.pathways"]]
	gmtlist <- LoadedData[["gmtlist"]]

	if (sel_geneset %in% names(kegg.pathways$human$kg.sets)) {
		geneset_genes <- kegg.pathways$human$kg.sets[[sel_geneset]]
	}	else {
		for (setname in names(gmtlist)) {
			if (sel_geneset %in% names(gmtlist[[setname]])) {
				geneset_genes <- gmtlist[[setname]][[sel_geneset]]
				break
			}
		}
	}

	geneset_genenames <- hgnc %>%
	dplyr::filter(entrez_id %in% geneset_genes) %>%
	dplyr::pull(symbol)
	return(geneset_genenames)
}

ProcessUploadGeneList <- function(gene_list) {
	if(grepl("\n", gene_list)) {
		gene_list <-  stringr::str_split(gene_list, "\n")[[1]]
	} else if (grepl(",",gene_list)) {
		gene_list <-  stringr::str_split(gene_list, ",")[[1]]
	}
	gene_list <- gsub(" ", "", gene_list, fixed = TRUE)
	gene_list <- unique(gene_list[gene_list != ""])
	return(gene_list)
}


homologs=readRDS("db/Homologs.rds")
homolog_mapping<-function(genelist, species1, species2, homologs) {
  if (species2=="human") {
    genelist2=toupper(genelist)
  } else {
    genelist2=str_to_title(genelist)
  }
  if (tolower(species1) %in% c("human", "mouse", "rat")) {
    species1=tolower(species1)
    lookup<-homologs%>%filter(Species1==species1, Species2==species2)
    df<-data.frame(genelist, genelist2)%>%left_join(lookup%>%transmute(genelist=symbol1, mapped_symbol=symbol2))%>%
      mutate(mapped_symbol=ifelse(is.na(mapped_symbol), genelist2, mapped_symbol))
    genelist2<-df$mapped_symbol
  }
  return(genelist2)
}


UserColorPlalette <- function (colpalette, items){
	n <- brewer.pal.info[colpalette,'maxcolors']
	if (length(unique(items)) <=  n) {
	user_color <- RColorBrewer::brewer.pal(length(unique(items)), colpalette)
	} else {
		user_color <- colorRampPalette(RColorBrewer::brewer.pal(n, colpalette))(length(unique(items)))
	}
	return(user_color)
}


mycss <- "select ~ .selectize-control .selectize-input { max-height: 200px;  overflow-y: auto;}
.shiny-notification{position: fixed; top: 33%; left: 33%; right: 33%; }"

config=NULL
server_dir=NULL
test_dir=NULL
gmt_file_info=NULL
if (file.exists("config.csv")) { #load optional configuration file
	config=read_csv("config.csv")
	N=match("server_dir", config$category)
	if (!is.na(N)) {server_dir=config$value[N]}
	N=match("test_dir", config$category)
	if (!is.na(N)) {test_dir=config$value[N]}
	N=match("gmt_file_info", config$category)
	if (!is.na(N)) {gmt_file_info=config$value[N]}
	#browser() #debug
}

footer_text = '
<hr>
<div align="center" style="font-size:11px">xOmicsShiny Version 1.0</div>
'
