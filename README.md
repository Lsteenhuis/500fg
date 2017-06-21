This repo contains scripts which were used for the identification of immune components which are important during disease.
R scripts were made to do different analyses, such as the calculation of a connectivity score, correlation analyses, and differential expression analyses.

# ROADMAP

## ./gcmap
These are the scripts which were used to calculate the connectivity score with two different methods: GCMAP[1] and kidd et al [2] 


  	connectivity plotter <- code to create plots for connectivity scores (kidd et al)
	scoreCalculator <- calculates connectivity score ( kidd et al method)
   	clusMap <- creates plot for calculated Pertubation combinations (gcmap)
	createHeatMapWilcox <- creates heatmap for the wilcox results (gcmap)
	cscore.heatmap <- heatmap of cmap connectivity scores (gcmap)
	drugDeOverlap <- calculates overlap between DE scores and drug gene profile (gcmap)
	calculateEnrichmentScore <- code which uses enrichmentFunctions
	enrichmentFunctions <- functions for calculations of gcMAP
	keggToEns.R <- sets kegg gene ids to ensembl gene ids
	plotProbPadj.R <- calculates probability for gcmap
	runDeSeq.R <- runs deSeq for IT samples
	hyperDistTest <- small test file for hyper distribution test
	plotProbPadj <- create plots for probability scores

## ./hala
The scripts in this folder were used to create the tables for the HALLA program and the plotting of some results.


	halla.data.proc <- create halla tables
	manuscript_plot_edit <- create plots of halla

## ./prs
The scripts in this folder were used to perform correlation analyses between the polygenic risk scores and immmune components.


	correlation.heatmap <- create heatmaps of correlations
 	correlation.multiplots <- create barplots of correations
	correlation_functions <- function to creation of correlations
	correlationMatrixMaker 
	create_halla_tables
	halla_multiple_thresholds <- create plots of each threshold
	heatmap.theses <- create plots for thesis
	nonImmuneDisease_test <- create plots for nonImmuneDiseases
	prsSignificantyTest.R <- small fisher test for prs testing
	prs_aid_vs_cytokine_pvalue_checker <- creating plots of enrichtment
	prs_aid_vs_cytokine_lotter <- plots aid vs cytokine
	violin_plots <- creating of violin_plot


[1] Lamb J, et al, The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease., Scinece, 2006,1929-1935
[2] kidd et al, Mapping the effects of drugs on the immune system, nature biotechnology,34, 2016, 47-5