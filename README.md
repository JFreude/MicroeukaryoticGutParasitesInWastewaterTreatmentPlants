# Microeukaryotic gut parasites in wastewater treatment plants: diversity, activity, and removal
This repository is a collection of scripts that guide you through the data analyses of the paper of [Freudenthal et al., 2022]().  
The original data (the raw OTU tables with taxonomic annotation, the metadata, and the trait database we used to assign the parasites) are provided in the folder [OriginalData](OriginalData/). Intermediate files, as well as results and figures, can be obtained by running the data analysis scripts.

## Data processing
In this [pipeline]( R/01_DataPreprocessing.md) you find all steps that were performed the process the raw OUT data. Within this pipeline, you also find links to the scripts for  
*	[Supplementary Figure 1]( R/02_MultivariateDispersionAndBetaDiversity.md) - Comparing microbial communities between WWTP locations to identify outliers  
*	[Supplementary Figure 2]( R/03_VariationsCausedBySampleProcessing.md) - Assessment of the variation caused by sampling processing (sequencing)  
*	[Supplementary Figure 3]( R/04_RarefactionCurves.md) - Rarefaction curves for metagenomic (rDNA) and metatranscriptomic (rRNA) data  
*	[Supplementary Table 1]( R/05_OverviewCommunityComposition.md) - Microbial community composition after quality filtering   
*	[Supplementary Table 2]( R/06_ParasiticGenera.md) Parasitic genera in WWTPs based on both metagenomic and metatranscriptomic data  

## Area plots, line plots, and box plots
In this section, we provide scripts for the visualization of the changes in microbial community composition during wastewater treatment ([Figure 1]( R/07_MicrobialCommunityComposition.md)), the total number of rDNA and rRNA sequences ([Supplementary Figure 4]( R/09_TotalNumberSequences.md)), the removal of parasitic protists from wastewater ([Figure 2]( R/10_RemovalOfParasiticProtistsFromWastewater.md)) and the differences between measurable presence and activity of parasitic protists in the inflow ([Figure 3]( R/11_DetectionOfParasiticProtistsInWastewater_AbundanceVersusActivity.md)). Further, we show how to test for significant differences between the treatment compartments ([PERMANOVA]( R/08_PERMANOVA.md)).

## Network inference
Here we show how we conducted network analyses using SPIEC-EASI and SparCC ([Figure 3]( R/14_CoOccurrenceNetworksOfParasiticOrdersInTheInflowAndDenitrificationBioreactor.md)). Furthermore, we show how the NMDS plots for the WWTPs as a whole ([Supplementary Figure 5]( R/12_MicrobialCommunityStructureAndEnvironmentalFactorsAcrossWWTPs.md)) and for the individual compartment types ([Supplementary Figure 6]( R/13_MicrobialCommunityStructureAndEnvironmentalFactorsInTheSeparateWWTPCompartments.md)) were calculated.
