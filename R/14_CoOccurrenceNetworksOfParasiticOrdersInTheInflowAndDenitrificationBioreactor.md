Co-occurrence networks of parasitic orders in the four WWTP compartments
================
Jule Freudenthal
2022-02-03

**R version:** 3.6.2 (2019-12-12), Dark and Stormy Night  
**Packages**

-   docstring v. 1.0.0
-   dplyr v. 1.0.7  
-   purrr v. 0.3.4  
-   reshape2 v. 1.4.4  
-   rlist v. 0.4.6.1  
-   SpiecEasi v. 1.1.1  
-   zeallot v. 0.1.0

``` r
# Load packages
library(docstring)
library(dplyr)
library(purrr)
library(reshape2)
library(rlist)
library(SpiecEasi)
library(zeallot)

# Load functions
source("./Functions/association.frequencies.R")
```

Co-occurrence network analyses are performed to assess the complexity of
the correlations between the microbial community and
parasites.Initially, two pre-processing steps are conducted to reduce
indirect associations (spurious edges).

To evaluate the influence of measured environmental data on the
microbial community, NMDS plots for the entire WWTPs as well as for each
compartment individually are computed. The environmental data are fitted
onto the ordinations using envfit.

# 01 Prevalence filter

It has been shown that co-absence can yield high correlation values.
Hence, taxa only present in up to three samples of one compartment
across locations are excluded.

``` r
# Create new directory for Network data
if(!file.exists("./DataToAnalyse")){dir.create("./DataToAnalyse")}
if(!file.exists("./DataToAnalyse/NetworkData")){dir.create("./DataToAnalyse/NetworkData")}
if(!file.exists("./DataToAnalyse/NetworkData/Raw")){dir.create("./DataToAnalyse/NetworkData/Raw")}

# Loop iterates over the data types (DNA/RNA data)
for(data_type in c("DNA", "RNA")){
  
  # Load count- and meta- and taxonomy data
  load(paste0("./DataToAnalyse/RData/", data_type, "_data.rds"))
  
  # Loop iterates over the WWTP compartments 
  for(compartment in c("INF", "DNF", "NFC", "EFF")){

    # Subset data, keep only samples belonging to one compartment
    counts_comp <- counts[,grepl(compartment, colnames(counts))]
   
    # Convert into presence/absence matrix
    counts_comp_copy <- counts_comp
    counts_comp_copy[counts_comp_copy > 0] <- 1
    
    # Row sum of each taxon = number of occurrences across samples
    to_filter <- which(rowSums(counts_comp_copy) <= round(1/3*ncol(counts_comp_copy)))
    print(paste("OTUs present in <=", round(1/3*ncol(counts_comp_copy)), "of", ncol(counts_comp), 
                "samples of the", compartment, "compartment are summed into one pseudo taxon."))
 
    # Filter count data
    print(paste("Reads of", length(to_filter)-nrow(counts_comp[rowSums(counts_comp) == 0,]), "out of", 
                nrow(counts_comp[rowSums(counts_comp) != 0,]), "OTUs are summed into one pseudo taxon."))
    counts_comp_filtered <- counts_comp[-to_filter,]
    counts_comp_filtered <- rbind(counts_comp_filtered,
                                  t(data.frame(filtered_taxa=colSums(counts_comp[to_filter,]))))

    # Export abundances as txt
    write.table(counts_comp_filtered,
                file = paste0("./DataToAnalyse/NetworkData/Raw/", data_type, "_", compartment, ".txt"), 
                sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
  }
}
```

    ## [1] "OTUs present in <= 3 of 9 samples of the INF compartment are summed into one pseudo taxon."
    ## [1] "Reads of 452 out of 969 OTUs are summed into one pseudo taxon."
    ## [1] "OTUs present in <= 3 of 10 samples of the DNF compartment are summed into one pseudo taxon."
    ## [1] "Reads of 785 out of 1486 OTUs are summed into one pseudo taxon."
    ## [1] "OTUs present in <= 3 of 10 samples of the NFC compartment are summed into one pseudo taxon."
    ## [1] "Reads of 750 out of 1439 OTUs are summed into one pseudo taxon."
    ## [1] "OTUs present in <= 3 of 9 samples of the EFF compartment are summed into one pseudo taxon."
    ## [1] "Reads of 788 out of 1392 OTUs are summed into one pseudo taxon."
    ## [1] "OTUs present in <= 3 of 10 samples of the INF compartment are summed into one pseudo taxon."
    ## [1] "Reads of 572 out of 1148 OTUs are summed into one pseudo taxon."
    ## [1] "OTUs present in <= 3 of 10 samples of the DNF compartment are summed into one pseudo taxon."
    ## [1] "Reads of 709 out of 1236 OTUs are summed into one pseudo taxon."
    ## [1] "OTUs present in <= 3 of 9 samples of the NFC compartment are summed into one pseudo taxon."
    ## [1] "Reads of 668 out of 1114 OTUs are summed into one pseudo taxon."
    ## [1] "OTUs present in <= 3 of 10 samples of the EFF compartment are summed into one pseudo taxon."
    ## [1] "Reads of 804 out of 1317 OTUs are summed into one pseudo taxon."

# 02 Sparse Inverse Covariance for Ecological Statistical Inference (SPIEC-EASI)

``` r
# Create new directory for Network data
if(!file.exists("./Results/Networks")){dir.create("./Results/Networks")}
if(!file.exists("./Results/Networks/SpiecEasi")){dir.create("./Results/Networks/SpiecEasi")}

sink('./Results/Networks/SpiecEasi/SpiecEasi_overview.txt')

# Loop iterates over the data types (DNA/RNA data)
for(data_type in c("DNA", "RNA")){

  # Loop iterates over the WWTP compartments 
  for(compartment in c("INF", "DNF", "NFC", "EFF")){
    
    print(paste(data_type, compartment))

    # Load count- and meta- and taxonomy data
    counts <- read.table(paste0("./DataToAnalyse/NetworkData/Raw/", data_type, "_", compartment, ".txt"), 
                         sep="\t", header = TRUE, row.names = 1)
    
    # Step 1: Calculate network
    spiec_out <- spiec.easi(t(counts), method='mb', lambda.min.ratio=0.001, 
                            nlambda=50, pulsar.params=list(rep.num=20, seed=1244,
                                                           thresh=0.05))
  
    # Step 2: Extract edge weights
    betaMat <- as.matrix(symBeta(getOptBeta(spiec_out)))
    colnames(betaMat) <- row.names(betaMat) <- colnames(spiec_out$est$data)
    betaMat[lower.tri(betaMat, diag = T)] <- 0
    
    # Step 3: Transform SPIEC-EASI results for Cytoscape
    graph.table <- melt(betaMat)
    graph.table <- graph.table[graph.table$value != 0,]
    colnames(graph.table) <- c("Node 1", "Node 2", "Edge weights")
    
    # Remove binned taxa from table (resulting from prevalence filter)
    graph.table <- graph.table[which(!(graph.table$`Node 1` == 'filtered_taxa' | 
                                         graph.table$`Node 2` == 'filtered_taxa')),]
    
    # Get number of negative and positive edges
    print(paste0("Number of positive edges: ", sum(graph.table$`Edge weights` > 0)))
    print(paste0("Number of negative edges: ", sum(graph.table$`Edge weights` < 0)))
    print(paste0("Total number of edges: ", nrow(graph.table)))
    #print(paste0("check how far we are from the target stability threshold (0.05): ", 
    #             getStability(spiec_out)))
    
    # Create new column specifying if correlation is positive or negative
    graph.table$Association <- "negative"
    graph.table$Association[graph.table$`Edge weights` > 0] <- "positive"

    # Export data as txt file
    write.table(graph.table, file = paste0("./Results/Networks/SpiecEasi/", data_type, "_",
                                           compartment, ".txt"), sep="\t", quote=FALSE, row.names=FALSE)
  }
}
sink()
```

# 03 Sparse Correlations for Compositional data (SparCC, python)

**Attention:** SparCC is a python package that runs on the command line.
For visualization the code is given.

-   **Step 1:** Compute correlations  
-   **Step 2:** Compute bootstraps  
-   **Step 3:** Compute pseudo p-values

<!-- -->

    ## echo OFF
    ## 
    ## :: Define pathways for input files, outputfiles and to python
    ## set SparCCPackage=C:\Users\Jule\anaconda3\pkgs\sparcc-0.1.0-0\python-scripts\
    ## set InputPath=C:\Users\Jule\sciebo\AG_Bonkowski\Projects\ContrastingActivityPatterns\DataToAnalyse\NetworkData\Raw\
    ## set OutputPath=C:\Users\Jule\sciebo\AG_Bonkowski\Projects\ContrastingActivityPatterns\DataToAnalyse\NetworkData\SparCC\
    ## 
    ## :: Create directory for output files
    ## if not exist "%OutputPath%" (
    ##    mkdir "%OutputPath%"
    ## )
    ## 
    ## :: Define Data type and WWTP compartments
    ## set data_type=DNA,RNA
    ## set compartments=INF,DNF,NFC,EFF
    ## 
    ## :: Loop iterates over the data types (DNA/RNA)
    ## for %%d in (%data_type:,= %) do (
    ##  echo data type: %%d
    ##  
    ##  :: Loop iterates over the WWTP compartments
    ##  for %%c in (%compartments:,= %) do (
    ##  
    ##      echo compartment: %%c
    ## 
    ##      echo Step 1: Compute correlations
    ##      
    ##      python %SparCCPackage%SparCC.py %InputPath%%%d_%%c.txt -i 20 -c %OutputPath%%%d_%%c_SparCC.txt > %OutputPath%%%d_%%c_Logfile_corr.log
    ## 
    ##      echo Step 2: Re-sample data
    ##      
    ##      if not exist %OutputPath%Resampling\%%d_%%c (
    ##          mkdir %OutputPath%Resampling\%%d_%%c
    ##      )
    ##      
    ##      python %SparCCPackage%MakeBootstraps.py %InputPath%%%d_%%c.txt -n 100 -t permutation_#.txt -p %OutputPath%Resampling\%%d_%%c\
    ##      
    ##      echo Step 3: Compute correlations with re-sampled data
    ## 
    ##      if not exist %OutputPath%Bootstraps\%%d_%%c (
    ##          mkdir %OutputPath%Bootstraps\%%d_%%c
    ##      )
    ## 
    ##      for /l %%a in (0,1,99) do (
    ##          python %SparCCPackage%SparCC.py %OutputPath%Resampling\%%d_%%c\permutation_%%a.txt -c %OutputPath%Bootstraps\%%d_%%c\bootstraps_%%a.txt >> %OutputPath%%%d_%%c_Logfile_bootstraps.log
    ##      )
    ## 
    ##      echo Step 4: Compute pseudo p-values
    ##      python %SparCCPackage%PseudoPvals.py %OutputPath%%%d_%%c_SparCC.txt %OutputPath%Bootstraps\%%d_%%c\bootstraps_#.txt 100 -o %OutputPath%%%d_%%c_pvals_two_sided.txt -t two_sided >> %OutputPath%%%d_%%c_Logfile_pvals.log
    ## 
    ##  )
    ## )

-   **Step 4:** Compute FDR for edges (R)

``` r
# Create new directory for Network data
if(!file.exists("./Results/Networks/SparCC")){dir.create("./Results/Networks/SparCC")}

# Get file names of associations and p_values
filenames_edges <- list.files("./DataToAnalyse/NetworkData/SparCC/", pattern=".*SparCC.txt", full.names=F)
filenames_pval <- list.files("./DataToAnalyse/NetworkData/SparCC/",  
                             pattern=".*pvals_two_sided.txt", full.names=F)

# Check if correlations and p-values are in corresponding order
if(!all.equal(gsub("_SparCC.txt","",filenames_edges),
              gsub("_pvals_two_sided.txt","",filenames_pval))){
  stop("Correlations and p-values are not in corresponding order")
}

sink('./Results/Networks/SparCC/SparCC_overview.txt')

# Loop iterates over the file names
for(filenames in seq(filenames_edges)){
  print(gsub("_SparCC.txt","",filenames_edges[filenames]))
    
  # Load SparCC associations and corresponding pseudo p-values
  edges <- read.table(paste0("./DataToAnalyse/NetworkData/SparCC/", filenames_edges[filenames]),
                      header = T, row.names = 1, sep="\t", check.names=FALSE)
  pvals <- read.table(paste0("./DataToAnalyse/NetworkData/SparCC/", filenames_pval[filenames]),
                      header = T, row.names = 1, sep="\t", check.names=FALSE)
  
  # Check if the row and column names of both tables match
  if(!all.equal(rownames(edges), rownames(pvals), colnames(edges), colnames(pvals))){
    stop("Row and Column names of the tables must match")
  }
  
  # Set p-values of 0 to a non-zero, small p-value to take the logarithm
  pvals[pvals==0] <- 0.000000001
  
  # Log-transform p-values
  log.pval <- -1*log10(pvals) 
  
  # Remove all edges with significance below/equal 1.3 = FDR 0.05
  log.pval[log.pval<=1.3] <- 0
  
  # Filter correlations
  edges[log.pval == 0] <- NA
  
  # Reduce symmetrical data set to upper triangle 
  edges[lower.tri(edges, diag = T)] <- NA
  
  # Transform Correlations for Cytoscape
  edges <- melt(as.matrix(edges))
  edges <- edges[!is.na(edges$value),]
  colnames(edges) <- c("Node 1", "Node 2", "Edge weights")
  
  # Remove binned taxa from table (from prevalence filter)
  edges <- edges[which(!(edges$`Node 1` == 'filtered_taxa' | edges$`Node 2` == 'filtered_taxa')),]
  
  # Get number of negative and positive edges
  print(paste0("Number of positive edges: ", sum(edges$`Edge weights` > 0)))
  print(paste0("Number of negative edges: ", sum(edges$`Edge weights` < 0)))
  print(paste0("Total number of edges: ", nrow(edges)))
  
  # Create new column specifying if correlation is positive or negative
  edges$Association <- "negative"
  edges$Association[edges$`Edge weights` > 0] <- "positive"
    
  # Export data as txt file
  write.table(edges, file = paste0("./Results/Networks/SparCC/", 
                                   gsub("_SparCC", "", filenames_edges[filenames])),
              sep="\t", quote=FALSE, row.names=FALSE)
}
sink()
```

# 04 Combine SpiecEasi and SparCC

``` r
# Create new directory for Network data
if(!file.exists("./Results/Networks/CompSparSpiec")){
  dir.create("./Results/Networks/CompSparSpiec")}

sink('./Results/Networks/CompSparSpiec/CompSparSpiec.txt')

# Loop iterates over the data types (DNA/RNA data)
for(data_type in c("DNA", "RNA")){
  
  # Loop iterates over the WWTP compartments 
  for(compartment in c("INF", "DNF", "NFC", "EFF")){

    print(paste(data_type, compartment))
    
    # Load edge tables
    spiec <- read.table(paste0("./Results/Networks/SpiecEasi/", data_type, "_",
                               compartment, ".txt"), header = T, sep = "\t", dec = ".")
    spar <- read.table(paste0("./Results/Networks/SparCC/", data_type, "_",
                              compartment, ".txt"), header = T, sep = "\t", dec = ".")
    
    # Combine Networks, keep SparCC correlation values
    edges1 <- which(paste0(spar$Node.1, spar$Node.2, spar$Association) %in% 
                    paste0(spiec$Node.1, spiec$Node.2, spiec$Association))
    edges2 <- which(paste0(spar$Node.2, spar$Node.1, spar$Association) %in% 
                    paste0(spiec$Node.1, spiec$Node.2, spiec$Association))
    spar_spiec <- spar[unique(c(edges1,edges2)),]
    colnames(spar_spiec) <- gsub("\\."," ",colnames(spar_spiec))
    rownames(spar_spiec) <- NULL
    
    # Get number of negative and positive edges
    print(paste0("Number of positive edges: ", sum(spar_spiec$`Edge weights` > 0)))
    print(paste0("Number of negative edges: ", sum(spar_spiec$`Edge weights` < 0)))
    print(paste0("Total number of edges: ", nrow(spar_spiec)))
    
    # Export edges
    write.table(spar_spiec, file = paste0("./Results/Networks/CompSparSpiec/", 
                                          data_type, "_", compartment, ".txt"), 
                sep="\t", quote=FALSE, row.names=FALSE)
  }
}
sink()
```

# 05 Transform taxonomy data for Cytoscape

**spiec.easi** uses the row names of the abundances as Source - and
Target node. To be able to combine taxonomy and network data in
**Cytoscape** one has to transform the taxonomy table by setting the row
names of the abundance data as first column. This column is used as
‘key’-to assign the taxonomy to the network table. Since the row names
of the abundance and taxonomy data match, one can also use the row names
of the taxonomy data as Key column.

``` r
# Loop iterates over the data types (DNA/RNA data)
for(data_type in c("DNA", "RNA")){
  
  # Load count- and meta- and taxonomy data
  load(paste0("./DataToAnalyse/RData/", data_type, "_data.rds"))
  
  # Loop iterates over the WWTP compartments 
  for(compartment in c("INF", "DNF", "NFC", "EFF")){
  
    # Subset data, keep only samples belonging to one compartment and taxa occurring in that compartment
    counts_comp <- counts[,grepl(compartment, colnames(counts))]
    counts_comp <- counts_comp[rowSums(counts_comp) != 0,]

    # Subset taxonomy
    taxonomy_comp <- taxonomy[rownames(taxonomy) %in% rownames(counts_comp),]
    
    # Assign new column to taxonomy with number of reads
    if(all.equal(rownames(counts_comp), rownames(taxonomy_comp))){
      taxonomy_comp$`Number of reads` <- rowSums(counts_comp)
    } else {
      stop("Rownames of count and taxonomy data are not equal")
    }

    # Change _X to undetermined
    taxonomy_comp <- data.frame(apply(taxonomy_comp, MARGIN = 2, 
                                      function(x) gsub("_X.*", " undetermined", x)))
  
    # Create key column 
    taxonomy_comp <- cbind(Key=rownames(taxonomy_comp), taxonomy_comp)
  
    # Export table for SpiecEasy and SparCC and comparisons
    write.table(taxonomy_comp, file = paste0("./Results/Networks/CompSparSpiec/", data_type,
                                             "_", compartment, "_taxonomy.txt"), 
                sep="\t", quote=FALSE, row.names=FALSE)
  }
}
```

# 06 Summaries taxonomy and edge tables at order level

The frequencies of correlations per taxonomic level is counted. As an
example, when calculating the frequencies of correlations per order, the
number of both negative and positive correlations at OTU level (here
Genus) are counted for each order.

``` r
# Create a reference taxonomy with unique order names 
# (we sometimes list the same orders twice as we distinguish between free-living and parasitic taxa)

# Load DNA and RNA taxonomies
load(paste0("./DataToAnalyse/RData/DNA_data.rds"))
DNA_taxonomy <- taxonomy
load(paste0("./DataToAnalyse/RData/RNA_data.rds"))
RNA_taxonomy <- taxonomy

# Merge both taxonomies
reference_taxonomy <- full_join(DNA_taxonomy, RNA_taxonomy)

# Change _X to undetermined
reference_taxonomy <- data.frame(apply(reference_taxonomy, MARGIN = 2, 
                                       function(x) gsub("_X.*", " undetermined", x)))

# Copy the taxonomy without family and genus level, then exclude all duplicated rows 
reference_taxonomy_copy <- reference_taxonomy[,!colnames(reference_taxonomy) %in% c("Family", "Genus")]
reference_taxonomy_copy <- unique(reference_taxonomy_copy)

# Assign numbers to the order names to make them unique
reference_taxonomy_copy <- reference_taxonomy_copy[order(reference_taxonomy_copy$Order),]
reference_taxonomy_copy$NewKey <- paste(seq(nrow(reference_taxonomy_copy)), 
                                        reference_taxonomy_copy$Order, sep = "_")
```

``` r
# Create new directory for Network data
if(!file.exists("./Results/Networks/CompSparSpiec/Order")){
  dir.create("./Results/Networks/CompSparSpiec/Order")}

# Loop iterates over the data types (DNA/RNA data)
for(data_type in c("DNA", "RNA")){

  # Loop iterates over the WWTP compartments 
  for(compartment in c("INF", "DNF", "NFC", "EFF")){

    # Load edges and taxonomy table of respective WWTP compartment and data set 
    taxonomy <- read.table(paste0("./Results/Networks/CompSparSpiec/", data_type, "_", 
                                  compartment, "_taxonomy.txt"), header = T, sep = "\t", dec = ".")
    spar_spiec <- read.table(paste0("./Results/Networks/CompSparSpiec/", data_type, "_", 
                                  compartment, ".txt"), header = T, sep = "\t", dec = ".")
    colnames(spar_spiec) <- gsub("\\."," ",colnames(spar_spiec))
    
    # Assign new unique key column to taxonomy
    taxonomy <- left_join(taxonomy, reference_taxonomy_copy)

    # Calculate frequencies per Order
    c(edges,taxonomy) %<-% 
      association.frequencies(edge.table=spar_spiec[,!colnames(spar_spiec) %in% "Association"], 
      taxonomy=taxonomy, old.key="Key", new.key="NewKey")
  
    # Delete NewKey column
    taxonomy$NewKey <- NULL
    
    # Export edges
    write.table(edges, file = paste0("./Results/Networks/CompSparSpiec/Order/", data_type, "_", 
                                     compartment, ".txt"), sep="\t", quote=FALSE, row.names=FALSE)
    
    # Export taxonomy
    write.table(taxonomy, file = paste0("./Results/Networks/CompSparSpiec/Order/", data_type, "_",
                                        compartment, "_taxonomy.txt"), 
                sep="\t", quote=FALSE, row.names=FALSE)
  }
}
```

# 07 Summaries taxonomy and edge tables at order level across all compartment

``` r
# Loop iterates over the data types (DNA/RNA data)
for(data_type in c("DNA", "RNA")){
  
  # Create open lists
  edges_tables <- list()
  taxonomy_tables <- list()

  # Loop iterates over the WWTP compartments 
  for(compartment in c("INF", "DNF", "NFC", "EFF")){
  
    # Load edges and taxonomy table of respective WWTP compartment and data set 
    taxonomy <- read.table(paste0("./Results/Networks/CompSparSpiec/Order/", data_type, "_", 
                                  compartment, "_taxonomy.txt"), header = T, sep = "\t", dec = ".")
    spar_spiec <- read.table(paste0("./Results/Networks/CompSparSpiec/Order/", data_type, "_", 
                                  compartment, ".txt"), header = T, sep = "\t", dec = ".")
    #colnames(spar_spiec) <- gsub("\\."," ",colnames(spar_spiec))
    
    # Change factors to charcters
    taxonomy <- data.frame(taxonomy %>% mutate_if(is.factor, as.character))
    spar_spiec <- data.frame(spar_spiec %>% mutate_if(is.factor, as.character))
    
    if(!all(c(spar_spiec$Node.1, spar_spiec$Node.2) %in% taxonomy$Key)){
      stop()
    }
    
    # Delete self loops
    spar_spiec <- spar_spiec[spar_spiec$Node.1 != spar_spiec$Node.2,]

    # Specify compartment in column names of edge and taxonomy table
    spar_spiec <- spar_spiec %>% rename(!!paste("Edge weights", data_type, compartment) := "Edge.weights",
                                        !!paste("Association", data_type, compartment) := "Association")
    taxonomy <- taxonomy %>% rename(!!paste("Number of reads", data_type, compartment) := "Number.of.reads")
    
    # Reduce the network, only keep nodes connected to parasites
    parasites <- taxonomy[taxonomy$Traits == "parasite","Key"]
    spar_spiec_parasites <- spar_spiec[unique(c(which(spar_spiec$Node.1 %in% parasites), 
                                                which(spar_spiec$Node.2 %in% parasites))),]
    taxonomy_parasites <- taxonomy[taxonomy$Key %in% 
                                     unique(c(spar_spiec_parasites$Node.1, spar_spiec_parasites$Node.2)),]
    
    # Save edge and taxonomy table in list
    edges_tables <- list.append(edges_tables, spar_spiec_parasites)
    taxonomy_tables <- list.append(taxonomy_tables, taxonomy_parasites)

  }
  
  # Merge edge and taxonomy table respectively, mutate factors to characters
  edges <- edges_tables %>% reduce(full_join)
  taxonomy <- taxonomy_tables %>% reduce(full_join)
  
  # Exclude WWTP compartments NFC & EFF from edge table
  # Afterward exclude rows with 4 NAs (means no entries for the respective association)
  edges_INF_DNF <- edges[,!grepl("NFC|EFF", colnames(edges))]
  edges_INF_DNF <- edges_INF_DNF[rowSums(is.na(edges_INF_DNF)) != 4,]
  taxonomy_INF_DNF <- taxonomy[,!grepl("NFC|EFF", colnames(taxonomy))]
  taxonomy_INF_DNF <- taxonomy_INF_DNF[taxonomy_INF_DNF$Key %in% 
                                         c(edges_INF_DNF$Node.1, edges_INF_DNF$Node.2),] 

  # Replace NAs with empty space, cytoscape has problems with NAs
  edges <- edges %>% replace(is.na(.), "")
  edges_INF_DNF <- edges_INF_DNF %>% replace(is.na(.), "")
  taxonomy <- taxonomy %>% replace(is.na(.), "")
  taxonomy_INF_DNF <- taxonomy_INF_DNF %>% replace(is.na(.), "")
  
  # Export tables
  colnames(edges) <- gsub("\\."," ",colnames(edges))
  write.table(edges, file = paste0("./Results/Networks/CompSparSpiec/Order/", data_type,
                                   "_all_parasites.txt"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  
  colnames(edges_INF_DNF) <- gsub("\\."," ",colnames(edges_INF_DNF))
  write.table(edges_INF_DNF, file = paste0("./Results/Networks/CompSparSpiec/Order/", data_type,
                                           "_INF_DNF_all_parasites.txt"), 
              sep="\t", quote=FALSE,row.names=FALSE)
  
  # Export taxonomy
  write.table(taxonomy, file = paste0("./Results/Networks/CompSparSpiec/Order/", data_type, 
                                      "_all_parasites_taxonomy.txt"), 
              sep="\t", quote=FALSE, row.names=FALSE)
  write.table(taxonomy_INF_DNF, file = paste0("./Results/Networks/CompSparSpiec/Order/", data_type, 
                                              "_INF_DNF_all_parasites_taxonomy.txt"), 
              sep="\t", quote=FALSE, row.names=FALSE)
}
```
