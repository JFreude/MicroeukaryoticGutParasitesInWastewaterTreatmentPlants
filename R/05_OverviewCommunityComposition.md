Microbial community composition after quality filtering
================
Jule Freudenthal
2022-02-03

**R version:** 3.6.2 (2019-12-12), Dark and Stormy Night  
**Packages**

-   dplyr v. 1.0.7  
-   kableExtra v. 1.3.4  
-   rlist v. 0.4.6.1

``` r
# Load packages
library(dplyr)
library(kableExtra)
library(rlist)
```

We’ll give an overview of the number of reads and OTUs of prokaryotes,
protists, fungi and microscopic metazoans.

``` r
# Load table with names of preprocessed count tables
filenames_count_tables <- read.table("./R/Lists/Preprocessing_data.csv", header = TRUE, 
                                     row.names = 1, sep = ";", colClasses = "character")

# Create open list for tables
overview <- list()

# Loop iterates over the data types (DNA/RNA data)
for(data_type in c("DNA", "RNA")){
  
  # Load count table and metadata
  count_table <- read.table(paste0("./DataToAnalyse/PreprocessedData/", data_type, "_" , 
                                   filenames_count_tables[8,]), 
                            header = TRUE, sep = ";", dec = ".", row.names = 1)
  
  # Specify the microbial community
  count_table$Microbial_Community <- "Prokaryotes"
  count_table[count_table$Database == "PR2",]$Microbial_Community <- "Protists"
  count_table[count_table$Database == "PR2" & count_table$Phylum == "Metazoa",]$Microbial_Community <- 
    "Metazoa" 
  count_table[count_table$Database == "PR2" & count_table$Phylum == "Fungi",]$Microbial_Community <- "Fungi"
  
  # Extract count data 
  counts <- select_if(count_table, is.numeric)
  
  # Number of reads
  number_of_reads <- c(Prokaryotes=sum(counts[count_table$Microbial_Community == "Prokaryotes",]),
                       Protists=sum(counts[count_table$Microbial_Community == "Protists",]),
                       Fungi=sum(counts[count_table$Microbial_Community == "Fungi",]),
                       Metazoa=sum(counts[count_table$Microbial_Community == "Metazoa",]))
  percentage_reads <- round(number_of_reads/sum(number_of_reads)*100,2)
  
  # Check if selecting was correct
  if(sum(number_of_reads) != sum(counts)){
    stop("Number of reads must match summed counts")}
  
  # Number of OTUs
  number_of_OTUs <- c(Prokaryotes=nrow(counts[count_table$Microbial_Community == "Prokaryotes",]),
                       Protists=nrow(counts[count_table$Microbial_Community == "Protists",]),
                       Fungi=nrow(counts[count_table$Microbial_Community == "Fungi",]),
                       Metazoa=nrow(counts[count_table$Microbial_Community == "Metazoa",]))
  percentage_OTUs <- round(number_of_OTUs/sum(number_of_OTUs)*100,2)
  
  # Check if selecting was correct
  if(sum(number_of_OTUs) != nrow(counts)){
    stop("Number of OTUs must match number of rows of counts data")}
  
  # Merge to data frame
  table <- setNames(data.frame(cbind(number_of_reads, percentage_reads, number_of_OTUs, 
                                     percentage_OTUs)), 
                    c("Number of reads", "[%]", "Number of OTUs", "[%]"))
  # Save the table in the list
  overview <- list.append(overview, table)
  
}

# Merge both tables to one
table <- do.call(cbind, overview)

# Visualize table
kable(table, booktabs = TRUE, align = "cccccccc") %>%
  kable_classic(full_width = F, html_font = "Calibri") %>%
  kable_styling(font_size = 12, bootstrap_options = "condensed", position = "center", 
                latex_options = "scale_down") %>%
  add_header_above(c(" " = 1, "Metagenomics" = 4, "Metatranscriptomics" = 4)) %>%
  as_image(width = 10) 
```

![](C:\Users\Jule\AppData\Local\Temp\Rtmp2peNa0\file5ff0824f79.png)<!-- -->
