Microbial community composition after quality filtering
================
Jule Freudenthal
2021-12-10

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
  add_header_above(c(" " = 1, "Metagenomics" = 4, "Metatranscriptomics" = 4))
```

<table class=" lightable-classic table table-condensed" style="font-family: Calibri; width: auto !important; margin-left: auto; margin-right: auto; font-size: 12px; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="empty-cells: hide;" colspan="1">
</th>
<th style="padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4">

<div style="border-bottom: 1px solid #111111; margin-bottom: -1px; ">

Metagenomics

</div>

</th>
<th style="padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4">

<div style="border-bottom: 1px solid #111111; margin-bottom: -1px; ">

Metatranscriptomics

</div>

</th>
</tr>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:center;">
Number of reads
</th>
<th style="text-align:center;">
\[%\]
</th>
<th style="text-align:center;">
Number of OTUs
</th>
<th style="text-align:center;">
\[%\]
</th>
<th style="text-align:center;">
Number of reads
</th>
<th style="text-align:center;">
\[%\]
</th>
<th style="text-align:center;">
Number of OTUs
</th>
<th style="text-align:center;">
\[%\]
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Prokaryotes
</td>
<td style="text-align:center;">
479487
</td>
<td style="text-align:center;">
94.34
</td>
<td style="text-align:center;">
1366
</td>
<td style="text-align:center;">
70.16
</td>
<td style="text-align:center;">
228569
</td>
<td style="text-align:center;">
42.49
</td>
<td style="text-align:center;">
1142
</td>
<td style="text-align:center;">
60.52
</td>
</tr>
<tr>
<td style="text-align:left;">
Protists
</td>
<td style="text-align:center;">
23198
</td>
<td style="text-align:center;">
4.56
</td>
<td style="text-align:center;">
388
</td>
<td style="text-align:center;">
19.93
</td>
<td style="text-align:center;">
294510
</td>
<td style="text-align:center;">
54.75
</td>
<td style="text-align:center;">
480
</td>
<td style="text-align:center;">
25.44
</td>
</tr>
<tr>
<td style="text-align:left;">
Fungi
</td>
<td style="text-align:center;">
3763
</td>
<td style="text-align:center;">
0.74
</td>
<td style="text-align:center;">
126
</td>
<td style="text-align:center;">
6.47
</td>
<td style="text-align:center;">
9980
</td>
<td style="text-align:center;">
1.86
</td>
<td style="text-align:center;">
196
</td>
<td style="text-align:center;">
10.39
</td>
</tr>
<tr>
<td style="text-align:left;">
Metazoa
</td>
<td style="text-align:center;">
1802
</td>
<td style="text-align:center;">
0.35
</td>
<td style="text-align:center;">
67
</td>
<td style="text-align:center;">
3.44
</td>
<td style="text-align:center;">
4876
</td>
<td style="text-align:center;">
0.91
</td>
<td style="text-align:center;">
69
</td>
<td style="text-align:center;">
3.66
</td>
</tr>
</tbody>
</table>
