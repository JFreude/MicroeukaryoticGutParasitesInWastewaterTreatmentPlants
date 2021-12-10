Parasitic genera in WWTPs based on both metagenomic and
metatranscriptomic data
================
Jule Freudenthal
2021-12-10

**R version:** 3.6.2 (2019-12-12), Dark and Stormy Night **Packages**

-   kableExtra v. 1.3.4  
-   rlist v. 0.4.6.1

``` r
# Load packages
library(kableExtra)
library(rlist)
```

We assign functional traits to the taxa identified, using published
reference databases. Based to the these trait databases, we classified
the following taxa as **parasites**:

1.  **Protists**: All protist genera associated with human and animal
    gut and/or feces.  
2.  **Prokaryotes, fungi, microscopic metazoa**: All prokaryote, fungal
    and microscopic-metazoan genera that include potentially pathogenic
    species to humans and animals.

Now, we’ll give an overview of all identified parasitic genera.

``` r
# Load table with names of preprocessed count tables 
filenames_count_tables <- read.table("./R/Lists/Preprocessing_data.csv", header = TRUE, 
                                     row.names = 1, sep = ";", colClasses = "character")

# Create open list for tables
parasitic_taxa <- list()

# Loop iterates over the data types (DNA/RNA data)
for(data_type in c("DNA", "RNA")){
  
  # Load count table and metadata
  count_table <- read.table(paste0("./DataToAnalyse/PreprocessedData/", data_type, "_" , 
                                   filenames_count_tables[9,]), 
                            header = TRUE, sep = ";", dec = ".", row.names = 1)
  
  # Specify the microbial community
  count_table$`Microbial community` <- "Prokaryotes"
  count_table[count_table$Database == "PR2",]$`Microbial community` <- "Protists"
  count_table[count_table$Database == "PR2" & 
                count_table$Phylum == "Metazoa",]$`Microbial community` <- "Metazoa" 
  count_table[count_table$Database == "PR2" & 
                count_table$Phylum == "Fungi",]$`Microbial community` <- "Fungi"
  
  # Subset data set, keep only parasitic genera
  parasites <- count_table[count_table$Traits == "parasite", c("Genus", "Microbial community")]
  
  # Save the table in the list
  parasitic_taxa <- list.append(parasitic_taxa, parasites)
}

# Merge both tables to one
parasitic_taxa <- do.call(rbind, parasitic_taxa)

# Unique gut-associated genera
parasitic_taxa <- unique(parasitic_taxa)

# Order table
parasitic_taxa <- parasitic_taxa[order(parasitic_taxa$`Microbial community`, 
                                       parasitic_taxa$Genus),]

# Genus names italic
parasitic_taxa$Genus <- cell_spec(parasitic_taxa$Genus,italic = T)

# Subset table
parasitic_taxa <- cbind(parasitic_taxa[seq(nrow(parasitic_taxa)/2),
                                       c("Genus", "Microbial community")],
                        parasitic_taxa[((nrow(parasitic_taxa)/2)+1):nrow(parasitic_taxa),
                                       c("Genus", "Microbial community")])

# Delete row names
rownames(parasitic_taxa) <- NULL

# Visualize table
kable(parasitic_taxa, escape = F, booktabs = TRUE, align = "lclc") %>%
  kable_classic(full_width = F, html_font = "Calibri")  %>%
  kable_styling(font_size = 12, bootstrap_options = "condensed", position = "center", 
                latex_options = "scale_down") 
```

<table class=" lightable-classic table table-condensed" style="font-family: Calibri; width: auto !important; margin-left: auto; margin-right: auto; font-size: 12px; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Genus
</th>
<th style="text-align:center;">
Microbial community
</th>
<th style="text-align:left;">
Genus
</th>
<th style="text-align:center;">
Microbial community
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Aspergillus</span>
</td>
<td style="text-align:center;">
Fungi
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Dientamoeba</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Candida</span>
</td>
<td style="text-align:center;">
Fungi
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Entamoeba</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Cryptococcus</span>
</td>
<td style="text-align:center;">
Fungi
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Enterobryus</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Fusarium</span>
</td>
<td style="text-align:center;">
Fungi
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Giardia</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Penicillium</span>
</td>
<td style="text-align:center;">
Fungi
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Gregarines\_XX</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Trichosporon</span>
</td>
<td style="text-align:center;">
Fungi
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Guttulinopsis</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Acremonium</span>
</td>
<td style="text-align:center;">
Fungi
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Hexamita</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Cladophialophora</span>
</td>
<td style="text-align:center;">
Fungi
</td>
<td style="text-align:left;">
<span
style="  font-style: italic;   ">Hexamitinae-Enteromonadida\_X</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Mucor</span>
</td>
<td style="text-align:center;">
Fungi
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Ichthyophonus</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Ochroconis</span>
</td>
<td style="text-align:center;">
Fungi
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Monocystis</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Rhodotorula</span>
</td>
<td style="text-align:center;">
Fungi
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Perkinsida\_XXX</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Ascaris</span>
</td>
<td style="text-align:center;">
Metazoa
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Rhinosporidium</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Arcobacter</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Rhynosporidae\_X</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Bacillus</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Rosculus</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Campylobacter</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Sappinia</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Clostridium</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Sphaerothecum</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Corynebacterium</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Trepomonas</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Enterococcus</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Trimitus</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Helicobacter</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Blechomonas</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Klebsiella</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Crithidia</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Legionella</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Enteromonas</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Leptospira</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Helkesimastix</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Listeria</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Herpetomonas</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Mycobacterium</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Lacusteria</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Pseudomonas</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Leishmania</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Salmonella</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Leptomonas</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Vibrio</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Paratrypanosoma</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Yersinia</span>
</td>
<td style="text-align:center;">
Prokaryotes
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Phytomonas</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Acanthamoeba</span>
</td>
<td style="text-align:center;">
Protists
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Pseudotrichomonas</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Anurofeca</span>
</td>
<td style="text-align:center;">
Protists
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Sainouron</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Blastocystis</span>
</td>
<td style="text-align:center;">
Protists
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Tetratrichomonas</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Blastodinium</span>
</td>
<td style="text-align:center;">
Protists
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Trichomitus</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Copromyxa</span>
</td>
<td style="text-align:center;">
Protists
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Trichomonadidae\_X</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Creolimax</span>
</td>
<td style="text-align:center;">
Protists
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Trichomonas</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
<tr>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Cryptosporidium</span>
</td>
<td style="text-align:center;">
Protists
</td>
<td style="text-align:left;">
<span style="  font-style: italic;   ">Trypanosomatidae\_X</span>
</td>
<td style="text-align:center;">
Protists
</td>
</tr>
</tbody>
</table>
