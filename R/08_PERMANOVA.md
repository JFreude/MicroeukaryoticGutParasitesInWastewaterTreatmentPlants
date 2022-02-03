Permutational Multivariate Analysis of Variance of WWTP compartments
================
Jule Freudenthal
2022-02-03

**R version:** 3.6.2 (2019-12-12), Dark and Stormy Night  
**Packages**

-   vegan v. 2.5-7

``` r
# Load packages
library(vegan)
```

We test for differences between the treatment compartments both in terms
of community composition (rDNA) and activity (rRNA), using Permutational
Multivariate Analysis of Variance (PERMANOVA, function adonis, package
vegan).

``` r
# Loop iterates over the data types (DNA/RNA data) 
for(data_type in c("DNA", "RNA")){
  
  # Load count- and meta- and taxonomy data
  load(paste0("./DataToAnalyse/RData/", data_type, "_data.rds"))

  # Normalize counts (relative counts)
  normalized_counts <- sweep(x = counts, MARGIN = 2, STATS = colSums(counts), FUN = '/')
  
  # Calculate PERMANOVA 
  permanova <- adonis(t(normalized_counts) ~ Sample_Place, data = metadata, method ="bray",
                      strata = metadata$Location, permutations=999)
  perm<-permanova$aov.tab
  
  print(data_type)
  print(perm)
}
```

    ## [1] "DNA"
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##              Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Sample_Place  3    3.9931 1.33103  16.484 0.59258  0.001 ***
    ## Residuals    34    2.7454 0.08075         0.40742           
    ## Total        37    6.7385                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## [1] "RNA"
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##              Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
    ## Sample_Place  3    3.3391 1.11304  6.3069 0.3509  0.001 ***
    ## Residuals    35    6.1768 0.17648         0.6491           
    ## Total        38    9.5160                 1.0000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
