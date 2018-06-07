---
title: "heatmap"
author: "Carlo Berg"
date: '2018-06-04'
output: html_document
---




```r
eggNOG.tpm <- read.delim(file = "../data/lmo2012_transect2014_redox2014.eggNOG.tpm.tsv")
```



```r
KEGG.tpm <- read.delim(file = "../data/lmo2012_transect2014_redox2014.KEGG-pathway-module.tpm.tsv")


rownames(KEGG.tpm) <- KEGG.tpm[,1]
KEGG.tpm <- KEGG.tpm[,-1]
KEGG.tpm <- apply(KEGG.tpm, 2, as.numeric)
KEGG.tpm <- as.matrix(KEGG.tpm)

d3heatmap(KEGG.tpm) #, colors = cols, scale = "column"
```

preservebc363deae7263bd8

