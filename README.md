# Baltic Sea metagenome dashboard
Carlo Berg and Anders Andersson

An R Shiny App dashboard to display metagenomic data from various samples of the Baltic Sea. An external sample can be uploaded and used to predict environmental parameters from its metagenomic gene abundance data and application of a random forest model trained with the Baltic Sea metagenomic data.

![Settings tab](img/settings.png)  
  
![Heatmap tab](img/heatmap.png)

## Usage

From the panel on the left side of the app several actions can be performed.  First, the samples can be selected, filtered, and an external data file uploaded that is used alongside with the Baltic Sea metagenome samples. The external data file should be in tab-delimited format. Then, you can view the metagenomic data in a heatmap and view or download the count data as a table. In the final step, the external data file can be used to predict values for several environmental parameters by application of a random forest model that was trained with the Baltic Sea metagenomic data. 

## Live version
A deployed live version is also available at [shinyapps.io](https://cberg.shinyapps.io/baltic-sea-metagenome-dashboard/)



The [BONUS BLUEPRINT](https://blueprint-project.org) project has received funding from BONUS (Art 185), funded jointly by the EU and the national funding institutions of Denmark, Sweden, Germany, Finland, and Estonia.