# SAMESEA Activity 1.1 Explorer

A Shiny application developed within the SAMESEA project to visualize predicted presence and interaction risk between marine sentinel species and human activities across the Adriatic–Ionian (EUSAIR) region.

## Prerequisites

Before running the app, users have to download all the files in this repository and install all required R packages.

> [!IMPORTANT]
> Because GitHub cannot host very large files, one required dataset must be downloaded manually [here](https://owncloud.uniroma2.it/index.php/s/PSGimNpyB6SK4Wk).
> Once downloaded, the file **df_pred.rds** must be placed in the folder **data/**.

To install the R packages run this in your RStudio console:

```         
install.packages(c(
  "AMORE", "dplyr", "leaflet", "markdown", "randomForest",
  "shiny", "shinyjs", "sf", "viridis"
))
```

To run the SAMESEA ShinyApp open **app.R** and click on the "Run App" button on the top right of the main RStudio panel, or run:

```         
shiny::runApp()
```

> [!NOTE]
> A lightweight version of the app featuring only the Default explorer mode is available at [shinyapps.io](https://rlab.shinyapps.io/samesea-app-d/).

## Folder structure
.
├── app.R
├── data/    # input data (df_pred.rds must be added manually)
├── meta/    # default presence & interaction maps
├── models/  # pre-trained RF and MPN models
└── www/     # images and help files        
