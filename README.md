# fragilitycurves R package
- This currently includes functions used for modelling fragility curves, specifically using ordinal models and accounting for spatial correlation in damage. Functions to estimate annual losses based on these are also included.
- A tutorial using simulated earthquake data in Haiti is provided ("demo/damage_spatial_corr_demo.Rmd").
- A similar analysis conducted on real-life damage data from the 2010 Haiti earthquake given in "demo/haiti_data_demo.pdf".
- This is work in progress: please inform Michele of any issues/suggestions and if you have any functions to add in.

Instructions for installing the package:
1. Download the zip (for Windows) or tar.gz (for Mac) file from GitHub.
2. If you do not have compiler tools already installed, follow [these instructions](https://github.com/kaskr/adcomp/wiki/Download) to install them. 
3. In RStudio, install the non-standard package dependencies using:

```
install.packages(c('TMB', 'RcppEigen', 'raster', 'ggplot2', 'geoR', 'gridExtra', 'rgdal', 'sp', 'MASS', 'rgeos'))
```

4. In RStudio, use "Tools" -> "Install Packages" ->  Install from: "Package Archive File" and find where the downloaded zip or tar.gz is.
