# fragilitycurves R package
- This currently includes functions used for modelling fragility curves, specifically using ordinal models and accounting for spatial correlation in damage. Functions to estimate annual losses based on these are also included.
- A tutorial using simulated earthquake data in Haiti is provided ("damage_spatial_corr_demo_user.Rmd"/"damage_spatial_corr_demo.pdf").
- This is work in progress: please inform Michele of any issues/suggestions and if you have any functions to add in.

Instructions for installing the package:
1. Download the tar.gz file from GitHub.
2. If you do not have compiler tools already installed, follow [these instructions](https://github.com/kaskr/adcomp/wiki/Download) to install them. 
3. In RStudio, install the non-standard package dependencies using:

```
install.packages(c('TMB', 'gridExtra', 'raster', 'geoR', 'lemon', 'gstat', 'rgeos', 'rgdal', 'RcppEigen'))
```

4. In RStudio, use "Tools" -> "Install Packages" ->  Install from: "Package Archive File" and find where the downloaded tar.gz is.
