
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gefcom17Model

<!-- badges: start -->

<!-- badges: end -->

This package contains the code I used to compete in GEFCom2017, a
hierarchical probabilistic load forecasting contest. More details
concerning the contest may be found at:

<https://en.wikipedia.org/wiki/GEFCom#GEFCom_2017>

<http://www.drhongtao.com/gefcom/2017>

Competing under the nom de plume ‘ch46’, I finished 10th out of about 70
academic and company teams enrolled.

The requirement of this contest was to produce estimated quantiles of
electricity demand for the eight load zones in ISO New England (ISONE),
as well as two aggregated zones: Massachusetts, and the entire ISO.

Load and weather data was from ISONE, although one had the option to
supplement it. I chose not to, competing in the defined track.

The evaluation metric was the so-called pinball loss, averaged over all
quantiles. It is defined below, where
<img src=https://private.codecogs.com/png.download?%5Ctau> is the target
quantile:

<img src =https://private.codecogs.com/png.download?%5Cbegin%7Bmultiline%7D%20%5Cbegin%7Baligned%7D%20F_%7B%5Ctau%7D%28y%2C%5Chat%7By%7D%29%20%3D%20%26%5Chspace%7B.1cm%7D%28y%20-%20%5Chat%7By%7D%29%20%5Ctau%20%5Chspace%7B1cm%7D%20%5Ctextrm%7B%20if%20%7D%20y%20%5Cgeq%20%5Chat%7By%7D%5C%5C%5C%20%3D%20%26%5Chspace%7B.1cm%7D%28%5Chat%7By%7D%20-%20y%29%20%281%20-%20%5Ctau%29%20%5Ctextrm%7B%20if%20%7D%20%5Chat%7By%7D%20%3E%20y%20%5Cend%7Baligned%7D%20%5Cend%7Bmultiline%7D>

Accurate temperature estimates are a requirement for accurate load
estimates. I estimated temperature quantiles by estimating mean and
variance of temperature, where variance was the squared residuals of the
mean prediction:

<img src=https://private.codecogs.com/png.download?%5Csigma_%7Btemp%7D%5E%7B2%7D%20%3D%20%28temp%20-%20%5Chat%7Btemp%7D%29%5E%7B2%7D>

From those predictions I estimated quantiles of temperature, e.g for the
90th quantile:

<img src=https://private.codecogs.com/png.download?q90_%7Btemp%7D%20%3D%20%5Chat%7Bmu_%7Btemp%7D%7D%20+%20qnorm%280.90%29*%5Chat%7B%5Csigma_%7Btemp%7D%7D>

This methodology assumes that conditional mean and conditional median
are the same. I did not check this assumption. It may be true for
temperature; it almost certainly isn’t for load.

I modeled temperature using xgboost and gamboost. For a still unknown
reason I got slightly better validation results using gamboost.

To model load I tried a few different algorithms; please see the source
code for details. What worked best was to model mean load with xgboost,
then feed that model the estimated temperature quantiles to generate
load quantiles. That model is specified as ‘gt’ in the source code.

I spent many hours trying to find explanatory variables, mostly
interaction terms, that would enhance performance. In the end both the
temperature and load models contained a modest number of variables. That
specification may be found in the yaml configuration files in extdata.

## Installation

``` r
devtools::install_github('gefcom17model')
```

## Example

``` r
library(gefcom17Model)

# Assuming you have a directory structure like mine
gen_cfg = '~/gefcom17/config/gen_and_dates_cfg.yaml'
wx_cfg = '~/gefcom17/config/wx_cfg.yaml'
ld_cfg = '~/gefcom17/config/load_cfg.yaml'

# one time
get_isone_smd_data(2017, 2020, '~/gefcom17/smd_data')
combine_isone_data(gen_cfg)
clean_isone_data(gen_cfg)

# lather, rinse, repeat
model_temperature(wx_cfg, gen_cfg)
model_load(ld_cfg, gen_cfg)
```