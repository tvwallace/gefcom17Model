% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod_load.R
\name{model_load}
\alias{model_load}
\title{Produce estimates of quantiles of electricity load for ISO New England (ISONE)}
\usage{
model_load(load_cfg = "", dts_cfg = "")
}
\arguments{
\item{load_cfg}{a yaml file containing modeling parameters}

\item{dts_cfg}{a yaml file containing general setup parameters}
}
\description{
\code{model_load} produces quantile estimates of electricity load.
This function is driven by yaml configuration files,
examples of which are included with this package.  The configuration
files also explain the meaning/import/use of the parameters.
}
\details{
If writeFiles = 1 (in the modeling configuration file), then predictions
and results are written to the directories specified in the configuration
file.
}
\examples{
\dontrun{
model_load('~/gefcom17/config/load_cfg.yaml',
'~/gefcom17/config/gen_and_dates_cfg.yaml')
}
}
