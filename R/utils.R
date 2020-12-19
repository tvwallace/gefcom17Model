`%notin%` <- function(x,y) !(x %in% y)

#' Get SMD data from ISO New England (ISONE)'s website
#'
#' Note: due to lack of file-naming consistency, this function only
#' works for years 2016 and later.  See the isone website referenced
#' in the code for more information.
#'
#' \code{get_isone_smd_data} retrieves hourly load and price data from
#' ISONE's website, and writes the yearly files to the directory specified
#'
#' @param minYr an integer, the first year for which you want data, defaults to 2011
#' @param maxYr an integer, the last year for which you want data, defaults to current year
#' @param data_dir a string, the directory in which to write files,
#' without the trailing '/'
#'
#' @export
#' @examples
#' \dontrun{
#' get_isone_smd_data(2016, 2018, '/home/frodo/isone/smd_data')
#' }

get_isone_smd_data<-function(minYr = 2016, maxYr = year(Sys.Date()), data_dir = '')
{
    if(data_dir == '') {
        print('get_isone_smd_data: need to specify directory for data files')
        return()
    }

    year_range = seq(minYr, maxYr)

    for(yr in year_range ) {

        url1 = 'https://www.iso-ne.com/static-assets/documents/'
        url2 = sprintf('%s/02/', yr)

        excStub = 'xlsx'
        yrStr = sprintf('%d_', yr)

        if( yr < 2017 ) {
            excStub = 'xls'
            yrStr = ''
        }

        infileName = sprintf('%ssmd_hourly.%s', yrStr, excStub)
        outfileName = sprintf('%d_smd_hourly.%s', yr, excStub)
        full_url = paste0(url1, url2, infileName)
        download.file(full_url, destfile = sprintf('%s/%s',  data_dir,
                                                   outfileName))
    }
}

#' Calculate the pinball loss between actual and the predicted quantile
#'
#' \code{pinball} calculates the goodness of fit between an actual quantity,
#' e.g. electricity demand, and a predicted quantile of that quantity.
#'
#' @param y_act a single number or numerical vector of actual quantities
#' @param y_pred  the predicted quantile of that quantity
#' @param qu the target quantile \eqn{0 < qu < 1}
#' @param dosum when 1, return the pinball loss as a single number
#' @return a vector of pinball values, or their sum, depending on value of dosum
#' @export
#' @examples
#' \dontrun{
#' act = rep(2700, 3)
#' pred = rep(2671, 3)
#' ploss = pinball( act, pred, 0.1)
#' # should return 8.7
#' }

pinball = function(y_act, y_pred, qu, dosum = 1)
{
    if( qu <= 0 | qu > 1) {
        print('pinball needs a proper quantile: 0 < tau < 1')
        return()
    }
    res = qu * (y_act - y_pred) * ((y_act - y_pred) >= 0) +
          (qu - 1) * (y_act - y_pred) * ((y_act - y_pred) < 0)

    if( dosum == 1 ) { res = sum(res) }
    res
}

#' Calculate the average pinball across multiple quantiles
#'
#' \code{avgPinball} calls \code{pinball} across a range of
#' quantiles and returns the average
#'
#' @param qdt a data table containing the actual results
#' and predicted quantiles.
#'
#' **Note:** \code{avgPinball} expects the names in qdt to be 'act'
#' for the actuals, and of the form 'q0.1' for the predicted quantiles
#'
#' @param qdt  a data table containing actuals and predicted quantiles
#' @param quants the target quantiles \eqn{0 < quantile < 1} to assess
#' @return the average pinball score across the samples and quantiles
#' @export
#' @examples
#' \dontrun{
#' qdt = data.table(act=rep(2700,3),
#'                  q0.1 = c(2600,2500,2500),
#'                  q0.5=c(2900,2800,2700),
#'                  q0.9 = c(3100, 3000, 2900))
#' qnts = c(0.1, 0.5, 0.9)
#' avgPinball( qdt, qnts) # should return 32.22
#' }

avgPinball = function(qdt, quants)
{
    y = qdt[, act]
    avgPB = 0

    for (q in quants)
    {
        tname = paste("q", q ,sep="")
        f = qdt[[tname]]
        avgPB = avgPB + pinball(y, f, q)
    }
    round(avgPB/(dim(qdt)[1] * length(quants)),4)
}

#' Convert factor variables into integers
#'
#' Because xgboost, at least the way I used it, does not work with non-numerical
#' data \code{fact_to_int} converts factors to integers.  This function was
#' more important when I thought inclusion of interactions, e.g.
#' hour x weekend/weekday, would help performance.  It turned out not to be the
#' case.
#'
#' @param dt a data table possibly containing factors variables
#' @return the data table with any factors converted in integers

fact_to_int<-function(dt) {

    if( dim(dt)[1] < 1) {
        print('fact_to_int: data table needs something in it')
        return()
    }

    featNames = names(dt)

    for (ftre in featNames)
    {
        if (class(dt[[ftre]])=="factor") dt[, (ftre) :=
                                              as.integer( dt[, get(ftre)])]
    }
    return(dt)
}

#' Add 'year fraction' (yrFrac) to a data table
#'
#' So-called 'year fraction' is the percentage of hours that have elapsed
#' given the date and hour. For example at 2 am on 3/1/17, yrFrac = 0.16187:
#' (59 days * 24 hours + 2 hours)/8760 hours in a year.
#'
#' @param dt a data table to which yrFrac will be added
#' @return the data table with yrFrac added

add_year_frac<-function(dt) {

    if( dim(dt)[1] < 1) {
        print('add_year_frac: data table needs something in it')
        return()
    }

    nhr = 8760
    uyrs = unique(dt[,year])
    dt[ ,yrFrac := 0]

    for ( yr in uyrs)
    {
        divvy = nhr
        if( yr %% 4 == 0) { divvy = 8784 }
        nrow = dim(dt[year == yr,])[1]
        adder = 0

        # ISONE data starts on 3/2/2003
        if( yr == 2003) { adder = 60*24 }
        dt[ year == yr, yrFrac := round(((1:nrow)+adder)/divvy,7)]
    }
    return(dt)
}
