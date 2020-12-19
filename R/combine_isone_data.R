#' Combine the xls(x)'s downloaded from ISONE's website
#'
#' \code{combine_isone_data} combines the yearly excel file retrieved
#' from ISONE's website into one RDS file. In the process it standardizes
#' column names, and deals with DST start/end.
#'
#' **Note** the output file for this function is the input file
#' for \code{clean_isone_data}
#'
#' @param yaml_cfg a yaml file with configuration settings
#'
#' @export
#' @examples
#' \dontrun{
#' combine_isone_data('/home/wally/gefcom17/config/gen_and_dates_cfg.yaml')
#' }

combine_isone_data<-function(yaml_cfg = '') {

    if(yaml_cfg == '') {
        print('combine_load_data needs a config file specifying')
        print('the data dictionary and file names')
        return()
    }

    cfg = yaml::read_yaml(yaml_cfg)

    shtNames = c("ISONE CA", "ME", "NH", "VT", "CT", "RI", "SEMASS", "WCMASS",
                 "NEMASSBOST")
    zns = data.table()

    # Combine the individual years of data. Only the ISO tab has: Reg_Capacity_Price
    # (RegCP), # Reg_Service_Price (RSP) and  System_Load (SYSLoad).  Add 0 values
    # for zones.  For now, None of these are used to predict future demand.
    for (i in 2003:2020)
    {
        fname = sprintf('%s/%d_smd_hourly.xls', cfg$data_dir, i)
        if( i >= 2017 ) fname = sprintf('%s/%d_smd_hourly.xlsx', cfg$data_dir, i)

        print(fname)
        #Make tab and variable names common
        for( sht in shtNames)
        {
            if ( i < 2005 & sht == "ISONE CA") sht = "NEPOOL"

            if ( i > 2015 )
            {
                if( sht == "ISONE CA")   { sht = "ISO NE CA" }
                if( sht == "SEMASS")     { sht = "SEMA" }
                if( sht == "WCMASS")     { sht = "WCMA" }
                if( sht == "NEMASSBOST") { sht = "NEMA" }
            }

            temp = as.data.table(readxl::read_excel(fname, sht))
            setnames(temp, c(2,3,4), c("he","da_dem", "rt_dem"))
            temp[, he := as.integer(he)]
            temp[, Date := as.Date(Date)]

            if ( i > 2015 ) setnames(temp, c("Dry_Bulb","Dew_Point"),
                                           c("DryBulb","DewPnt"))

            if( substr(sht,1,3) == "ISO" | sht == "NEPOOL")
            {
                if ( i > 2015)
                {
                    setnames(temp, c("System_Load", "Reg_Service_Price",
                                     "Reg_Capacity_Price"), c("SYSLoad","RSP","RegCP"))
                } else
                {
                    # Reg_Service_Price was added in 2016
                    temp[, RSP := 0]
                }
                temp[, zone := "ISO"]
            } else
            {
                # Pad zone tables with empty data
                temp[, SYSLoad := 0]
                temp[, RSP := 0]
                temp[, RegCP := 0]
                temp[, zone := sht]
            }

            if( sht == 'ISO NE CA' & i > 2016) {
                temp[, `:=` (Min_5min_RSP = NULL, Max_5min_RSP = NULL,
                             Min_5min_RCP = NULL, Max_5min_RCP = NULL)]
            }
            zns = rbind(zns, temp)
        }
    }

    zns[ zone == "SEMASS", zone := "SEMA"]
    zns[ zone == "WCMASS", zone := "WCMA"]
    zns[ zone == "NEMASSBOST", zone := "NEMA"]

    #There doesn't appear to be bad data, but just in case
    zns = zns[ !is.na(Date),]
    zns[ , he := as.integer(he)]

    # Through 2015 the load files had a value of 0 for HE 2 on the DST start date
    # and 2x the load for HE 2 on the DST end date. Starting in 2016, HE 2 on
    # the DST start date was the average of HEs 1 and 3, and HE 2 on the DST
    # end date was the average of the two HE 2's
    znames = fread(sprintf('%s/%s', cfg$data_dir, cfg$zoneNamesFile))
    znames = znames[zone != 'TOTAL', zone] # TOTAL used for submissions only
    dstDate = fread(sprintf('%s/%s', cfg$data_dir, cfg$dstDatesFile))
    dst_strt = dstDate$st
    flds = c("rt_dem","da_dem","DA_LMP", "DA_EC", "DA_CC","DA_MLC","RT_LMP",
             "RT_EC", "RT_MLC")

    resIso = data.table()
    resZn = data.table()

    # DST start fix, ~20x faster than triple loop I used originally
    adt = zns[Date %in% dst_strt & (he ==1 | he == 3),]
    sdcol = names(adt)[names(adt) %notin% c('Date','zone')]
    ad2 = adt[, lapply(.SD, mean, na.rm=TRUE), by=c('Date', 'zone'),
                .SDcols = sdcol]
    zn2 = zns[ !(Date %in% dst_strt & he == 2),]
    zns = rbind(zn2, ad2)
    zns = zns[order(Date,he,zone)]

    # DST end fix
    dst_end = dstDate[year <= 2015, end] # see above
    zns[ Date %in% dst_end & he == 2, (flds) := lapply(.SD, function(x) x/2),
         .SDcols = flds]
    # output file for read is infile for clean
    saveRDS(zns, sprintf('%s/%s', cfg$data_dir, cfg$infile))
}
