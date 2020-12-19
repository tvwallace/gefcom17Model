#' Clean the ISONE data
#'
#' \code{clean_isone_data} removes variables that aren't needed, then adds some
#' derived weather variables, and potentially useful interactions. It then
#' writes an RDS file named and to a directory specified in the configuration
#' file.
#'
#' @param cln_yaml_cfg a yaml file with configuration settings
#'
#' @export
#' @examples
#' \dontrun{
#' clean_isone_data('/home/wally/gefcom17/config/gen_and_dates_cfg.yaml')
#' }

clean_isone_data<-function(cln_yaml_cfg = '') {

    if( cln_yaml_cfg == '' ) {
        print('cln_data_add_vars needs a yaml cfg file')
        return()
    }

    cfg = yaml::read_yaml(cln_yaml_cfg)

    if( cfg$infile == '' | cfg$outfile == '' | cfg$holsfile == ''|
        cfg$notneedfile == '') {
        print('cln_data_add_vars, need to provide:')
        print('input, output and holiday, and not needed file names')
        return()
    }

    zns = readRDS(sprintf('%s/%s', cfg$data_dir, cfg$infile))
    not_nd = fread(sprintf('%s/%s', cfg$data_dir, cfg$notneedfile))
    toKp = names(zns)[names(zns) %notin% not_nd$no_need]
    zns = subset(zns, select=toKp)
    zns = setnames(zns, "he", "hour")
    zones = unique(zns$zone)

    # some derived wx variables
    zns[, tc   := round((DryBulb - 32)*0.55555,1)]
    zns[, tdc  := round((DewPnt - 32)*0.55555,1)]
    zns[, rh   := 100*(exp((17.625*tdc)/(243.04+tdc))/exp((17.625*tc)/(243.04+tc)))]
    zns[, rh   := round(rh,2)]
    zns[, thi  := round(tc + (0.36*tdc) + 41.2,1)]
    zns[, cdd  := pmax(DryBulb-65,0)]
    zns[, cdd2 := round(pmax((0.4*(DryBulb+DewPnt)+15)-65,0),1)]
    zns[, DewPnt := NULL]
    zns[, hdd  := round(pmax(65-DryBulb,0),1)]

    # Add holidays
    hols = fread(sprintf('%s/%s', cfg$data_dir, cfg$holsfile))
    hols[, Date := as.Date(date)]
    hols[, isHol := 1]
    hols = subset(hols, select=c(Date,isHol,code))

    zns = merge(zns, hols, by="Date", all.x=T)
    zns[ is.na(isHol), isHol := 0]
    zns[ is.na(code), code := ""]

    # Add day of week, useful for prediction.  1 is Sunday
    zns[, dow   := wday(Date)]
    zns[, year  := year(Date)]
    zns[, month := month(Date)]
    zns[, day   := mday(Date)]

    #Add some potentially useful interactions, one example only
    zns[, dow  := as.factor(dow)]
    zns[, hour := as.factor(hour)]
    zns[, xhw := as.integer(interaction(hour,dow))]
    # ...

    zn2 = data.table()
    for (znName in zones)
    {
        zn1 = zns[zone == znName,]
        zn1 = add_year_frac(zn1)
        zn2 = rbind(zn2, zn1)
    }
    saveRDS(zn2, sprintf('%s/%s', cfg$data_dir, cfg$outfile))
}
