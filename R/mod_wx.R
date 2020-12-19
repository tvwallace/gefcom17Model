#' Produce estimates of quantiles of temperature
#'
#' \code{model_temperature} produces quantile estimates of
#' temperature.  This function is driven by yaml configuration files,
#' examples of which are included with this package.  The configuration
#' files also explain the meaning/import/use of the parameters.
#'
#' **Note:** The gamboost specifications are largely hard-coded.  I could
#' not figure out how to parameterize the function call with respect to
#' the formula: gamboost( tc ~ bbs(...) + ...)
#'
#' If writeFiles = 1 (in the modeling configuration file), then predictions
#' and results are written to the directories specified in the configuration
#' file.
#'
#' @param  wx_cfg  a yaml file containing modeling parameters
#' @param  dts_cfg a yaml file containing general setup parameters
#' @export
#' @examples
#' \dontrun{
#' model_weather('~/gefcom17/config/wx_cfg.yaml',
#' '~/gefcom17/config/gen_and_dates_cfg.yaml')
#' }

model_temperature<-function(wx_cfg  = '', dts_cfg = '') {

    if( wx_cfg == '' | dts_cfg == '') {
        print('model_weather needs 2 config files:')
        print('1) model parameters')
        print('2) general parameters and dates')
        return()
    }

    cfg = yaml::read_yaml(wx_cfg)
    dcfg = yaml::read_yaml(dts_cfg)

    doInt  = cfg$doInt
    doXGB  = cfg$doXGB
    doImp  = cfg$doImp
    doGam  = cfg$doGam
    doTest = cfg$doTest

    quants = (1:(cfg$num_quants-1))/cfg$num_quants

    if (cfg$rnd_seed > 0) { set.seed(cfg$rnd_seed) }

    ozns = readRDS(sprintf('%s/%s', dcfg$data_dir, dcfg$outfile))

    if (doInt == 1) ozns = fact_to_int(ozns)

    resDT = data.table()

    zones = unique(ozns[,zone])
    if(cfg$zones[1] != 'all') zones = cfg$zones

    for( tzne in zones)
    {
        print(sprintf("Zone: %s*****", tzne), quote=F)
        zns = ozns[zone == tzne,]
        zns = zns[order(Date,hour)]
        shftDays = 0

        valdt = data.table( dst = dcfg$val_dts$strt, dend = dcfg$val_dts$end)
        lstdt = data.table( eod = dcfg$val_dts$lst)

        mths = data.table(mbeg = dcfg$mth_data_beg, mend = dcfg$mth_data_end)

        print(mths)

        if( dim(mths)[1] == 1) bmme = paste(mths$mbeg, mths$mend, sep='')

        if (doTest == 1)
        {
            valdt = data.table( dst = dcfg$tst_dts$strt, ddend = cfg$tst_dts$end)
            lstdt = data.table( eod = dcfg$tst_dts$lst)
        }

        nmth = dim(mths)[1]
        nrun = dim(valdt)[1]

        for ( k in  1:nmth)
        {
            for ( i in 1:nrun)
            {
                trng = valg = gam = gamV = gpred = gredv = qdt = qdtg = NULL
                mbeg = mths[k, mbeg]
                mend  = mths[k,mend]
                eDate = lstdt[i,eod]

                fstDate = as.Date(eDate) - 365*dcfg$n_yrs_data
                zns = zns[ Date >= fstDate,]

                valg = zns[ Date > valdt[i, dst] & Date < valdt[i, dend],]

                trun = sprintf("RUN: %s - %d",tzne, unique(valg$year))
                print(trun, quote=F)

                vlmn = valg[1,]$month
                vlyr = valg[1,]$year

                # Get correct data range, e.g. Oct '15 - Feb '16
                if ( mbeg == mend)     { trng = zns[ month == mbeg & Date < eDate,] }
                else if ( mbeg > mend) { trng = zns[ (month >= mbeg | month <= mend)
                                                      & Date < eDate,] }
                else                    { trng = zns[ (month >= mbeg & month <= mend)
                                                      & Date < eDate,] }

                if( doXGB == 1)
                {
                    valLab = valg[[cfg$wx_mn_tgt]] #needs to be vector
                    valt = subset(valg, select = cfg$wx_mn_pred)
                    dval<-Matrix::sparse.model.matrix(~.-1, data = valt)
                    dvald<-xgboost::xgb.DMatrix(data=dval, label=valLab)

                    trnLab = trng[[cfg$wx_mn_tgt]]
                    trnt = subset(trng, select=cfg$wx_mn_pred)
                    dtrn<-Matrix::sparse.model.matrix(~.-1, data = trnt)
                    dtrnd<-xgboost::xgb.DMatrix(data=dtrn, label=trnLab)

                    watchlist<-list(val=dvald, train=dtrnd)

                    nrndm = switch(tzne,
                                  "NH"   = cfg$nrnd_nh,
                                  "ME"   = cfg$nrnd_me,
                                  "VT"   = cfg$nrnd_vt,
                                  "CT"   = cfg$nrnd_ct,
                                  "RI"   = cfg$nrnd_ri,
                                  "SEMA" = cfg$nrnd_sema,
                                  "WCMA" = cfg$nrnd_wcma,
                                  "NEMA" = cfg$nrnd_nema,
                                  "ISO"  = cfg$nrnd_iso,
                                  cfg$nrnd_def)

                    # For test, have more data
                    if( doTest == 1) nrndm = round(cfg$tstRndMult*nrndm,0)

                    paramm <- list( objective  = cfg$obj, booster = cfg$bster,
                                    eta = cfg$eta, max_depth = cfg$mxdpth,
                                    subsample = cfg$subsamp,
                                    colsample_bytree = cfg$colsamp,
                                    min_child_weight = cfg$minchwt)

                    mod <- xgboost::xgb.train( params = paramm, data = dtrnd,
                                               nrounds = nrndm, verbose = cfg$verb,
                                               watchlist = watchlist,
                                               print_every_n = cfg$prnt_n,
                                               maximize = FALSE)

                    vpred = predict(mod, dvald)

                    rmse  = sqrt(sum((vpred-valLab)^2)/length(valLab))

                    if (doImp == 1) {
                        nimp = xgboost::xgb.importance(dimnames(dtrn)[[2]],
                                                       model=mod)
                    }

                    # Variance training
                    print("Var regr", quote=F)
                    trnt2 = subset(trng, select=cfg$wx_var_pred)
                    valt2 = subset(valg, select=cfg$wx_var_pred)

                    dval<-Matrix::sparse.model.matrix(~.-1, data = valt2)
                    dtrn<-Matrix::sparse.model.matrix(~.-1, data = trnt2)

                    trnLab2 = (trnLab - predict(mod, dtrnd))^2
                    valLab2 = (valLab - vpred)^2

                    dvald<-xgboost::xgb.DMatrix(data = dval,label = valLab2)
                    dtrnd<-xgboost::xgb.DMatrix(data = dtrn,label = trnLab2)
                    watchlist<-list(val = dvald,train = dtrnd)

                    nrndv = switch(tzne,
                                   "NH"   = cfg$nrnd_nhv,
                                   "ME"   = cfg$nrnd_mev,
                                   "VT"   = cfg$nrnd_vtv,
                                   "CT"   = cfg$nrnd_ctv,
                                   "RI"   = cfg$nrnd_riv,
                                   "SEMA" = cfg$nrnd_semav,
                                   "WCMA" = cfg$nrnd_wcmav,
                                   "NEMA" = cfg$nrnd_nemav,
                                   "ISO"  = cfg$nrnd_isov,
                                   cfg$nrnd_defv)

                    if( doTest == 1) nrndv = round(cfg$tstRndMult*nrndv,0)

                    paramv <- list( objective  = cfg$obj, booster = cfg$bster,
                                    eta = cfg$eta_v, max_depth = cfg$mxdpthv,
                                    subsample = cfg$subsampv,
                                    colsample_bytree = cfg$colsampv,
                                    min_child_weight = cfg$minchwtv)

                    modV <- xgboost::xgb.train( params = paramv, data = dtrnd,
                                                nrounds = nrndv, verbose = 1,
                                                watchlist = watchlist,
                                                print_every_n = cfg$prnt_n,
                                                maximize = FALSE)

                    predVar = predict(modV, dvald)
                    rmsev  = sqrt(sum((predVar-valLab2)^2)/length(valLab2))

                    if (doImp == 1) nimpv = xgboost::xgb.importance(model=modV)

                    # Hack quantiles, assumes conditional mean == conditional median
                    # Also, use of qnorm probably not so correct
                    sds = sqrt(predVar)
                    qdt = valg
                    qdt = setnames(qdt, "tc", "act")

                    for (quant in quants)
                    {
                        qnm = paste("q",quant,sep="")
                        temp = vpred + qnorm(quant)*sds
                        qdt[[qnm]] = round(temp,1)
                    }

                    apb = avgPinball(qdt, quants)
                    qdt[ , Date := as.character(Date)]

                    if( cfg$writeFiles) {
                        fname = sprintf('%s/%sX-%d-%d-%s.csv', dcfg$temp_dir,
                                        tzne,vlmn,vlyr, bmme)
                        fwrite(qdt, fname)
                    }
                } #doXGB

                if ( doGam == 1)
                {
                    print('GAM boost', quote=F)
                    gam = mboost::gamboost(tc~ bbs(yrFrac,
                                                   degree = cfg$gb_mn_yf_dg,
                                                   df = cfg$gb_mn_yf_df,
                                                   knots = cfg$gb_mn_yf_kt,
                                                   differences = cfg$gb_mn_yf_dlt) +
                                                bbs(hour, degree = cfg$gb_mn_hr_dg,
                                                    df = cfg$gb_mn_hr_df,
                                                    knots = cfg$gb_mn_hr_kt),
                                            data = trng, family = Gaussian(),
                                            control = boost_control (mstop = cfg$gb_mn_stp,
                                                                nu = cfg$gb_mn_nu,
                                                                risk = cfg$gb_mn_rsk,
                                                                trace = cfg$gb_mn_trc,
                                                                stopintern = cfg$gb_mn_stpi,
                                                                center = cfg$gb_mn_cntr))

                    gpred = predict(gam, newdata=as.data.frame(valg))
                    rmse2 = sqrt(sum((valg[,act]-gpred)^2)/dim(valg)[1])

                    trng[, resi := (trng[, tc] - gam$fitted())^2]
                    valg[, resi := (valg[, act] -gpred)^2]

                    gamv = gamboost(resi ~ bbs(yrFrac, degree = cfg$gb_sd_yf_dg,
                                               df = cfg$gb_sd_yf_df,
                                               knots = cfg$gb_sd_yf_kt,
                                               differences = cfg$gb_sd_yf_dlt) +
                                            bbs(hour),
                                    data = trng, family= Gaussian(),
                                    control=boost_control(mstop = cfg$gb_sd_stp,
                                                          nu = cfg$gb_sd_nu,
                                                          risk = cfg$gb_sd_rsk,
                                                          trace = cfg$gb_sd_trc,
                                                          stopintern = cfg$gb_sd_stpi,
                                                          center = cfg$gb_sd_cntr))

                    gpredv = predict(gamv, newdata=as.data.frame(valg))
                    rmse2v = sqrt(sum((valg[,resi]-gpredv)^2)/dim(valg)[1])
                    sdsg = sqrt(gpredv)
                    qdtg = subset(valg, select=c(Date,hour,month,day,year,act))

                    for (quant in quants)
                    {
                        qnm = paste("q",quant,sep="")
                        temp = gpred + qnorm(quant)*sdsg
                        qdtg[[qnm]] = round(temp,1)
                    }

                    apb2 = avgPinball(qdtg, quants)
                    qdtg[ , Date := as.character(Date)]

                    if( cfg$writeFiles) {
                        fname = sprintf('%s/%sG-%d-%d-%s.csv', dcfg$temp_dir,
                                        tzne, vlmn, vlyr, bmme)
                        fwrite(qdtg, fname)
                    }
                } #doGam

                tmpdt = data.table(zone=tzne,valMth = vlmn, valYr = vlyr, monBeg = mbeg,
                                   monEnd = mend, endDate = eDate, rmse_xgb_mn = rmse,
                                   rmse_gam_mn = rmse2, rmse_xgb_var = rmsev,
                                   rmse_gam_var = rmse2v, pinball_xgb = apb,
                                   pinball_gam = apb2)
                resDT = rbind(resDT, tmpdt)
            } # val dates

            # avg of all validation periods
            bmval = paste(mbeg, mend, sep="")
            tmp2 = resDT[zone == tzne & bmme == bmval,]
            tmpdt = tmp2[, lapply(.SD, mean), .SDcols = names(tmp2[, 7:12])]
            tmpdt[, `:=` (zone=tzne, valMth = 999, valYr = 999, monBeg = mbeg,
                          monEnd = mend, endDate = eDate)]
            resDT = rbind(resDT, tmpdt)
        } # month pairs
    } # zones for loop
    print(resDT)

    if(cfg$writeFiles) {
        st1 = month.abb[vlmn]
        resfn = sprintf('%s/res_temps_%s.csv', dcfg$res_dir, st1)
        fwrite(resDT, resfn)

        if(doImp) {
            nmfn = sprintf('%s/nimp_mean_temp_%s.csv', dcfg$res_dir, st1)
            fwrite(nimp, nmfn)
            nvfn = sprintf('%s/nimp_var_temp_%s.csv',  dcfg$res_dir, st1)
            fwrite(nimpv, nvfn)
        }
    }
}
