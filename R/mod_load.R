#' Produce estimates of quantiles of electricity load for ISO New England (ISONE)
#'
#' \code{model_load} produces quantile estimates of electricity load.
#' This function is driven by yaml configuration files,
#' examples of which are included with this package.  The configuration
#' files also explain the meaning/import/use of the parameters.
#'
#' If writeFiles = 1 (in the modeling configuration file), then predictions
#' and results are written to the directories specified in the configuration
#' file.
#'
#' @param  load_cfg  a yaml file containing modeling parameters
#' @param  dts_cfg a yaml file containing general setup parameters
#' @export
#' @examples
#' \dontrun{
#' model_load('~/gefcom17/config/load_cfg.yaml',
#' '~/gefcom17/config/gen_and_dates_cfg.yaml')
#' }

model_load<-function(load_cfg  = '', dts_cfg = '') {

    if( load_cfg == '' | dts_cfg == '') {
        print('model_load needs 2 config files:')
        print('1) model parameters')
        print('2) validation/test dates')
        return()
    }

    cfg = yaml::read_yaml(load_cfg)
    dcfg = yaml::read_yaml(dts_cfg)

    doInt = cfg$doInt
    doImp = cfg$doImp
    doTest = cfg$doTest
    tstRunNum = cfg$tstRunNum

    bmme = sprintf('%d%d', dcfg$mth_data_beg, dcfg$mth_data_end)

    # 10 here, 100 for temps
    quants = (1:(cfg$num_quants-1))/cfg$num_quants

    if (cfg$rnd_seed > 0) set.seed(cfg$rnd_seed)

    ozns = readRDS(sprintf('%s/%s', dcfg$data_dir, dcfg$outfile))
    if(cfg$doInt)  ozns = fact_to_int(ozns)

    resDT = data.table()
    predDT = data.table()
    xtms_res = data.table()
    gt_res = data.table()

    zones = unique(ozns[,zone])
    if(cfg$zones[1] != 'all') zones = cfg$zones

    for( tzne in zones)
    {
        wxzn = tzne
        zns = ozns[zone == tzne,]
        zns = subset(zns, select=c(cfg$ld_mn_pred,cfg$ld_mn_tgt))

        valdt = data.table( dst = dcfg$val_dts$strt, dend = dcfg$val_dts$end)
        lstdt = data.table( eod = dcfg$val_dts$lst)

        # unlike weather code, don't have ability to iterate on month ranges
        mbeg=dcfg$mth_data_beg[1]
        mend=dcfg$mth_data_end[1]

        if (doTest == 1) {
            valdt = data.table( dst = dcfg$tst_dts$strt, ddend = cfg$tst_dts$end)
            lstdt = data.table( eod = dcfg$tst_dts$lst)
        }

        nrun = dim(valdt)[1]

        for ( i in  1:nrun) {

            eDate = lstdt[i,eod]
            fstDate = as.Date(eDate) - 365*dcfg$n_yrs_data
            zns = zns[ Date >= fstDate,]

            if( valdt[i, dend] < lstdt[i, eod]) {
                print('WARNING: possible overfitting, fitting IN SAMPLE')
            }

            val = zns[ Date > valdt[i, dst] & Date < valdt[i, dend],]
            trun = sprintf('RUN: %s - %d',tzne, year(val[1,Date]))
            print(trun, quote=F)

            vlmn = val[1,]$month
            vlyr = year(val[1,Date])

            if ( mbeg == mend)     { trn = zns[ month == mbeg & Date < eDate,] }
            else if ( mbeg > mend) { trn = zns[ (month >= mbeg | month <= mend) &
                                                  Date < eDate,] }
            else                   { trn = zns[ (month >= mbeg & month <= mend) &
                                                  Date < eDate,] }

            if( mbeg == 3 & mend == 11) { trn = trn[ month %in% c(3,4,5,10,11)] }

            vtact  = val[, tc]
            valLab = val[[cfg$ld_mn_tgt]]

            valt = copy(val) #need to keep val around
            valt[, c(cfg$rem_cols) := NULL]
            dval<-Matrix::sparse.model.matrix(~.-1, data = valt)
            dvald<-xgboost::xgb.DMatrix(data=dval,label=valLab)

            trnLab = trn[[cfg$ld_mn_tgt]]
            trnt = copy(trn)
            trnt[, c(cfg$rem_cols) := NULL]
            dtrn<-Matrix::sparse.model.matrix(~.-1, data = trnt)
            dtrnd<-xgboost::xgb.DMatrix(data=dtrn,label=trnLab)
            watchlist<-list(val=dvald,train=dtrnd)

            nrndm = switch(tzne,
                           'NH'   = cfg$nrnd_nh,
                           'ME'   = cfg$nrnd_me,
                           'VT'   = cfg$nrnd_vt,
                           'CT'   = cfg$nrnd_ct,
                           'RI'   = cfg$nrnd_ri,
                           'SEMA' = cfg$nrnd_sema,
                           'WCMA' = cfg$nrnd_wcma,
                           'NEMA' = cfg$nrnd_nema,
                           'ISO'  = cfg$nrnd_iso,
                           cfg$nrnd_def)

            if( doTest == 1) nrndm = round(nrndm*cfg$tstRndMult, 0)

            paramm <- list( objective  = cfg$obj, booster = cfg$bster,
                            eta = cfg$eta, max_depth = cfg$mxdpth,
                            subsample = cfg$subsamp,
                            colsample_bytree = cfg$colsamp,
                            min_child_weight = cfg$minchwt)

            mod <- xgboost::xgb.train( params = paramm, data = dtrnd,
                                      nrounds = nrndm,
                                      verbose = cfg$verb,
                                      maximize = FALSE,
                                      watchlist = watchlist,
                                      print_every_n = cfg$prnt_n)

            vpred = predict(mod, dvald)
            rmse  = sqrt(sum((vpred-valLab)^2)/length(valLab))

            print('Var regr', quote=F)

            dvalv<-Matrix::sparse.model.matrix(~.-1, data = valt)
            dtrnv<-Matrix::sparse.model.matrix(~.-1, data = trnt)

            trnLab2 = (trnLab - predict(mod, dtrnd))^2
            valLab2 = (valLab - vpred)^2

            dvaldv<-xgboost::xgb.DMatrix(data=dvalv, label=valLab2)
            dtrndv<-xgboost::xgb.DMatrix(data=dtrnv, label=trnLab2)
            watchlist<-list(val=dvaldv, train=dtrndv)

            nrndv = switch(tzne,
                           'NH'   = cfg$nrnd_nhv,
                           'ME'   = cfg$nrnd_mev,
                           'VT'   = cfg$nrnd_vtv,
                           'CT'   = cfg$nrnd_ctv,
                           'RI'   = cfg$nrnd_riv,
                           'SEMA' = cfg$nrnd_semav,
                           'WCMA' = cfg$nrnd_wcmav,
                           'NEMA' = cfg$nrnd_nemav,
                           'ISO'  = cfg$nrnd_isov,
                           cfg$nrnd_defv)

            if( doTest == 1) { nrndv = round(nrndv*cfg$tstRndMult, 0)}

            paramv <- list( objective  = cfg$obj, booster = cfg$bster, eta = cfg$eta_v,
                            max_depth = cfg$mxdpthv, subsample = cfg$subsampv,
                            colsample_bytree = cfg$colsampv,
                            min_child_weight = cfg$minchwtv)

            modV <- xgboost::xgb.train( params = paramv, data = dtrnd,
                                       nrounds = nrndv,
                                       verbose = 1,
                                       watchlist = watchlist,
                                       print_every_n = cfg$prnt_n,
                                       maximize = FALSE)

            predVar = predict(modV, dvaldv)
            rmsev  = sqrt(sum((predVar-valLab2)^2)/length(valLab2))

            #Hack quantiles assumes conditional mean == conditional median
            sds1 = sqrt(predVar)
            qdt = val
            qdt = setnames(qdt, 'rt_dem', 'act')

            for (quant in quants)
            {
                qnm = paste('q',quant,sep='')
                temp = vpred + qnorm(quant)*sds1
                qdt[, (qnm) := temp]
            }

            apb_xms = avgPinball(qdt, quants)

            tdt2 = subset(qdt, select=c(Date,hour,act,q0.1,q0.2,q0.3,q0.4,q0.5,
                                        q0.6,q0.7,q0.8,q0.9))
            tdt2[ , zone := tzne]
            tdt2[, alg := 'xms']
            predDT = rbind(predDT,tdt2)

            #Read in XGB est temps. Don't need VAR reg above, var comes from temp
            #quant. However, should check var reg w/ temps

            fnameX = sprintf('%s/%sX-%d-%d-%s.csv', dcfg$temp_dir, wxzn, vlmn,
                             vlyr, bmme)
            qtmpX = fread(fnameX)

            qdtX = subset(qdt, select=c(Date,act,hour))

            for (quant in quants) {

                if( quant == 0.1) { quant = 0.05}
                if( quant == 0.9) { quant = 0.95}

                qnm = paste('q',quant,sep='')
                valt[, tc := round(qtmpX[[qnm]],1)]
                dval<-Matrix::sparse.model.matrix(~.-1, data = valt)
                preds = predict(mod, dval)
                qdtX[, (qnm) := preds]

            }

            #Fake quantiles, there's just n where n is number of quantiles
            tqdt = subset(qdtX, select=c(Date,hour,act))
            qfoo = subset(qdtX, select=-c(Date,hour,act))
            qmat = t(apply(qfoo, 1, sort))
            qdt1 = as.data.table(qmat)
            setnames(qdt1, paste('q', c(1:9)*0.1,sep=''))

            qdtX = cbind(tqdt, qdt1)
            apb_xt = avgPinball(qdtX, quants)

            qdtX[ , zone := tzne]
            qdtX[, alg := 'xt']
            predDT = rbind(predDT,qdtX)

            #temp msig
            qdtXs = subset(qdtX, select = c(Date,hour,act))
            qfoo  = subset(qdtX, select = -c(Date,hour,act))

            valt2 = copy(valt)
            valt2[, tc := round(qtmpX[, q0.5], 1)]
            setcolorder(valt2, names(valt))
            dval<-Matrix::sparse.model.matrix(~.-1, data = valt2)
            predv = predict(modV, dval)
            sds2 = sqrt(predv)

            for (quant in quants){
                qnm = paste('q', quant,sep='')
                temp = qdtX[, q0.5] + qnorm(quant)*sds2
                qdtXs[, (qnm) := temp]
            }

            apb_xtms = avgPinball(qdtXs, quants)
            qdtXs[, `:=` (zone = tzne, alg = 'xtms')]
            predDT = rbind(predDT,qdtXs)

            #Error checking for xgb temp msig (predict mean/var on med predicted temp)
            qdtXs[, `:=` (pred_temp = round(qtmpX$q0.5, 1), act_temp = round(vtact, 1),
                          ovrfcst = as.integer( act < q0.1),
                          undfcst = as.integer( act > q0.9))]
            qdtXs[, dTemp := round(act_temp - pred_temp, 1)]

            # fix non use of apply
            tmpdt = data.table(q0.1 = pinball(qdtXs[, act], qdtXs[, q0.1], 0.1, 0))
            tqnts = quants[2:length(quants)]

            for ( quant in tqnts)
            {
                qnm = paste('q',quant,sep='')
                tmpdt[, (qnm) := pinball(qdtXs[, act], qdtXs[[qnm]], quant, 0)]
            }

            qdtXs[, avg_pnbl := rowMeans(tmpdt)]
            xtms_res = rbind(xtms_res, qdtXs)
            # Error checking for xgb temp msig

            #GAM temps COPY PASTA 1
            fnameG = sprintf('%s/%sG-%d-%d-%s.csv', dcfg$temp_dir, wxzn, vlmn,
                             vlyr, bmme)
            qtmpG = fread(fnameG)
            qdtG = subset(qdt, select=c(Date,hour,act))

            for (quant in quants)
            {
                qnm = paste('q',quant,sep='')
                valt[, tc := round(qtmpG[[qnm]],1)]
                dval<-Matrix::sparse.model.matrix(~.-1, data = valt)
                preds = predict(mod, dval)
                qdtG[, (qnm) := preds]
            }

            tqdt = subset(qdtG, select=c(Date,hour,act))
            qfoo = subset(qdtG, select=-c(Date,hour,act))
            qmat = t(apply(qfoo, 1, sort))
            qdt1 = as.data.table(qmat)
            setnames(qdt1, paste('q', c(1:9)*0.1,sep=''))

            qdtG = cbind(tqdt, qdt1)
            apb_gt = avgPinball(qdtG, quants)
            qdtG[ , zone := tzne]
            qdtG[, alg := 'gt']
            predDT = rbind(predDT,qdtG)

            # Error checking for gam temp
            qdtG[, `:=` (pred_temp = round(qtmpG$q0.5, 1), act_temp = round(vtact, 1),
                          ovrfcst = as.integer( act < q0.1),
                          undfcst = as.integer( act > q0.9))]
            qdtG[, dTemp := round(act_temp - pred_temp, 1)]

            tmpdt = data.table(q0.1 = pinball(qdtG[,act], qdtG[,q0.1], 0.1, 0))
            tqnts = quants[2:length(quants)]
            for ( quant in tqnts){
                qnm = paste('q',quant,sep='')
                tmpdt[, (qnm) := pinball(qdtG[,act], qdtG[[qnm]], quant, 0)]
            }
            qdtG[, avg_pnbl := rowMeans(tmpdt)]
            gt_res = rbind(gt_res, qdtG)

            #Avg temp
            qdtA1 = subset(qdt, select=c(Date,hour,act))
            for (quant in quants)
            {
                qnm = paste('q',quant,sep='')
                valt[, tc := round((qtmpX[[qnm]] + qtmpG[[qnm]])/2, 1)]

                dval<-Matrix::sparse.model.matrix(~.-1, data = valt)
                temp = predict(mod, dval)
                qdtA1[, (qnm) := temp]
            }

            tqdt = subset(qdtA1, select=c(Date, hour, act))
            qfoo = subset(qdtA1, select=-c(Date, hour, act))
            qmat = t(apply(qfoo, 1, sort))
            qdt1 = as.data.table(qmat)
            setnames(qdt1, paste('q', c(1:9)*0.1,sep=''))
            qdtA1 = cbind(tqdt, qdt1)
            apb_at = avgPinball(qdtA1, quants)
            qdtA1[ , zone := tzne]
            qdtA1[, alg := 'at']
            predDT = rbind(predDT,qdtA1)

            # Avg pred
            qdtA2 = subset(qdt, select=c(Date,hour,act))
            for (quant in quants)
            {
                qnm = paste('q',quant,sep='')
                qdtA2[, (qnm) := (qdtX[[qnm]] + qdtG[[qnm]])/2]
            }

            apb_ap = avgPinball(qdtA2, quants)
            qdtA2[, zone := tzne]
            qdtA2[, alg := 'ap']
            predDT = rbind(predDT,qdtA2)

            mnload = mean(val[,act])

            tdt = data.table(zone = tzne, val_mth = vlmn, val_yr = vlyr, mth_data_beg = mbeg,
                             mth_data_end = mend, rmse_mn = rmse, rmse_var = rmsev,
                             apb_xms = apb_xms, apb_xt = apb_xt, apb_xtms= apb_xtms,
                             apb_gt = apb_gt, apb_at = apb_at, apb_ap= apb_ap,
                             mean_load = mnload)

            resDT = rbind(resDT, tdt)
        } # nrun, val periods

        if (doImp == 1) {
            nimp  = xgboost::xgb.importance(dimnames(dtrn)[[2]], model=mod)
            nimpv = xgboost::xgb.importance(dimnames(dtrnv)[[2]], model=modV)
        }

        scols = names(resDT[, 6:14])
        tdt = as.data.table(dplyr::summarise( resDT[ zone == tzne],
                                              dplyr::across(scols, mean)))
        tdt[, `:=` (zone = tzne, val_mth = 999, val_yr = 999,
                    mth_data_beg = mbeg, mth_data_end = mend)]
        resDT = rbind(resDT, tdt)
    } #zones for loop

    # Had to submit MA as well. For ISO-level prediction, had choice
    # of using ISO-level model, or summing models of constituent zones.
    aggies = c('ISO2', 'MA')  #do in this order
    summ_cols = c('act', paste0('q0.', seq(1:9)))
    options(dplyr.summarise.inform = F)

    for (agg in aggies)
    {
        nuz = length(unique(predDT[,zone]))
        if( nuz < 2) {
            agg_msg = sprintf('skipping agg: %s', agg)
            print(agg_msg, quote = F)
            next
        }

        if( agg == 'ISO2') { toAgg = predDT[ zone != 'ISO']}
        else if( agg == 'MA') { toAgg = predDT[ zone %in% c('SEMA', 'WCMA', 'NEMA'),] }

        if ( dim(toAgg)[1] < 1) { print( 'nothing to agg'); next }

        yrs = unique(year(toAgg[,Date]))
        nyr = length(yrs)
        mth = unique(month(toAgg[,Date]))
        #assumes only 1 month predicted at a time

        for( yr in yrs) {

            tmpDT = toAgg[year(Date) == yr]
            tmpDT = dplyr::group_by(tmpDT, Date, hour, alg)
            qAll  = dplyr::summarise(tmpDT,
                                     dplyr::across(all_of(summ_cols), sum))
            qAll = as.data.table(qAll) # hate tibble
            qAll[, zone := agg]
            predDT = rbind(predDT, qAll)

            #calculate pinball for aggregations
            algs = unique(qAll[, alg])
            tlst = data.table()
            for (algo in algs)  tlst[, (algo) := avgPinball(qAll[alg == algo], quants)]

            mnload = mean(toAgg[ year(Date) == yr, act]) *
                          length(unique( toAgg[year(Date) == yr, zone] ))

            tdt = data.table(zone = agg, val_mth = vlmn, val_yr = yr,
                             mth_data_beg = mbeg, mth_data_end = mend, rmse_mn = 0,
                             rmse_var = 0, apb_xms = tlst$xms, apb_xt = tlst$xt,
                             apb_xtms= tlst$xtms, apb_gt = tlst$gt, apb_at = tlst$at,
                             apb_ap= tlst$ap, mean_load = mnload)

            resDT = rbind(resDT, tdt)
        } #iterate over years in study

        scols = names(resDT[, 6:14])
        tdt = as.data.table(dplyr::summarise( resDT[ zone == agg],
                                              dplyr::across(scols, mean)))
        tdt[, `:=` (zone = agg, val_mth = 999, val_yr = 999, mth_data_beg = mbeg,
                   mth_data_end = mend)]
        resDT = rbind(resDT, tdt)

    } #agg loop for MA and ISO

    for ( i in 1:8)
    {
        a = paste('q', quants[i], sep='')
        b = paste('q', quants[i+1], sep='')

        if( sum(predDT[[b]] < predDT[[a]]) > 0) {
            print( sprintf('QUANT CROSSING: %s < %s', b, a))
        }
    }

    print(resDT)

    if (cfg$writeFiles) {

        tgt_mth = month.abb[vlmn]
        xtms_res[, dow := wday(Date)]
        xtms_res[, Date := as.character(Date)]
        gt_res[, dow := wday(Date)]
        gt_res[, Date := as.character(Date)]
        xtms_res = rbind(xtms_res, gt_res)

        gfn = sprintf('%s/xtms_gt_validation_res_%s_%s.csv', dcfg$res_dir,
                      tgt_mth, bmme)
        fwrite(xtms_res, gfn)

        val_test = 'validation'

        if(doTest) {
            val_test = 'test'
        } else {
            resfn = sprintf('%s/summ_val_results_%s_%s.csv', dcfg$res_dir,
                            tgt_mth, bmme)
            fwrite(resDT, resfn)
        }

        predfn = sprintf('%s/all_preds_%s_%s_%s.csv', dcfg$sub_dir, tgt_mth,
                         val_test, bmme)
        fwrite(predDT, predfn)

        if(doImp) {
          nmfn = sprintf('%s/nimp_mean_load_%s.csv', dcfg$res_dir, tgt_mth)
          fwrite(nimp, nmfn)
          nvfn = sprintf('%s/nimp_var_load_%s.csv',  dcfg$res_dir, tgt_mth)
          fwrite(nimpv, nvfn)
        }

    }

    if( doTest == 1) {

        if( tmth == 3) {
            newPred2 = predDT[ (alg == 'gt' & zone != 'CT') | (alg == 'xtms' &
                                zone == 'CT'),]
            newPred2 = newPred2[ zone != 'ISO']
            newPred2[ zone == 'ISO2', zone := 'ISO']
        }
        else if ( tmth == 4) {
            newPred2 = predDT[ alg == 'xtms' & zone != 'ISO2', ]
        }
        else{
            newPred2 = predDT[alg == 'gt',]
            newPred2 = newPred2[zone != 'ISO2']
        }

        #Write to file
        newPred3 = subset(newPred2, select=-c(act, alg))
        newPred3[, Date := as.character(Date)]
        newPred3[ , Date := gsub('-','/', Date, fixed=T)]
        setnames(newPred3, 'hour', 'Hour')
        oldn = paste('q',quants,sep='')
        newn = paste('Q',c(1:9)*10,sep='')
        setnames(newPred3,oldn,newn)
        pfn = sprintf('%s/test_%d_%d.csv', dcfg$subDir, tmth, cfg$tstRunNum)
        fwrite(newPred3, pfn)
    }
}
