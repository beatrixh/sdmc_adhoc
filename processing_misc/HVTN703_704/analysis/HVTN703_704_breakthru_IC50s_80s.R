rm(list=ls())
library(survival)

# read in data 
dat703 = read.csv("T:/vaccine/p703/s573/qdata/VTN703_breakthrough_NAb_20200722.txt", sep="\t")
dat704 = read.csv("T:/vaccine/p704/s670/qdata/VTN704_breakthrough_NAb_20200910.txt", sep="\t")
dat = rbind(dat703, dat704)

dat$isolate = factor(dat$isolate)
dat = subset(dat, select=c("isolate","poscrit","initdilution","titer","labid"))
dat$lod = ""
dat$lod[which(grepl("^>", dat$titer))] = ">"
dat$lod[which(grepl("^<", dat$titer))] = "<"
dat$titer = as.numeric(sub("^<|>", "", dat$titer))


# generate the concordance results
res = c()

for( iso in levels(dat$isolate) ) {
  for( crit in c(50, 80) ){
    ss = subset(dat, poscrit==crit & isolate==iso)   
    ss5 = subset(ss, initdilution==5)
    ss100 = subset(ss, initdilution==100)
    
    if( all(ss$lod=="") ){
      mu = exp(mean(log(ss$titer)))
      nconcord = sum(abs(log(ss$titer) - log(mu)) <= log(3))
      n = nrow(ss)
      res = rbind(res, data.frame(isolate=iso, poscrit=crit, mu=mu, min=min(ss$titer), max=max(ss$titer), nconcord=nconcord, n=n, cat="all"))
    } else if( all(ss5$lod=="") ) {
      ss = ss5
      mu = exp(mean(log(ss$titer)))
      nconcord = sum(abs(log(ss$titer) - log(mu)) <= log(3))
      n = nrow(ss)
      res = rbind(res, data.frame(isolate=iso, poscrit=crit, mu=mu, min=min(ss$titer), max=max(ss$titer), nconcord=nconcord, n=n, cat="start conc 5 ug/ml"))
    } else if( all(ss100$lod=="") && nrow(ss100) > 0 ) {
      ss = ss100
      mu = exp(mean(log(ss$titer)))
      nconcord = sum(abs(log(ss$titer) - log(mu)) <= log(3))
      n = nrow(ss)
      res = rbind(res, data.frame(isolate=iso, poscrit=crit, mu=mu, min=min(ss$titer), max=max(ss$titer), nconcord=nconcord, n=n, cat="start concentration 100 ug/ml"))
      
    } else if(all(ss$lod==">")) {
      if(nrow(ss100) > 0){
        ss = ss100
      }
      ss$surv = Surv(log(ss$titer), ifelse(ss$lod==">", 0, 1))
      fit = survreg(ss$surv~1, dist="gaussian")
      mu = exp(coef(fit))
      nconcord = n = nrow(ss)
      res = rbind(res, data.frame(isolate=iso, poscrit=crit, mu=mu, min=min(ss$titer), max=max(ss$titer), nconcord=nconcord, n=n, cat="all >"))
      
    } else if(all(ss5$lod==">") && !any(ss100$lod=="<")) {
      ss = ss100
      ss$surv = Surv(log(ss$titer), ifelse(ss$lod==">", 0, 1))
      fit = survreg(ss$surv~1, dist="gaussian")
      mu = exp(coef(fit))
      n = nrow(ss)
      nconcord = sum( abs(log(ss$titer) - log(mu)) <= log(3) )
      
      res = rbind(res, data.frame(isolate=iso, poscrit=crit, mu=mu, min=min(ss$titer), max=max(ss$titer), nconcord=nconcord, n=n, cat="mix at 100 ug/ml"))
    } else {
      print(iso)
      print(crit)
      print(table(ss$initdilution, ss$lod))
    }  
  }
}


res$qc_flag = FALSE
res$qc_flag[(res$isolate == "H704_1739_210_RE_p002s_CtoTat8") & (res$poscrit == 80)] = TRUE

res$qc_note = ""
res$qc_note[(res$isolate == "H704_1739_210_RE_p002s_CtoTat8") & (res$poscrit == 80)] = "No estimates with initdilution==100. All estimates from initdilution==5 were >5"

savedir = 'N:/vtn/lab/SDMC_labscience/studies/HVTN/HVTN703_704/analysis/nAb/ic_estimates/'
write.table(res, paste0(savedir,"HVTN703_704_IC50_80_estimates_2026-05-07.txt"), sep="\t", row.names=FALSE)

q(save="no")

