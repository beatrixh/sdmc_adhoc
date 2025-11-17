## ------------------------------------------------------------------------------------------##
# Author: Beatrix Haddock
# Date: 2025-11-14
# Purpose: Calculate one IC50/80 per mAb-isolate pair from replicate estimates
## ------------------------------------------------------------------------------------------##

library(survival)
library(dplyr)
library(tidyr)
library(psych)

# read in data
df = read.csv(
  "N:/cavd/Studies/cvd914/qdata/LabData/NAB/CAVD914_NAB07_U6_20250924.txt",
  sep='\t'
)

# formatting
colnames(df) <- tolower(colnames(df))
df <- df %>%
  rename(
    tetO_HIV_Env_IMC_name = isolate,
    mab_name = guspec,
    val = value
  )


# parse out isolate name
df$isolate <- sapply(strsplit(df$tetO_HIV_Env_IMC_name, "\\."), function(x) substr(x[2], 2, nchar(x[2])))

# parse out pubid
df$pubid <- sapply(strsplit(df$isolate, "_"), function(x) paste(x[1], x[2], sep = "-"))

# id vars
df$mab_isolate_titer <- paste(df$mab_name, df$isolate, df$value_label, sep = "|")
df$mab_isolate_pubid_titer <- paste(df$mab_name, df$isolate, df$pubid, df$value_label, sep = "|")

# subset data

df <- df[df$value_label %in% c("TITER_50_CURVE", "TITER_80_CURVE"),
         c("dilution_min", "mab_isolate_pubid_titer", "mab_isolate_titer",
           "mab_name", "isolate", "pubid", "value_label", "val")]

# categorical vals
df$isolate <- factor(df$isolate)
df$mab_name <- factor(df$mab_name)
df$pubid <- factor(df$pubid)

# create flag for rows above/below cutoff
df$lod = ""
df$lod[which(grepl("^>", df$val))] = ">"
df$lod[which(grepl("^<", df$val))] = "<"

# remove < and > from 'val' column
df$val = as.numeric(sub("^<|>", "", df$val))

# drop zeros; lab meant to exclude these
df <- df[df$val > 0, ]

# calculate mean titers using tobit as needed
titers = c()

for(iso in levels(df$isolate)){ #isolates/pubids map 1:1
  for(mab in levels(df$mab_name)){
    for(t in c("TITER_80_CURVE", "TITER_50_CURVE")){
      # subset to one isolate-mab pair and either IC50 or 80
      ss = subset(df, isolate==iso & mab_name==mab & value_label==t)
      
      if(all(ss$lod=="")){
        # all-numeric case
        titers = rbind(
          titers,
          data.frame(
            isolate=iso,
            mab_name=mab,
            pct_neut=t,
            titer=geometric.mean(ss$val),
            case="all numeric",
            n_numeric=nrow(ss),
            n_above=0,
            n_below=0
          )
        )
      } else if(all(ss$lod==">")){
        # all above upper threshold case
        titers = rbind(
          titers,
          data.frame(
            isolate=iso,
            mab_name=mab,
            pct_neut=t,
            titer=">25",
            case="all above",
            n_numeric=0,
            n_above=nrow(ss),
            n_below=0
          )
        )
      } else if(all(ss$lod=="<")){
        # all below lower threshold case
        titers = rbind(
          titers,
          data.frame(
            isolate=iso,
            mab_name=mab,
            pct_neut=t,
            titer="<0.011431184270690446",
            case="all below",
            n_numeric=0,
            n_above=0,
            n_below=nrow(ss)
          )
        )
      } else if(any(ss$lod==">") && !any(ss$lod=="<")){
        # some above threshold, some numeric, none below
        n_above <- sum(ss$lod == ">")
        ss$surv <- Surv(log(ss$val), ifelse(ss$lod==">", 0, 1))
        fit = survreg(ss$surv~1, dist="gaussian")
        mu = exp(coef(fit))
        titers = rbind(
          titers,
          data.frame(
            isolate=iso,
            mab_name=mab,
            pct_neut=t,
            titer=mu,
            case="some above",
            n_numeric=nrow(ss)-n_above,
            n_above=n_above,
            n_below=0
          )
        )
      } else if (any(ss$lod=="<") && !any(ss$lod==">")){
        # some below threshold, some numeric, none above
        n_below <- sum(ss$lod == "<")
        ss$surv <- Surv(log(ss$val), ifelse(ss$lod=="<", 0, 1), type="left")
        fit = survreg(ss$surv~1, dist="gaussian")
        mu = exp(coef(fit))
        titers = rbind(
          titers,
          data.frame(
            isolate=iso,
            mab_name=mab,
            pct_neut=t,
            titer=mu,
            case="some below",
            n_numeric=nrow(ss)-n_below,
            n_above=0,
            n_below=n_below
          )
        )
      } else {
        # if any rows have both above and below
        n_above <- sum(ss$lod == ">")
        n_below <- sum(ss$lod == "<")
        titers = rbind(
          titers,
          data.frame(
            isolate=iso,
            mab_name=mab,
            pct_neut=t,
            titer="CAUTION, BOTH ABOVE AND BELOW",
            case="CAUTION, BOTH ABOVE AND BELOW",
            n_numeric=nrow(ss)-n_above-n_below,
            n_above=n_above,
            n_below=n_below
          )
        )
      }
    }
  }
}

savedir <- "N:/vtn/lab/SDMC_labscience/studies/VISC/Ochsenbauer_914/assays/C-C_nAb/misc_files/interim_data/"
today <- as.character(Sys.Date())

write.table(
  titers,
  paste0(savedir, "tobit_titers_", today, ".txt"),
  sep='\t',
  row.names=FALSE,
  quote = FALSE
)
