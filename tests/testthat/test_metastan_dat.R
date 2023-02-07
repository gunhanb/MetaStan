################################################################################
# actual data example: RCTs (two treatments, no covariables):

require("bayesmeta")
data("CrinsEtAl2014")

crins.ar <- metastan_data(id=c(publication, publication),
                          treat=rep(c("control", "treatment"),
                                    each=nrow(CrinsEtAl2014)),
                          type="binomial",
                          count=c(cont.AR.events, exp.AR.events),
                          total=c(cont.total, exp.total),
                          data=CrinsEtAl2014)
crins.ar$outcome




# 2 studies, normal endpoint, no treatment/control, no covariable:
xx <- metastan_data(id=c("Mueller","Mueller", "Meier","Meier","Meier"),
                    type="normal",
                    y=c(1,2,3,3,4),
                    se=c(1,1,2,2,2))
xx

# 2 studies, binary endpoint, treatment + control, "dose" covariable:
xx <- metastan_data(id=c("Mueller","Mueller", "Meier","Meier","Meier"),
                   treatment=c("c","t","c","t","t"),
                   x=cbind("dose"=c(10,20,10,20,50)),
                   type="binomial",
                   count=c(1,2,2,2,4),
                   total=c(10,10,15,20,30))
xx$outcome


dat_long = create_MetaStan_dat(dat = dat.Berkey1995,
                               armVars = c(responders = "r",
                                           sampleSize = "n"),
                               nArmsVar = "nd")
head(dat_long$data_long)

# check out the actual "plain" returned object:
print.default(xx)




# 2 studies, normal endpoint, one ("dose") covariable:
xx <- MetaStanData(id=c("Mueller","Mueller", "Meier","Meier","Meier"),
                   type="normal",
                   y=c(1,2,3,3,4),
                   se=c(1,1,2,2,2),
                   x=cbind("dose"=c(10,50,10,20,50)))
xx

# 2 studies, normal endpoint, one covariable (vector, named "x" by default):
xx <- MetaStanData(id=c("Mueller","Mueller", "Meier","Meier","Meier"),
                   type="normal",
                   y=c(1,2,3,3,4),
                   se=c(1,1,2,2,2),
                   x=c(10,50,10,20,50))
xx

# 2 studies, normal endpoint, two covariables:
xx <- MetaStanData(id=c("Mueller","Mueller", "Meier","Meier","Meier"),
                   type="normal",
                   y=c(1,2,3,3,4),
                   se=c(1,1,2,2,2),
                   x=matrix(c(10,50,10,20,50, 2,2,3,3,3), nrow=5, ncol=2,
                            dimnames=list(NULL, c("dose", "duration"))))
xx

# 2 studies, binary endpoint, no treatment/control, no covariable:
xx <- MetaStanData(id=c("Mueller","Mueller", "Meier","Meier","Meier"),
                   type="binomial",
                   count=c(1,2,2,2,4),
                   total=c(10,10,15,20,30))
xx




################################################################################
# actual data examples: dose-finding studies (one treatment, several doses):

#require("devtools")
#devtools::install_github("wviechtb/metadat",  lib="~/temp/test")
library("metadat")

data("dat.ursino2021")
head(dat.ursino2021)

# specify arguments as vectors:
sorafenib <- MetaStanData(id=dat.ursino2021$study,
                          type="binomial",
                          count=dat.ursino2021$events,
                          total=dat.ursino2021$total,
                          x=cbind("dose"=dat.ursino2021$dose))

# specify arguments as columns of a data frame:
sorafenib <- MetaStanData(id=study,
                          type="binomial",
                          count=events,
                          total=total,
                          x=cbind("dose"=dose),
                          data=dat.ursino2021)

sorafenib

print.default(sorafenib)


####################
# 2nd example:

data("dat.roever2022")
head(dat.roever2022)

# specify arguments as vectors:
irinotecan <- metastan_data(id=dat.roever2022$study,
                           type="binomial",
                           count=dat.roever2022$events,
                           total=dat.roever2022$total,
                           x=cbind("dose"=dat.roever2022$dose))

# specify arguments as columns of a data frame:
irinotecan <- MetaStanData(id=study,
                           type="binomial",
                           count=events,
                           total=total,
                           x=cbind("dose"=dose, "year"=year),
                           data=dat.roever2022)
irinotecan



################################################################################
# actual data example: proportions (no treatments, no covariables,
# log-odds estimate to be used as a MAP prior):

neuenschwander2010 <- MetaStanData(id=c("Feagan2005", "Rutgeerts2005a",
                                        "Rutgeerts2005b", "VanAssche2006"),
                                   type="binomial",
                                   count=c(9, 18, 7, 6),
                                   total=c(63, 121, 123, 56))


################################################################################
# actual data example: prevalence study (no treatments, no covariables):

source("CrisafulliEtAl2020.R")
head(dat.crisafulli2020)

# specify arguments as vectors:
dmd <-  MetaStanData(id=dat.crisafulli2020$study,
                     type="binomial",
                     count=dat.crisafulli2020$cases,
                     total=dat.crisafulli2020$total)

# specify arguments as columns of a data frame:
dmd <-  MetaStanData(id=study,
                     type="binomial",
                     count=cases, total=total,
                     data=dat.crisafulli2020)
dmd
dmd$call
print(dmd, n=10)
print.default(dmd)


################################################################################
# actual data example: RCTs (two treatments, no covariables):

require("bayesmeta")
data("CrinsEtAl2014")

crins.ar <- metastan_data(id=c(publication, publication),
                         treat=rep(c("control", "treatment"),
                                   each=nrow(CrinsEtAl2014)),
                         type="binomial",
                         count=c(cont.AR.events, exp.AR.events),
                         total=c(cont.total, exp.total),
                         data=CrinsEtAl2014)
crins.ar$outcome
print(crins.ar, n=12)

# as above, including a covariable:
crins.ar <- MetaStanData(id=c(publication, publication),
                         treat=rep(c("control", "treatment"),
                                   each=nrow(CrinsEtAl2014)),
                         x=cbind("followup"=rep(followup, 2)),
                         type="binomial",
                         count=c(cont.AR.events, exp.AR.events),
                         total=c(cont.total, exp.total),
                         data=CrinsEtAl2014)
crins.ar
print.default(crins.ar)

# as above, using "PTLD" outcome:
crins.ptld <- MetaStanData(id=c(publication, publication),
                           treat=rep(c("control", "treatment"), each=3),
                           type="binomial",
                           count=c(cont.PTLD.events, exp.PTLD.events),
                           total=c(cont.total, exp.total),
                           data=CrinsEtAl2014[complete.cases(CrinsEtAl2014[,c(12,17)]),])
crins.ptld

