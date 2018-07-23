# Crins et al PTLD (post-transplant lymphoproliferative disease)
# and mortality data

require("bayesmeta")
data("CrinsEtAl2014")

ptld.exp <- c(1,1,NA,NA,0,NA)
ptld.cont <- c(1,0,NA,NA,0,NA)

death.exp  <- c(4,NA,1,NA,4,2)
death.cont <- c(3,NA,3,NA,3,3)

CrinsEtAl2014 <- cbind(CrinsEtAl2014,
                       "exp.PTLD.events"  = ptld.exp[c(4,3,5,2,1,6)],
                       "cont.PTLD.events" = ptld.cont[c(4,3,5,2,1,6)],
                       "exp.death.events"  = death.exp[c(1,2,4,3,5,6)],
                       "cont.death.events" = death.cont[c(1,2,4,3,5,6)])

print(CrinsEtAl2014[,c(1,16,12,17,15)])

print(CrinsEtAl2014[,c(1,18,12,19,15)])


CrinsEtAl2014$IL2RA <- NULL
CrinsEtAl2014$CNI <- NULL
CrinsEtAl2014$MMF <- NULL
CrinsEtAl2014$exp.AR.events <- NULL
CrinsEtAl2014$exp.SRR.events <- NULL
CrinsEtAl2014$cont.AR.events <- NULL
CrinsEtAl2014$cont.SRR.events <- NULL


write.table(CrinsEtAl2014, file = "data-raw/dat.Crins2014.txt", sep = ",")
require("metafor")

es.ptld <- escalc(measure="OR",
                  ai=exp.PTLD.events,  n1i=exp.total,
                  ci=cont.PTLD.events, n2i=cont.total,
                  slab=publication,
                  data=CrinsEtAl2014)
forest(rma(es.ptld), main="PTLD")


es.death <- escalc(measure="OR",
                   ai=exp.death.events,  n1i=exp.total,
                   ci=cont.death.events, n2i=cont.total,
                   slab=publication,
                   data=CrinsEtAl2014)
forest(rma(es.death), main="mortality")
