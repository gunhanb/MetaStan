library(MetaStan)
data("dat.Boucher2016.pairwise")
set.seed(12234)
dat.topi <- create_MetaStan_dat(dat = dat.Boucher2016.pairwise,
                              armVars = c(responders="r", sampleSize="n"))

ma.topi <- meta_stan(data=dat.Boucher,
                     likelihood = "binomial",
                     re=TRUE,
                     ncp=TRUE, mu_prior=c(0,10),
                     theta_prior=c(0,2.5),
                     tau_prior=0.5,
                     tau_prior_dist="half-normal",
                     chains=4, iter=4000,
                     warmup=2000)
ma.topi

myplot = forest_plot(x = ma.topi,
                     xlab="log-OR",
                     labels=dat.Boucher2016.pairwise$study)

# initialize plot
pdf("myforestplot.pdf", width=480, height=480)

# make plot
myplot

# save plot
dev.copy(pdf, "myforestplot.pdf")
dev.off()


posterior <- as.data.frame(rstan::extract(ma.topi$fit, pars="theta"))
hist(posterior$theta)


mreg.topi <- meta_stan(data=dat.topi, likelihood="binomial",
                       mu_prior=c(0,10),
                       theta_prior=c(0,2.5),
                       tau_prior=0.5,
                       tau_prior_dist="half-normal", mreg=TRUE,
                       cov=as.numeric(dat.Boucher2016.pairwise$duration=="long"),
                       beta_prior=c(0,10))
mreg.topi


data("dat.Boucher2016")
datMBMA <- create_MetaStan_dat(dat=dat.Boucher2016,
                              armVars=c(dose="d", responders="r",
                                        sampleSize="n"),
                              nArmsVar="nd")


MBMA.emax <- MBMA_stan(data=datMBMA, likelihood="binomial",
                       dose_response="emax", Emax_prior=c(0, 10),
                       ED50_prior_dist="functional",
                       tau_prior_dist="half-normal", tau_prior=0.5)


plot(MBMA.emax) + ggplot2::xlab("topiramate dose (mg)") +
  ggplot2::ylab("paresthesia probability")


