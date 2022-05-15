# run the NT glyphis TMB code - this one has mtDNA included

require(TMB)

compile("gnt1.cpp") # ,"-O1 -g",DLLFLAGS="")
dyn.load("gnt1.so")

###################
# set up the data #
###################

ckmr.df <- read.table("glyph_NT_ckmr.dat",header=T)

# code for kin types: up = 0, hsp = 1, fsp = 2

ckmr.df$kcode <- rep(0)
ckmr.df$kcode[ckmr.df$kin.type == 'hsp'] <- 1
ckmr.df$kcode[ckmr.df$kin.type == 'fsp'] <- 2
data <- as.list(ckmr.df)

# tmax: number of cohorts including min(c1,c2)-1

crng <- c(min(c(data$c1,data$c2)),max((c(data$c1,data$c2))))
tmax <- length(crng[1]:crng[2])+1
tzero <- crng[1]-1
data$tzero <- tzero
data$tmax <- tmax
data$nriv <- data$nsex <- 2
data$p_hsp <- 0.92

#########################
# set up the parameters #
#########################

pars <- list()

## xi - logit-scake omega s * r of them

pars$xi <- matrix(c(3,3,3,3),ncol=2,byrow=T)

# test it works

xi <- pars$xi
omega <- array(dim=c(2,2,2))
for(s in 1:2) {

  omega[s,1,1] <- 1/(1+exp(-xi[s,1]))
  omega[s,1,2] <- 1-omega[s,1,1]
  omega[s,2,2] <- 1/(1+exp(-xi[s,2]))
  omega[s,2,1] <- 1-omega[s,2,2] 

}

## lambda - log-scale growth rate
pars$lambda <- c(0,0)

## log(N0) for each river (total males and females)

pars$lN0 <- c(log(700),log(700))

## iphi - logit-scale survival (r)

pars$iphi <- c(2,2)

## lnu - log-scale river-specific litter effect

pars$lnu <- c(log(2),log(2))

## izeta - logit-scale sex ratio for each river

pars$izeta <- c(0,0)

## itheta - logit-scale multiple paternity parameter (maternal par)

pars$itheta <- -4

## lgamma - log-scale multiple female breeding parterns per year (paternal par)

pars$lgamma <- -5

########################
# create the AD object #
########################

inv.logit <- function(x){return(1/(1+exp(-x)))}

# 1. Just estimate abundance and survival

map.obj1 <- list(xi=rep(factor(NA),length(pars$xi)),
    lambda=rep(factor(NA),length(pars$lambda)),
    izeta=rep(factor(NA),length(pars$izeta)),
    lnu= rep(factor(NA),length(pars$lnu)),
    itheta=factor(NA ),lgamma=factor(NA ))

obj1 <- MakeADFun(data=data,parameters=pars,map=map.obj1,DLL='gnt1')

fit1 <- do.call("optim", obj1)
sdfit1 <- sdreport(obj1)
sdfit1
exp(fit1$par[1:2])
inv.logit(fit1$par[3:4])

# 2. 
# a. fixed survival at same values
# b. fix within-cohort parameters
# c. fix the spatial reproductive parameters

# bring in the sex ratio parameter now

map.obj2 <- list(xi=rep(factor(NA),length(pars$xi)),
    lambda=rep(factor(NA),length(pars$lambda)),
    iphi=rep(factor(NA),length(pars$iphi)),
    lnu= rep(factor(NA),length(pars$lnu)),
    itheta=factor(NA ),lgamma=factor(NA ))

obj2 <- MakeADFun(data=data,parameters=pars,map=map.obj2,DLL='gnt1')

fit2 <- do.call("optim", obj2)
sdfit2 <- sdreport(obj2)
sdfit2

exp(fit2$par[1:2])
inv.logit(fit2$par[3:4])


# 3.
# a. same as 2 but estimate spatial reproductive stuff

map.obj3 <- list(lambda=rep(factor(NA),length(pars$lambda)),
    iphi=rep(factor(NA),length(pars$iphi)),
    lnu= rep(factor(NA),length(pars$lnu)),
    itheta=factor(NA ),lgamma=factor(NA ))


obj3 <- MakeADFun(data=data,parameters=pars,map=map.obj3,DLL='gnt1')

fit3 <- do.call("optim", obj3)
sdfit3 <- sdreport(obj3)
sdfit3
exp(fit3$par[c(5,6)])
inv.logit(fit3$par[-c(1:6)])

xi <- matrix(fit3$par[1:4],ncol=2,byrow=T)

omega <- array(dim=c(2,2,2))
for(s in 1:2) {

  omega[s,1,1] <- 1/(1+exp(-xi[s,1]))
  omega[s,1,2] <- 1-omega[s,1,1]
  omega[s,2,2] <- 1/(1+exp(-xi[s,2]))
  omega[s,2,1] <- 1-omega[s,2,2] 

}

# 4. 
# a. fix sex ratio
# b. estimate xi

map.obj4 <- list(lambda=rep(factor(NA),length(pars$lambda)),
    iphi=rep(factor(NA),length(pars$iphi)),
    lnu= rep(factor(NA),length(pars$lnu)),
    izeta= rep(factor(NA),length(pars$izeta)),
    itheta=factor(NA ),lgamma=factor(NA ))

obj4 <- MakeADFun(data=data,parameters=pars,map=map.obj4,DLL='gnt1')

fit4 <- do.call("optim", obj4)
sdfit4 <- sdreport(obj4)
sdfit4

xi <- matrix(fit4$par[1:4],ncol=2,byrow=T)

omega <- array(dim=c(2,2,2))
for(s in 1:2) {

  omega[s,1,1] <- 1/(1+exp(-xi[s,1]))
  omega[s,1,2] <- 1-omega[s,1,1]
  omega[s,2,2] <- 1/(1+exp(-xi[s,2]))
  omega[s,2,1] <- 1-omega[s,2,2] 

}3.599073e-01 

# 5. 
# a. fix sex ratio
# b. estimate xi
# c. estimate nu

map.obj5 <- list(lambda=rep(factor(NA),length(pars$lambda)),
    iphi=rep(factor(NA),length(pars$iphi)),
    izeta= rep(factor(NA),length(pars$izeta)),
    itheta=factor(NA ),lgamma=factor(NA ))

obj5 <- MakeADFun(data=data,parameters=pars,map=map.obj5,DLL='gnt1')

fit5 <- do.call("optim", obj5)
sdfit5 <- sdreport(obj5)
sdfit5

xi <- matrix(fit5$par[1:4],ncol=2,byrow=T)

omega <- array(dim=c(2,2,2))
for(s in 1:2) {

  omega[s,1,1] <- 1/(1+exp(-xi[s,1]))
  omega[s,1,2] <- 1-omega[s,1,1]
  omega[s,2,2] <- 1/(1+exp(-xi[s,2]))
  omega[s,2,1] <- 1-omega[s,2,2] 

}

summary(sdfit5)
