library(relsurv)
data(slopop)
try(source("/home/klemen/work/logrank_test/gen_net2.r"))
source("/home/klemen/work/Pseudovalues/Simulations/rssurvPseudo.r")
dyn.load("/home/klemen/work/Pseudovalues/Simulations/pseudoWH_asy.so")
dyn.load("/home/klemen/work/Pseudovalues/Simulations/pseudoWH_pre.so")
dyn.load("/home/klemen/work/Pseudovalues/Simulations/sindex.so")
dyn.load("/home/klemen/work/Pseudovalues/Simulations/pseudoWH_none.so")
dyn.load("/home/klemen/work/Pseudovalues/Simulations/pseudoWH_asy2.so")
dyn.load("/home/klemen/work/Pseudovalues/Simulations/pseudoWH_pre2.so")
dyn.load("/home/klemen/work/Pseudovalues/Simulations/pseudoWH_none_LT.so")
dyn.load("/home/klemen/work/Pseudovalues/Simulations/pseudoWH_asy2_LT.so")
dyn.load("/home/klemen/work/Pseudovalues/Simulations/pseudoWH_pre2_LT.so")
# dyn.load("/home/klemen/work/Pseudovalues/Simulations/testna.so")
relsurv.rs.survPseudo <- function(formula = formula(data), data = parent.frame(), ratetable = slopop, na.action, eval.times, method = "KM", var.method = "precise", varSk.method = "precise",conf.type = "plain", conf.int = 0.95)
{ 
  call <- match.call()
  rform <- relsurv:::rformulate(formula, data, ratetable, na.action)
#   print(names(rform))
#   print(rform$data)
#   print(rform$X)
#   print(rform$m)
#   print(rform$Terms)
#   print(rform$mvalue)
  
  if (missing(eval.times)){eval.times <- rform$Y}
  eval.times <- sort(unique(eval.times))
  conf.type <- match.arg(conf.type, c("plain","log","log-log"))
  method <- match.arg(method,c("KM","lifetable"))
  var.method <- match.arg(var.method,c("precise","asymptotic (ind.)","none"))
  varSk.method <- match.arg(varSk.method,c("precise","reduced sample (KM)"))
  varSk.methodN <- match(varSk.method,c("precise","reduced sample (KM)"))
  if (max(eval.times) > max(rform$Y))
    warning("You want to estimate survival S(t) for t that is strictly bigger than all observed times in your data")
  ### DO KAM LAHKO RACUNAS, PSEVDO VREDNOSTI DO PREDZADNJEGA CASA, varianco pa do predpredzadnjega casa
  obstimes <- rform$Y
  inx <- order(obstimes,1-rform$status)
  # print(inx)
  obstimes <- obstimes[inx]
  status <- rform$status[inx]
  # print(obstimes)
  # print(status)
  N <- rform$n # = length(obstimes)
  time <- unique(obstimes)
  Y <- rep(N,length(time)) # number at risk
  D <- rep(0,length(time)) # number of events
  C <- rep(0,length(time)) # number of censored
  for (i in 1:length(time)){
    Y[i] <- N - length(which(obstimes < time[i])) # number at risk
    ind_i <- which(obstimes == time[i])
    D[i] <- sum(status[ind_i])
    C[i] <- sum(1-status[ind_i])
  }
  # print(rbind(time,Y,C,D))
  Y <- Y[D > 0]
  C <- C[D > 0]
  # print(rform$R[inx,])
  time <- time[D > 0] # unique event times
  D <- D[D > 0]
  # print(rbind(time,Y,C,D))
  NT <- length(time)
  # the part below is extracted from the function sindex (package prodlim)
  neval <- length(eval.times)
  ind <- .C("sindex", index = integer(neval), as.double(time),
            as.double(eval.times), as.integer(NT),
            as.integer(neval), as.integer(FALSE))$index
  ind = ind - 1
  # pseudo observations have to be calculated only at the indices above
  # IMPORTANT
  # S_p are calculated outside - can improve the speed
  Sp <- rep(NA, N * neval) # vector of dimension (number of patients x length(eval.times))
  for (t in 1:neval){
    # first N elements at eval.times[1], the next N elements at eval.times[2],...
    Sp[((t-1)*N + 1):(t*N)] <- relsurv:::exp.prep(rform$R[inx,],eval.times[t],ratetable)
  }
  # print(matrix(Sp,nrow=N,byrow=FALSE))
  # print(Sp)
  if (method == "KM"){
    if (var.method == "none"){
      # estimates only survival, but not the variance
      temp_res <- .C("pseudoWH_none", Y = as.double(Y), D = as.double(D), time = as.double(time), obsT = as.double(obstimes), status = as.double(status), S_E = double(neval), Sp = as.double(Sp), N = as.integer(N), NT = as.integer(NT), neval = as.integer(neval), ind = as.integer(ind))
    }
    else if (var.method == "precise"){
      # precise formula for the variance
      temp_res <- .C("pseudoWH_pre2", Y = as.double(Y), D = as.double(D), time = as.double(time), obsT = as.double(obstimes), status = as.double(status), S_E = double(neval), varS_E = double(neval), Sp = as.double(Sp), N = as.integer(N), NT = as.integer(NT), neval = as.integer(neval), ind = as.integer(ind), varSkMethod = as.integer(varSk.methodN), covSkSl = double(N * (N + 1) / 2))
    }
    else {
      # approximate formula for the variance
      # var.method == "asymptotic (ind.)"
      temp_res <- .C("pseudoWH_asy2", Y = as.double(Y), D = as.double(D), time = as.double(time), obsT = as.double(obstimes), status = as.double(status), S_E = double(neval), varS_E = double(neval), Sp = as.double(Sp), N = as.integer(N), NT = as.integer(NT), neval = as.integer(neval), ind = as.integer(ind))
    }
  }
  else if (method == "lifetable"){
    if (var.method == "none"){
      # estimates only survival, but not the variance
      temp_res <- .C("pseudoWH_none_LT", Y = as.double(Y), D = as.double(D), C = as.double(C), time = as.double(time), obsT = as.double(obstimes), status = as.double(status), S_E = double(neval), Sp = as.double(Sp), N = as.integer(N), NT = as.integer(NT), neval = as.integer(neval), ind = as.integer(ind))
    }
    else if (var.method == "precise"){
      # precise formula for the variance
      temp_res <- .C("pseudoWH_pre2_LT", Y = as.double(Y), D = as.double(D), C = as.double(C), time = as.double(time), obsT = as.double(obstimes), status = as.double(status), S_E = double(neval), varS_E = double(neval), Sp = as.double(Sp), N = as.integer(N), NT = as.integer(NT), neval = as.integer(neval), ind = as.integer(ind), varSkMethod = as.integer(varSk.methodN), covSkSl = double(N * (N + 1) / 2))
    }
    else {
      # approximate formula for the variance
      # var.method == "asymptotic (ind.)"
      temp_res <- .C("pseudoWH_asy2_LT", Y = as.double(Y), D = as.double(D), C = as.double(C), time = as.double(time), obsT = as.double(obstimes), status = as.double(status), S_E = double(neval), varS_E = double(neval), Sp = as.double(Sp), N = as.integer(N), NT = as.integer(NT), neval = as.integer(neval), ind = as.integer(ind))
    }
  }
  
  out <- NULL
  out$time=eval.times
  out$surv=temp_res$S_E
  if (var.method != "none"){out$std.err=sqrt(temp_res$varS_E)}
  rownames(out) <- NULL
  # confidence intervals
  se.fac <- sqrt(qchisq(conf.int, 1))
  if (conf.type == "plain") {
    out$lower <- as.vector(out$surv - se.fac * out$std.err)
    out$upper <- as.vector(out$surv + se.fac * out$std.err)
  }
  else if (conf.type == "log") {
    out$lower <- exp(as.vector(log(out$surv) - se.fac * out$std.err / out$surv))
    out$upper <- exp(as.vector(log(out$surv) + se.fac * out$std.err / out$surv))
  }
  else if (conf.type == "log-log") {
    out$lower <- exp(-exp(as.vector(log(-log(out$surv)) - 
                                      se.fac * out$std.err/out$surv/log(out$surv))))
    out$upper <- exp(-exp(as.vector(log(-log(out$surv)) + 
                                      se.fac * out$std.err/out$surv/log(out$surv))))
  }
  
  out$n <- temp_res$N
  out$conf.type <- conf.type
  out$method <- method
  out$conf.int <- conf.int
  out$call <- call
  # class(out) <- c("survfit","rs.surv")
  rm(temp_res)
  out
}

# # testni primer
# data(colrec)
# td <- colrec[990:999,1:5] # test data
# td$time <- c(70,70,458,475,475,861,980,1196,1360,3652)
# td$stat <- c(1,0,1,1,1,1,0,1,1,0)
# td <- td[order(td$time),]
# td
# 
# casi <- c(50,800,td$time)
# 
# P <- relsurv.rs.survPseudo(Surv(time,stat)~ratetable(age=age,sex=sex,year=diag),data=td,ratetable=slopop, var.method = "pre",eval.times=casi,varSk.method="re")
# P_old <- rs.survPseudo(time = td$time, status = td$stat, age = td$age/365.241, sex = td$sex, year = td$diag, times = casi, var.method = "asy", varSk.method = "pre")


# set.seed(23)
# N <- 100 #precise vresion: 10000 = 30min, 20000=4.4h
# x <- sample(c(0,1),N,repl=T)
# ages <- runif(N,min=50,max=65)
# spol <- sample(c(1,2),N,repl=TRUE)
# myAgedis <- function(n){ages};age_ef <- sex_ef <- 1; l_start=1990; l_stop=1999;lambda_0 <- 1e-4
# data <- gen_net2(n=N,sex=spol,agedis=myAgedis,start=l_start,stop=l_stop,base_haz=function(t){lambda_0},betas=c(0,age_ef/365.241,0,sex_ef),vars=matrix(x,ncol=1))
# # beta1 pripada motecemu dejavniku x, beta2 starosti, beta3 letu in beta4 spolu
# data$t_obs <- pmin(data$t_pop_ev,data$t_ex_ev)
# data$obs_stat <- (data$ex_stat|data$pop_stat)
# data$t_cens <- runif(N,min=0,max=365.241*(l_stop-l_start+1)*4)
# data$cens_stat <- (data$t_cens < 365.241*(l_stop-l_start+1))
# data$t_obs2 <- pmin(data$t_obs,data$t_cens)
# # data$t_obs2 <- round(pmin(data$t_obs,data$t_cens))
# data$obs_stat2 <- (data$t_obs <= data$t_cens)*data$obs_stat
# evalT <- data$t_obs2#[which(data$obs_stat2==1)] # seq(100,3500,100)
# system.time(outC <- relsurv.rs.survPseudo(Surv(t_obs2,obs_stat2)~sex+x+ratetable(age=age*365.241,sex=sex,year=year),data=data,ratetable=slopop, var.method = "pre",eval.times=evalT,varSk.method="pre"))
# # system.time(out <- rs.survPseudo(time = data$t_obs2, status = data$obs_stat2, age = data$age, sex = data$sex, year = data$year, times = evalT, var.method = "p", varSk.method = "re"))
# system.time(PP <- rs.surv(Surv(t_obs2,obs_stat2)~ratetable(age=age*365.241,sex=sex,year=year),data=data,ratetable=slopop,conf.type = "p"))
# # outC$surv
# # summary(PP,times=evalT)$surv
# plot(PP,mark.time=FALSE,ylim=c(0.4,1))
# lines(outC$time,outC$surv,col="red",type="s")
# lines(outC$time,outC$lower,col="blue",type="s",lty=2)
# lines(outC$time,outC$upper,col="blue",type="s",lty=2)

# # ce vezi
# N <- 2000 #precise version: 10000 = 30min, 20000=4.4h
# x <- sample(c(0,1),N,repl=T)
# ages <- runif(N,min=95,max=100)
# spol <- sample(c(1,2),N,repl=TRUE)
# myAgedis <- function(n){ages};age_ef <- sex_ef <- 0.5; l_start=1990; l_stop=1990;lambda_0 <- 1e-2
# data <- gen_net2(n=N,sex=spol,agedis=myAgedis,start=l_start,stop=l_stop,base_haz=function(t){lambda_0},betas=c(0,age_ef/365.241,0,sex_ef),vars=matrix(x,ncol=1))
# # beta1 pripada motecemu dejavniku x, beta2 starosti, beta3 letu in beta4 spolu
# data$t_ex_ev <- round(data$t_ex_ev)
# data$t_pop_ev <- round(data$t_pop_ev)
# data$t_obs <- pmin(data$t_pop_ev,data$t_ex_ev)
# data$obs_stat <- (data$ex_stat|data$pop_stat)
# data$t_cens <- round(runif(N,min=0,max=365.241*(l_stop-l_start+1)*5))
# data$cens_stat <- (data$t_cens < 365.241*(l_stop-l_start+1))
# data$t_obs2 <- pmin(data$t_obs,data$t_cens)
# # data$t_obs2 <- round(pmin(data$t_obs,data$t_cens))
# data$obs_stat2 <- (data$t_obs <= data$t_cens)*data$obs_stat
# barplot(table(data$t_obs2))
# evalT <- data$t_obs2[which(data$obs_stat2==1)] # seq(100,3500,100)
# system.time(outC <- relsurv.rs.survPseudo(Surv(t_obs2,obs_stat2)~ratetable(age=age*365.241,sex=sex,year=year),data=data,ratetable=slopop, var.method = "pre",eval.times=evalT,varSk.method="pre"))
# # system.time(out <- rs.survPseudo(time = data$t_obs2, status = data$obs_stat2, age = data$age, sex = data$sex, year = data$year, times = evalT, var.method = "p", varSk.method = "re"))
# system.time(PP <- rs.surv(Surv(t_obs2,obs_stat2)~ratetable(age=age*365.241,sex=sex,year=year),data=data,ratetable=slopop,conf.type = "p"))
# # outC$surv
# # summary(PP,times=evalT)$surv
# plot(PP,mark.time=FALSE,ylim=c(0,1))
# lines(outC$time,outC$surv,col="red",type="s")
# lines(outC$time,outC$lower,col="red",type="s",lty=2)
# lines(outC$time,outC$upper,col="red",type="s",lty=2)
# P_lt <- relsurv.rs.survPseudo(Surv(t_obs2,obs_stat2)~ratetable(age=age*365.241,sex=sex,year=year),data=data,ratetable=slopop, var.method = "pre",eval.times=evalT,varSk.method="pre",method="l")
# # P_lt <- pseudoLifeTabEst(data$t_obs2,data$obs_stat2,data$age,data$year,data$sex,evalT)
# points(P_lt$time,P_lt$surv,pch=16,col="blue",cex=0.3)
# points(P_lt$time,P_lt$lower,pch=16,col="blue",cex=0.3)
# points(P_lt$time,P_lt$upper,pch=16,col="blue",cex=0.3)
