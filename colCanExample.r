# pseudoPaper example
library(relsurv)
data("colrec")
data("slopop")
source("/home/klemen/work/Pseudovalues/Simulations/relsurv.rs.survPseudo.r")
head(colrec)
d_col <- colrec[colrec$site=="colon",]
dim(colrec)
dim(d_col)
table(d_col$site)
d_col$site <- factor(d_col$site)
levels(d_col$site)
head(d_col)
dim(d_col)
d_col <- d_col[d_col$sex==2,]
summary(d_col$age/365)
d_col$agegr <- cut(d_col$age/365.241,seq(10,100,5))
head(d_col)
nessie(Surv(time,stat)~agegr+ratetable(age=age,sex=sex,year=diag),ratetable=slopop,data=d_col,times=c(5,10,15,20))
# d_col <- d_col[which(d_col$age/365 <= 75 & d_col$age/365 >= 20),]
d_col <- d_col[which(d_col$age/365 <= 85 & d_col$age/365 >= 20),]
nessie(Surv(time,stat)~agegr+ratetable(age=age,sex=sex,year=diag),ratetable=slopop,data=d_col,times=c(5,10,15,20))
dim(d_col)
summary(d_col$age/365)
summary(d_col$time)
summary(as.date(d_col$diag))
table(d_col$sex)


PP<- rs.surv(Surv(time,stat)~1+ratetable(sex=sex,age=age,year=diag),ratetable=slopop,data=d_col,conf.type = "p")

P <- relsurv.rs.survPseudo(Surv(time,stat)~1+ratetable(sex=sex,age=age,year=diag),ratetable=slopop,data=d_col,var.method = "p",varSk.method = "p")

plot(PP,mark.time=F,ylab="",xlab="Time (years)",ylim=c(0.3,1),lwd=4.5,cex.axis=1.5,cex.lab=1.5,xaxt="n",xlim=c(0,20*365.241))
x_os <- seq(0,20,5)
axis(1,at=x_os*365,labels = x_os,cex.axis=1.5)
lines(P$time,P$surv,col="grey",type="s",lwd=2.5,lty=1)
lines(P$time,P$lower,col="grey",type="s",lwd=3,lty=2)
lines(P$time,P$upper,col="grey",type="s",lwd=3,lty=2)
legend("topright",c("PP","new"),fill=c("black","grey"),cex=1.5,bty="n")

d_col2 <- d_col[which(d_col$age / 365 <= 70),] 
nrow(d_col2)
PP2<- rs.surv(Surv(time,stat)~1+ratetable(sex=sex,age=age,year=diag),ratetable=slopop,data=d_col2,conf.type = "p")

P2 <- relsurv.rs.survPseudo(Surv(time,stat)~1+ratetable(sex=sex,age=age,year=diag),ratetable=slopop,data=d_col2,var.method = "p",varSk.method = "p")

plot(PP2,mark.time=F,ylab="",xlab="Time (years)",ylim=c(0.3,1),lwd=4.5,cex.axis=1.5,cex.lab=1.5,xaxt="n",xlim=c(0,20*365.241))
x_os <- seq(0,20,5)
axis(1,at=x_os*365,labels = x_os,cex.axis=1.5)
lines(P2$time,P2$surv,col="grey",type="s",lwd=2.5,lty=1)
lines(P2$time,P2$lower,col="grey",type="s",lwd=3,lty=2)
lines(P2$time,P2$upper,col="grey",type="s",lwd=3,lty=2)
legend("topright",c("PP","new"),fill=c("black","grey"),cex=1.5,bty="n")

par(mfrow=c(1,2))
# na vseh
plot(PP,mark.time=F,ylab="",xlab="Time (in years)",ylim=c(0.3,1),lwd=4.5,cex.axis=1.5,cex.lab=1.5,xaxt="n",xlim=c(0,20*365.241))
x_os <- seq(0,20,5)
axis(1,at=x_os*365,labels = x_os,cex.axis=1.5)
lines(P$time,P$surv,col="grey",type="s",lwd=2.5,lty=1)
lines(P$time,P$lower,col="grey",type="s",lwd=3,lty=2)
lines(P$time,P$upper,col="grey",type="s",lwd=3,lty=2)
legend("topright",c("PP","new"),fill=c("black","grey"),cex=1.5,bty="n")
# na primernih
plot(PP2,mark.time=F,ylab="",xlab="Time (in years)",ylim=c(0.3,1),lwd=4.5,cex.axis=1.5,cex.lab=1.5,xaxt="n",xlim=c(0,20*365.241))
x_os <- seq(0,20,5)
axis(1,at=x_os*365,labels = x_os,cex.axis=1.5)
lines(P2$time,P2$surv,col="grey",type="s",lwd=2.5,lty=1)
lines(P2$time,P2$lower,col="grey",type="s",lwd=3,lty=2)
lines(P2$time,P2$upper,col="grey",type="s",lwd=3,lty=2)
legend("topright",c("PP","new"),fill=c("black","grey"),cex=1.5,bty="n")
par(mfrow=c(1,1))

utezi <- relsurv:::exp.prep(data.frame(age=d_col$age,year=d_col$diag,sex=d_col$sex),20*365.241,slopop)

inx_min <- which(utezi == min(utezi))
# 
# data.frame(age=d_col$age,year=d_col$diag,sex=d_col$sex)[inx_min,1] / 365.241#age
# as.date(data.frame(age=d_col$age,year=d_col$diag,sex=d_col$sex)[inx_min,2]) #year
utezi[inx_min]
# 1 / utezi[inx_min]
# summary(1/utezi)
# hist(1/utezi)
# slopop["89","1996",2]
# 
# for (l in 1:20){
#   print(l)
#   print(1/relsurv:::exp.prep(data.frame(age=d_col$age,year=d_col$diag,sex=d_col$sex)[inx_min,],l*365.241,slopop));
#   print(-log(relsurv:::exp.prep(data.frame(age=d_col$age,year=d_col$diag,sex=d_col$sex)[inx_min,],l*365.241,slopop)))
# }

# utezi po 17 letih za tiste, ki so se at risk at that time
inx_w <- which(d_col$time >= 17*365.241)
length(inx_w)
# sort(d_col$age[inx_w]/365.241)
weigh <- 1/relsurv:::exp.prep(data.frame(age=d_col$age[inx_w],year=d_col$diag[inx_w],sex=d_col$sex[inx_w]),17*365.241,slopop)
sort(weigh)
summary(weigh)

# Za SdS za plakat
par(mar=c(4,4.2,0,1)+0.1)
plot(PP,mark.time=F,ylab="Čisto preživetje",xlab="Čas (v letih)",ylim=c(0.3,1),lwd=4.5,cex.axis=1.3,cex.lab=1.3,xaxt="n",xlim=c(0,20*365.241),conf.int=F)
x_os <- seq(0,20,5)
axis(1,at=x_os*365,labels = x_os,cex.axis=1.3)
lines(P$time,P$surv,col="blue",type="s",lwd=2.5,lty=1)
legend("topleft",c("cenilka PP","nova cenilka"),fill=c("black","blue"),cex=1.5,bty="n")
par(mar=c(5,4,4,2)+0.1)
