## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width=6, fig.height=5)
library(dexter)

## ----get_data, results='hide'--------------------------------------------
db = start_new_project(verbAggrRules, "verbAggression.db", covariates=list(gender=""))
add_booklet(db, verbAggrData, "data")
add_item_properties(db, verbAggrProperties)

## ------------------------------------------------------------------------
dif_gender=DIF(db,"gender")
print(dif_gender)

## ---- echo=FALSE, results='hide', fig.height=3---------------------------
cc=fit_enorm(db,(gender=="Male")&(item_id%in%c("S1DoCurse","S1WantCurse","S3WantScold")))
plot(c(-2,3),c(0,1),xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='ability-scale')
lines(-1.4:2.6,rep(0,5),lty=2,col="gray")
lines(-1.4:2.6,rep(0.8,5),lty=2,col="gray")
for (i in 1:6)
text(cc$est$beta.cml[i],0.8,rownames(cc$est$beta.cml)[i],cex=0.6,adj=1,srt=90,pos = 3, xpd=NA)

cc=fit_enorm(db,(gender=="Female")&(item_id%in%c("S1DoCurse","S1WantCurse","S3WantScold")))
for (i in 1:6)
text(cc$est$beta.cml[i],0,rownames(cc$est$beta.cml)[i],cex=0.6,adj=1,srt=90,pos = 3, xpd=NA)
text(-1.3,0.8,"Males", pos=2, xpd=NA)
text(-1.3,0,"Females", pos=2, xpd=NA)

## ------------------------------------------------------------------------
plot(dif_gender)

## ------------------------------------------------------------------------
D=dif_gender$DIF_pair
o=c(grep("Do",rownames(D)),grep("Want",rownames(D)))
dif_gender$DIF_pair=D[o,o]
plot(dif_gender)

## ---- echo=FALSE, results='hide'-----------------------------------------
close_project(db)

