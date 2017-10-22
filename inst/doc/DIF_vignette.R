## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,fig.align='centre', fig.width=6, fig.height=5)
library(dexter)

## ----get_data, results='hide'--------------------------------------------
db = start_new_project(verbAggrRules, "verbAggression.db", covariates=list(gender=""))
add_booklet(db, verbAggrData, "data")
add_item_properties(db, verbAggrProperties)

## ------------------------------------------------------------------------
dif_gender=DIF(db,"gender")
print(dif_gender)

## ---- echo=FALSE, results='hide'-----------------------------------------
cc=fit_enorm(db,(gender=="Male")&(item_id%in%c("S1DoCurse","S1WantCurse","S3WantScold")))
plot(-2:3,c(rep(0,3),rep(1,3)),xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='ability-scale')
lines(-1.5:2.5,rep(0.2,5),lty=2,col="gray")
lines(-1.5:2.5,rep(0.8,5),lty=2,col="gray")
for (i in 1:6)
text(cc$est$beta.cml[i],0.8,rownames(cc$est$beta.cml)[i],cex=0.6,adj=1,srt=90,pos = 3)

cc=fit_enorm(db,(gender=="Female")&(item_id%in%c("S1DoCurse","S1WantCurse","S3WantScold")))
for (i in 1:6)
text(cc$est$beta.cml[i],0.2,rownames(cc$est$beta.cml)[i],cex=0.6,adj=1,srt=90,pos = 3)
text(-1.9,0.8,"Males")
text(-1.87,0.2,"Females")

## ------------------------------------------------------------------------
plot(dif_gender)

## ------------------------------------------------------------------------
D=dif_gender$DIF_pair
o=c(grep("Do",rownames(D)),grep("Want",rownames(D)))
dif_gender$DIF_pair=D[o,o]
plot(dif_gender)

## ---- echo=FALSE, results='hide'-----------------------------------------
dbDisconnect(db)

