## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(dexter)
library(dplyr)

## ------------------------------------------------------------------------
sim_Rasch = function(theta, delta) {
  n = length(theta)
  m = length(delta)
  data.frame(
    person_id = rep(paste0('p',1:n), m),
    item_id = rep(paste0('i',1:m), each=n),
    item_score = as.integer(rlogis(n*m, outer(theta, delta, "-")) > 0)
  )
}

simulated = sim_Rasch(rep(0.5, 2000), runif(20, -2, 2))

## ---- fig.align='center', fig.width=7------------------------------------
ss= simulated %>% 
  group_by(person_id) %>% 
  summarise(sumscore=sum(item_score)) 
par(mfrow=c(1,2))
hist(ss$sumscore)
plot(ecdf(ss$sumscore))
mm = fit_inter(simulated)

## ---- fig.align='center', fig.width=7------------------------------------
mm = fit_inter(simulated)
par(mfrow=c(1,1))
plot(mm, show.observed = TRUE, 
     items = c('i1','i2'),
     nr=1, nc=2)

## ---- fig.align='center', fig.height=4, fig.width=4----------------------
dd = individual_differences(simulated)
plot(dd)

## ------------------------------------------------------------------------
print(dd)

## ---- fig.align='center', fig.height=4, fig.width=4----------------------
db2 = start_new_project(verbAggrRules, "verbAggression.db")
add_booklet(db2, verbAggrData, "data")
dd = individual_differences(db2, booklet_id=="data")
plot(dd)

## ---- show=FALSE---------------------------------------------------------
dbDisconnect(db2)

