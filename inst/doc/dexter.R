## ------------------------------------------------------------------------
library(dexter)
head(verbAggrRules)
db = start_new_project(verbAggrRules, "verbAggression.db", covariates=list(gender="<unknown>"))

## ------------------------------------------------------------------------
head(verbAggrData)
add_booklet(db, verbAggrData, "agg")

## ------------------------------------------------------------------------
show_booklets(db)
show_items(db)

## ------------------------------------------------------------------------
data("verbAggrProperties")
head(verbAggrProperties)
add_item_properties(db, verbAggrProperties)
show_item_properties(db)
show_person_properties(db)

## ------------------------------------------------------------------------
tt = tia_tables(db)
knitr::kable(tt$itemStats, digits=3)
knitr::kable(tt$testStats, digits=3)

## ------------------------------------------------------------------------
m = fit_inter(db, booklet_id=='agg')
print(m)

## ------------------------------------------------------------------------
plot(m, "S1DoScold", show.observed=TRUE)

## ------------------------------------------------------------------------
plot(m, 'S1DoCurse', summate=FALSE)

## ------------------------------------------------------------------------
mSit = fit_domains(db, item_property= "situation")
plot(mSit)

## ------------------------------------------------------------------------
profile_plot(db, item_property='mode', covariate='gender')

## ------------------------------------------------------------------------
parms = fit_enorm(db)

## ------------------------------------------------------------------------
parms_gibbs = fit_enorm(db, method='Bayes')

## ------------------------------------------------------------------------
pv = plausible_values(db, parms)
plot(density(pv$PV1))

## ------------------------------------------------------------------------
pv = merge(pv, dbReadTable(db, 'dxPersons'))
boxplot(PV1~gender, data=pv)

## ---- show=FALSE---------------------------------------------------------
dbDisconnect(db)

