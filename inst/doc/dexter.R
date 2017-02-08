## ------------------------------------------------------------------------
library(dexter)
head(verbAggrRules)
db = start_new_project(verbAggrRules, "verbAggression.db")

## ------------------------------------------------------------------------
head(verbAggrData)
add_booklet(db, verbAggrData, "data")

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
m = fit_models(db, 1)
print(m)

## ------------------------------------------------------------------------
plot(m, 1, show.observed=TRUE)

## ------------------------------------------------------------------------
plot(m, 1, summate=FALSE)

## ------------------------------------------------------------------------
mSit = fit_domains(db, 1, "situation")
plot(mSit)

## ------------------------------------------------------------------------
profile_plot(db, 1, 'mode', 'Gender')

