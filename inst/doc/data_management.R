## ------------------------------------------------------------------------
library(dexter)
db = start_new_project(verbAggrRules, "verbAggression.db", 
                       covariates = list(gender="<unknown>"))

## ------------------------------------------------------------------------
add_booklet(db, x=verbAggrData, booklet_id="agg")

## ------------------------------------------------------------------------
add_item_properties(db, verbAggrProperties)

## ------------------------------------------------------------------------
knitr::kable(get_variables(db))

## ---- eval=FALSE---------------------------------------------------------
#  par = fit_enorm(db, gender=='female' & !(booklet_id == 'pretest' & item_position == 3))

## ---- eval=FALSE---------------------------------------------------------
#  bkl = 'pretest'
#  par = fit_enorm(db, gender=='female' & !(booklet_id == bkl & item_position == 3))

## ---- eval=FALSE---------------------------------------------------------
#  booklet_id = 'pretest' # local variable
#  par = fit_enorm(db, gender=='female' & !(booklet_id == local(booklet_id) & item_position == 3))

## ---- eval=FALSE---------------------------------------------------------
#  # assuming an item property called `cefr_level` exists in the project
#  design = design_as_network(db, booklet_id %in% c('bookletA','bookletX','bookletY') & cefr_level == 'B1')
#  design_is_connected(design)
#  ## [1] TRUE

## ---- eval=FALSE---------------------------------------------------------
#  par = fit_enorm(db, response != 'NA')

## ---- eval=FALSE---------------------------------------------------------
#  # goal: fit the extended nominal response model using only persons without any missing responses
#  library(dplyr)
#  
#  # to select on an aggregate level we extract the data from the Dexter project database
#  # and manipulate it using dplyr
#  data = get_responses(db, columns=c('person_id','item_id','item_score','response')) %>%
#      group_by(person_id) %>%
#      mutate(any_missing = any(response == 'NA')) %>%
#      filter(!any_missing)
#  
#  # the manipulated data can be fed back to the analysis function
#  par = fit_enorm(data)

## ---- show=FALSE---------------------------------------------------------
dbDisconnect(db)

