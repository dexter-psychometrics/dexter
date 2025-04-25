
library(dplyr)


RcppArmadillo::armadillo_throttle_cores(1)

verbAggCopy = function(pth = test_path('testdata/verbAggression.db'))
{
  con = DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  db = open_project(pth)
  
  RSQLite::sqliteCopyDatabase(db, con)
  
  DBI::dbDisconnect(db)
  return(con)
}

test_that('populations work',{
  db = verbAggCopy()
  
  #artificially create two overlapping booklets
  f2 = fit_enorm(db, (gender=='Male' & item_position<=16) | 
                     (gender=='Female' & item_position>8))
  
  expect_true(n_distinct(f2$inputs$design$booklet_id)==2)
  
  set.seed(123)
  pv=plausible_values(db, f2, covariates='gender')
  
  expect_true(mean(pv[pv$gender=='Male',]$PV1) > mean(pv[pv$gender=='Female',]$PV1))
  
  
  
  abl = ability(db,f2,method='WLE')
  
  test = inner_join(abl,pv,by='person_id')
  
  expect_gt(cor(test$PV1,test$theta), .9,'verb agg correlation ability and plausible value should be larger than .9')
  
  
  # see if sanity checks work
  p=get_persons(db) |>
    select(person_id) |>
    mutate(x1=rnorm(n()),x2 = sample(rnorm(3),n(),replace=TRUE))
  
  add_person_properties(db,p)
  
  expect_warning(plausible_values(db, f2, predicate=startsWith(item_id,'S1'), covariates='x1'),regexp='ignoring covariates',ignore.case=TRUE)
  expect_warning(plausible_values(db, f2, predicate=startsWith(item_id,'S1'), covariates='x2'),regexp='decimal',ignore.case=TRUE)
  #to do: what happens if one group? SHould we have a message or not?
  
  dbDisconnect(db)
})
  
RcppArmadillo::armadillo_reset_cores()