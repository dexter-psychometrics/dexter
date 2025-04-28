
library(dplyr)


RcppArmadillo::armadillo_throttle_cores(1)


test_that('populations work',{
  db = verbAggCopy()
  
  #artificially create two overlapping booklets
  f2 = fit_enorm(db, (gender=='Male' & item_position<=16) | 
                     (gender=='Female' & item_position>8))
  
  expect_true(n_distinct(f2$inputs$design$booklet_id)==2)
  
  set.seed(123)
  pv=plausible_values(db, f2, covariates='gender')
  
  expect_true(mean(pv[pv$gender=='Male',]$PV1) > mean(pv[pv$gender=='Female',]$PV1))
  
  expect_true(df_join_equal(get_testscores(db), select(pv,'person_id','booklet_id','booklet_score'), join_by='person_id'))
  
  # see that designs are not mangled
  x = get_responses(db, predicate=(gender=='Male' & item_position>=16) | 
                                  (gender=='Female' & item_position<8), 
    columns=c('person_id','item_id','item_score','gender')) |>
    rename(booklet_id='gender')
  
  
  pv2 = plausible_values(x, f2)
  
  expect_true(df_join_equal(get_testscores(x), select(pv2,'person_id','booklet_id','booklet_score'), join_by='person_id'))
  expect_true(mean(pv2[pv2$booklet_id=='Male',]$PV1) > mean(pv2[pv2$booklet_id=='Female',]$PV1))
  
  
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

  dbDisconnect(db)
})
  
RcppArmadillo::armadillo_reset_cores()