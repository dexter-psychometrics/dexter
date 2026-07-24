

RcppArmadillo::armadillo_throttle_cores(1)


test_that("pv's work",{
  db = verbAggCopy()
  
  #artificially create two overlapping booklets
  f2 = fit_enorm(db, (gender=='Male' & item_position<=16) | 
                     (gender=='Female' & item_position>8))
  
  expect_true(n_distinct(f2$inputs$design$booklet_id)==2)
  
  # this seed works with both 1 (mac) and 2 cores. 
  # Inelegant to use a seed but necessary, statistically this test only succeeds in ~94% of cases, which is entirely correct
  # to do: change for the next version since I do not want tot test everything with 1 and 2 cores but I do want some pv tests on cran
  set.seed(723)
  pv = plausible_values(db, f2, covariates='gender',nPV=10)
  
  expect_true(mean(pv[pv$gender=='Male',]$PV1) > mean(pv[pv$gender=='Female',]$PV1))
  
  expect_true(df_join_equal(get_testscores(db), select(pv,'person_id','booklet_id','booklet_score'), join_by='person_id'))
  
  
  pv_nocovar = plausible_values(db, f2,nPV=10) |>
    inner_join(get_persons(db), by='person_id')
  
  diff_nocovar = pv_nocovar |>
    pivot_longer(starts_with('PV'), names_to='iter', values_to='pv') |>
    group_by(gender) |>
    summarise(mu=mean(pv)) |>
    summarise(d=abs(diff(mu))) |>
    pull(d)
  
  diff_covar = pv |>
    pivot_longer(starts_with('PV'), names_to='iter', values_to='pv') |>
    group_by(gender) |>
    summarise(mu=mean(pv)) |>
    summarise(d=abs(diff(mu))) |>
    pull(d)
  
  # difference should be larger if we use a relevant covariate
  expect_gt(diff_covar, diff_nocovar)
  
  
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