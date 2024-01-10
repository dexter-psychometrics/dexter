context('check equating')

library(dplyr)
library(tidyr)

RcppArmadillo::armadillo_throttle_cores(1)

test_that('complete case works',{
  
  db = open_project('../verbAggression.db')

  ref_items = get_items(db) |> 
    filter(behavior == 'Scold') |> 
    pull(item_id)
  
  f = fit_enorm(db,method='Bayes')
  
  p2p = probability_to_pass(db,f,ref_items,pass_fail=7)
  
  expect_true(!anyNA(coef(p2p)),'no NA values')
  expect_true(all(between(as.matrix(select(coef(p2p),-(1:2))),0,1)), 'coef contains probabilities')
  
  pass = get_testscores(db, behavior == 'Scold') |>
    mutate(pass=booklet_score>=7) |>
    select(person_id,pass) |>
    inner_join(get_testscores(db)) |>
    group_by(booklet_score) |>
    summarise(mean_pass=mean(pass))
  
  tst = inner_join(pass, coef(p2p), by=c('booklet_score'='score_new'))
  
  expect_gt(cor(tst$probability_to_pass,tst$mean_pass),.9, 'relation between manifest prob pass and model prob pass')
  expect_lt(mean(abs(tst$probability_to_pass-tst$mean_pass)),.1, 'manifest prob pass estimated reasonably accurately')
  

})

RcppArmadillo::armadillo_reset_cores()