context('Test profile analysis')

library(dplyr)
library(DBI)
library(RSQLite)

RcppArmadillo::armadillo_throttle_cores(1)

verbAggCopy = function(pth = '../verbAggression.db')
{
  con = dbConnect(SQLite(), ":memory:")
  db = open_project(pth)
  
  sqliteCopyDatabase(db, con)
  
  dbDisconnect(db)
  return(con)
}


# to do: check for proper number of rows

test_that('profile analysis verb agg',{

  db = verbAggCopy()
  
  f = fit_enorm(db)
  p = profiles(db, f, 'behavior')
  
  expect_gt(cor(p$domain_score,p$expected_domain_score), 0.6, 
            'expected score should have a relation with observed score')
  
  expect_gt(cor(p$domain_score,p$expected_domain_score), cor(p$booklet_score,p$expected_domain_score),
            'domain should add extra information')
  
  expect_true(all(p |> 
                    group_by(person_id) |> 
                    summarise(sum_dif = abs(sum(expected_domain_score) - first(booklet_score))) |>
                    ungroup() |>
                    pull(sum_dif) < 1e-10), 
              'expected domains scores need to sum to total test score')
  
  # check inputs work with just parms
  pt = profile_tables(f, get_items(db),'situation')
  

  f = fit_enorm(db, method='Bayes')
  p = profiles(db, f, 'behavior')
  
  expect_gt(cor(p$domain_score,p$expected_domain_score), 0.6, 
            'expected score should have a relation with observed score (Bayes)')
  
  
  expect_gt(cor(p$domain_score,p$expected_domain_score), cor(p$booklet_score,p$expected_domain_score),
            'domain should add extra information (Bayes)')
  
  expect_true(all(p |> 
                    group_by(person_id) |> 
                    summarise(sum_dif = abs(sum(expected_domain_score) - first(booklet_score))) |>
                    ungroup() |>
                    pull(sum_dif) < 1e-10), 
              'expected domains scores need to sum to total test score (Bayes)')
  
  close_project(db)
})


RcppArmadillo::armadillo_reset_cores()
