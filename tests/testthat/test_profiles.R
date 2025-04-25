
library(dplyr)


test_that('profile analysis verb agg',{

  db = open_project(test_path('testdata/verbAggression.db'))
  
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
  
  r = get_responses(db,columns=c('person_id','item_id','situation','item_score')) |>
    mutate(p=dense_rank(person_id)%%2) |>
    filter(!(situation=='Call' & p==1))
  
  p = profiles(r,f,'situation') |>
    mutate(p=dense_rank(person_id)%%2) |>
    count(p,situation)
  
  expect_true(n_distinct(p$n)==1 && sum(p$situation=='Call')==1 && nrow(p)==7,
              label='profiles, unequal categories correctly handled')
  

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
