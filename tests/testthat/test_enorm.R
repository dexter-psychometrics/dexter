library(dplyr)
library(purrr)

# better names
expect_no_error = function(object) expect_error(object, regexp=NA)


context('Check fit_enorm')

test_that('calibration of verbal aggression dataset matches oplm results, with fixed and unfixed',{

  db = open_project('../verbAggression.db')
  
  #free calibration
  ff = fit_enorm(db)
  
  # check free calibration
  # mean absolute difference between oplm and dexter parameters should be negligible
  # but will not be 0 because of platforms, rounding, floating point, etc.
  expect_lt(
      mean((coef(ff) %>%
             mutate(item_id=substr(item_id,1,8)) %>%
             inner_join(read_oplm_par('../verbal_oplm/VERBAL.PAR'), by=c('item_id','item_score')) %>%
             mutate(difference=abs(beta.x-beta.y)))$difference),
      1e-15)
  
  #calibration with fixed_parameters
  
  # determine which are fixed from the cml file but use the parameters from the par file
  # since they are less rounded
  oplm_params = read_oplm_par('../verbal_oplm/VERBAL_FX.PAR') %>%
    inner_join(read_oplm_par('../verbal_oplm/VERBAL_FX.CML'), by=c('item_id','item_score')) %>%
    mutate(is_fixed=is.na(se.b)) %>%
    select(oplm_lab=item_id, item_score, beta=beta.x, se.b, is_fixed)
  
  # fixing id's etc.
  oplm_params[is.na(oplm_params$se.b),'se.b'] = 0
  
  items = show_items(db) %>%
    mutate(oplm_lab = substr(item_id,1,8))
  
  oplm_params = oplm_params %>% inner_join(items, by='oplm_lab')
  
  # calibration should not give errors or warnings
  fx = fit_enorm(db, fixed_params=oplm_params %>% filter(is_fixed))
  
  # beta correct
  expect_lt(
    mean((coef(fx) %>%
            inner_join(oplm_params, by=c('item_id','item_score')) %>%
            mutate(difference=abs(beta.x-beta.y)))$difference),
    1e-15)
  
  # se_b correct, less strict because se_b from cml file is severely rounded
  expect_lt(
    mean((coef(fx) %>%
            inner_join(oplm_params, by=c('item_id','item_score')) %>%
            mutate(difference=abs(SE_b-se.b)))$difference),
    1e-3)
  
  # # check that Bayesian is reasonably close.
  fx = fit_enorm(db, fixed_params=oplm_params %>% filter(is_fixed), method='Bayes')
  
  # not yet gotten a proper value from timo
  # expect_lt(as.data.frame(fx) %>%
  #            inner_join(oplm_params, by=c('item_id','item_score')) %>%
  #             mutate(difference = abs(mean_beta - beta)) %>%
  #             pull(difference) %>%
  #             mean(),
  #           1e-15)
  
  # check that omitting some score category gives correct output and an error
  # mentioning the correct item
  expect_output(   
    expect_error(fit_enorm(db, 
                   fixed_params=oplm_params %>% 
                     filter(is_fixed & !(item_id=='S1DoShout' & item_score == 1))),
    regexp='missing.+categories.+fixed'),
  regexp='S1DoShout')
  
  # check for correct error message
  expect_error(fit_enorm(db, item_score == 1), regexp='score.+variation')
  expect_no_error(fit_enorm(db, item_score <= 1))
  
  
  
  dbDisconnect(db)
})

