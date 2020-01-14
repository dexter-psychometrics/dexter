context('Check fit_enorm')

library(dplyr)

expect_no_error = function(object) expect_error(object, regexp=NA)




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
  
  items = get_items(db) %>%
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
            mutate(difference=abs(SE_beta-se.b)))$difference),
    1e-3)
  
  #check that Bayesian is reasonably close.
  fxb = fit_enorm(db, fixed_params=oplm_params %>% filter(is_fixed), method='Bayes')
  
  # this is a simulation based test so result differ each time
  # but typical value between 0.01-0.015, never observed higher than 0.05
  expect_lt(coef(fxb) %>%
              inner_join(oplm_params, by=c('item_id','item_score')) %>%
               mutate(difference = abs(mean_beta - beta)) %>%
               filter(!is_fixed) %>%
               pull(difference) %>%
               mean(),
             0.1)
  
  # check that omitting some score category gives correct output and an error
  # mentioning the correct item
  expect_output(   
    expect_error(fit_enorm(db, 
                   fixed_params=oplm_params %>% 
                     filter(is_fixed & !(item_id=='S1DoShout' & item_score == 1))),
    regexp='missing.+categories.+fixed'),
  regexp='S1DoShout')
  
  # check for correct error message
  # no score variation
  expect_error(fit_enorm(db, item_score == 1), regexp='score.+variation')
  # no 0 score category
  expect_output(
    expect_error(fit_enorm(db, item_score > 0 | item_id!='S1DoShout'), regexp='minimum.+must.+(zero)|0', ignore.case=TRUE),
    regexp='S1DoShout')
  expect_no_error(fit_enorm(db, item_score <= 1))
  
  
  
  dbDisconnect(db)
})

