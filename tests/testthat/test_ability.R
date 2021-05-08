context('check ability')

library(dplyr)

expect_no_error = function(object, info=NULL) expect_error(object, regexp=NA, info=info)



test_that('inconsistencies between data and parms are handled correctly',{

  db = open_project('../verbAggression.db')
  
  f1 = fit_enorm(db)
  f2 = fit_enorm(db, item_id != 'S4DoShout')
  f3 = fit_enorm(db, !(item_id=='S4DoShout' & item_score == 1))
  
  # params must cover all item, item_score combinations 
  expect_no_error({p1 = ability(db,f1)})
  expect_error({p2 = ability(db,f2)}, regexp='parameters.+items')
  expect_error({p3 = ability(db,f3)}, regexp='parameters.+scores')
  
  # of course the reverse is not necessary
  expect_no_error({p4 = ability(db,f1,item_score !=2 )})
  
  dbDisconnect(db)
})  


test_that('verbAgg abilities are monotone increasing', {

  db = open_project('../verbAggression.db')
  f = fit_enorm(db)
  
  # check ability mle is correct
  es = expected_score(f)
  expect_lt(
    ability_tables(f) %>%
      filter(is.finite(theta)) %>%
      mutate(error = abs(booklet_score - es(theta))) %>%
      pull(error) %>%
      mean(),
    0.001,
    label = "ability_tables mle on average estimated to within .001 of test_score")
  
  
  nscores = get_rules(db) %>%
    group_by(item_id) %>%
    summarize(m=max(item_score)) %>%
    ungroup() %>%
    pull(m) %>%
    sum() + 1
  
  for(method in eval(formals(ability_tables)$method))
  {
    for(prior in eval(formals(ability_tables)$prior))
      {
      expect_no_error({abl = ability_tables(f, method = method)}, info = paste('fit_enorm verbAgg -', method))
      expect_true(nrow(abl) == nscores)
      
      abl = abl %>%  
        filter(is.finite(theta)) %>%
        arrange(booklet_score) %>%
        mutate(p = lag(theta, default = -Inf)) %>%
        filter(theta < p)
  
      expect_true(nrow(abl) == 0, info = paste('abilities not increasing verbAgg -',method))  
    }
  }

  dbDisconnect(db)
})






