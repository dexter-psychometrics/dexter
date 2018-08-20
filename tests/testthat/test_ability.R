context('check ability')

library(dplyr)
library(purrr)

expect_no_error = function(object, info=NULL) expect_error(object, regexp=NA, info=info)



test_that('inconsistencies between data and parms are handled correctly',{
  db = open_project('../verbAggression.db')
  
  f1 = fit_enorm(db)
  f2 = fit_enorm(db, item_id != 'S4DoShout')
  f3 = fit_enorm(db, !(item_id=='S4DoShout' & item_score == 1))
  
  # params must cover all item, item_score combinations 
  expect_no_error({p1 = ability(db,f1)})
  expect_error({p2 = ability(db,f2)}, regexp='items.+parameters')
  expect_error({p3 = ability(db,f3)}, regexp='item_scores.+parameters')
  
  # of course the reverse is not necessary
  expect_no_error({p4 = ability(db,f1,item_score !=2 )})
  
  dbDisconnect(db)
})  


test_that('verbAgg abilities are monotone increasing', {
  db = open_project('../verbAggression.db')
  f = fit_enorm(db)
  
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
        arrange(sumScore) %>%
        mutate(p = lag(theta, default = -Inf)) %>%
        filter(theta < p)
  
      expect_true(nrow(abl) == 0, info = paste('abilities not increasing verbAgg -',method))  
    }
  }
  dbDisconnect(db)
})


## tests on big datasets
test_that('different estimation methods of ability converge on a large dataset',{
  skip_on_cran()
  if(dir.exists('../skip_on_cran'))
  {
    db = open_project('../skip_on_cran/pisa/pisa_math.db')
    
    expect_no_error({f = fit_enorm(db, booklet_id %in% c('1','3','4'))}, info='pisa fit_enorm')

    meth = eval(formals(ability_tables)$method)
    prior = eval(formals(ability_tables)$prior)
    
    # run ability for each available method
    est = reduce(apply(expand.grid(meth, prior),1,function(m)
    {
      expect_no_error({abl = ability_tables(f, method = m[1], prior=m[2])}, info=paste('ability pisa -', m[1], m[2]))
      colnames(abl)[3:ncol(abl)] = paste(m[1],m[2],colnames(abl)[3:ncol(abl)],sep='.')
      abl
    }), 
    inner_join, by=c('booklet_id','sumScore'))
    
    # inspect a regression of the theta's
    combn(paste0(apply(expand.grid(meth, prior),1,paste0,collapse='.'),'.theta'), 2, function(pair)
    {
      test = est %>% 
        select(one_of(pair)) %>%
        filter_all( all_vars(is.finite(.)))
      # some are exactly equal
      hh = suppressWarnings(summary(lm(test[[1]] ~ test[[2]])))
      
      expect_true(hh$r.squared > 0.98, 
                  info = paste('lm theta ',pair[1], '~', pair[2], paste0('expect r.squared > .99, found:',hh$r.squared )))


    })

    dbDisconnect(db)
  }
})


