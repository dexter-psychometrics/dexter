library(dplyr)


context('test data selection')

expect_no_error = function(object, info=NULL) expect_error(object, regexp=NA, info=info)


test_that('predicates work as expected',
{
  options(dexter.debug = TRUE)

  
  db = open_project('../verbAggression.db')
  
  #two ways to do the same
  
  r1 = get_responses(db, item_id %like% 'S1%')
  
  expect_false('get_responses_error_caught' %in% dexter:::debug.log$retrieve(), info='predicate is do-able with sql')
  
  r2 = get_responses(db, grepl('S1', item_id))
  
  expect_true(dexter:::debug.log$retrieve()$get_responses_filter_success, info='succesfully fall back and evaluate predicate with dplyr filter')
  
  expect_true(dexter:::df_identical(r1, r2))
  
  

  options(dexter.debug = FALSE)  
})
   