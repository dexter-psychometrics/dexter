library(dplyr)


context('test data selection')

expect_no_error = function(object, info=NULL) expect_error(object, regexp=NA, info=info)


test_that('predicates work as expected',
{
  
  db = open_project('../verbAggression.db')
  
  #two ways to do the same
  
  r1 = get_responses(db, item_id %like% 'S1%')

  r2 = get_responses(db, grepl('S1', item_id))
 
  expect_true(dexter:::df_identical(r1, r2))

})
   