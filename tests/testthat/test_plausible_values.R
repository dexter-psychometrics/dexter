context('test plausible_values')

test_that('populations work',{
  db = open_project('../verbAggression.db')
  
  #artificially create two overlapping booklets
  f2 = fit_enorm(db, (grepl('dxP\\d\\d$',person_id,perl=TRUE) & item_position<=16) | 
                     (!grepl('dxP\\d\\d$',person_id,perl=TRUE) & item_position>8))
  
  expect_true(length(f2$inputs$bkList)==2)
  
  pv=plausible_values(db, f2, 
                        predicate={(grepl('dxP\\d\\d$',person_id,perl=TRUE) & item_position<=16) | 
                                  (!grepl('dxP\\d\\d$',person_id,perl=TRUE) & item_position>8)},
                        covariates='gender')
  
  expect_true(mean(pv[pv$gender=='Male',]$PV1) > mean(pv[pv$gender=='Female',]$PV1))
  
  
  dbDisconnect(db)
})
  