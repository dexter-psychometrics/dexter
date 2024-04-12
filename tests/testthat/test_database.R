context('test database')

library(dplyr)
library(DBI)
library(RSQLite)

RcppArmadillo::armadillo_throttle_cores(1)

verbAggCopy = function(pth = test_path('verbAggression.db'))
{
  con = dbConnect(SQLite(), ":memory:")
  db = open_project(pth)
  
  sqliteCopyDatabase(db, con)
  
  dbDisconnect(db)
  return(con)
}

test_that('rule updates and sanity checks',
{
  db = verbAggCopy()
  
  rules = get_rules(db)
  expect_true(nrow(rules) == 72, 'expect 72 rules in verbal aggression database')
  
  ## test no-op
  expect_output(touch_rules(db, slice(rules,1:10)),
                 '(changed\\D*0\\D*added\\D*0\\D*)|(addeded\\D*0\\D*changed\\D*0\\D*)$',
                 info = 'expect message saying 0 changes when resupplying an existing rule')
  
  expect_true(dexter:::df_identical(arrange(rules,item_id,response),arrange(get_rules(db),item_id,response)),
              info='expect 0 actual changes when resupplying existing rule')
  
  
  ## test update score
  expect_output(touch_rules(db, filter(rules, item_id=='S1DoCurse' & response ==1) |> mutate(item_score=2) ),
                 '(changed\\D*1\\D*added\\D*0\\D*)|(added\\D*0\\D*changed\\D*1\\D*)$',
                 info = 'expect message saying 1 change and 0 added when modyfying item_score')
  
  expect_true(filter(get_rules(db), item_id=='S1DoCurse' & response ==1)$item_score == 2,
              info = 'expect item score changed to 2')
  
  ##### test sanity checks #####
  
  ## less_than_two_scores
  
  # update set does not have to pass this sanity check on it's own
  expect_no_error(touch_rules(db, filter(rules, item_id=='S1DoCurse' & response ==1) |> mutate(item_score=0)),
                  message='expect allowed to set an item score to 0')
  
  # but it should shout when every score is set to 0 for an item
  expect_output(
    expect_error(touch_rules(db, filter(rules, item_id=='S1DoCurse' & response ==2) |> mutate(item_score=0))),
    regexp='S1DoCurse\\s+TRUE',
    info='should complain about all scores set to 0')

  expect_false(all(filter(get_rules(db), item_id == 'S1DoCurse') |> pull(item_score) == 0),
               info='expect faulty update not to have occurred')
  
  ## duplicated_responses
  
  # resupplying an existing response should be a no-op, as tested above
  # but supplying an option two times in the new rules should be stopped
  expect_output(
    expect_error(touch_rules(db, tibble(item_id='S2DoCurse',response=c(5,5),item_score=c(3,5)))),
    regexp='S2DoCurse\\s+FALSE\\s+TRUE',
    info = 'should complain about duplicate responses')
  
  ## min_score_not_zero
  # change min score of existing rule away from 0
  expect_output(
    expect_error(touch_rules(db, tibble(item_id='S3DoCurse',response=0,item_score=5))),
    regexp='S3DoCurse\\s+FALSE\\s+FALSE\\s+TRUE',
    info='should complain about min score not 0')
  
  # we should be allowed to add a new rule with a score not 0
  # because this sanity check should be done in conjunction with the data already in the db
  expect_output(touch_rules(db, tibble(item_id='S3DoCurse',response=5,item_score=5)),
                 '(changed\\D*0\\D*added\\D*1\\D*)|(added\\D*1\\D*changed\\D*0\\D*)$',
                 info = 'expect message saying 0 change and 1 added when adding a new rule for an item')
  
  
  # but for a new item that should be a problem
  expect_output(
    expect_error(touch_rules(db, tibble(item_id='new',response=0,item_score=5))),
    regexp='new\\s+TRUE\\s+FALSE\\s+TRUE',
    info='should complain about min_score_not_zero and less_than_two_scores')
  
  expect_output(
    expect_error(touch_rules(db, tibble(item_id='new',response=c(0,1),item_score=c(2,5)))),
    regexp='new\\s+FALSE\\s+FALSE\\s+TRUE',
    info='should complain about min_score_not_zero and less_than_two_scores')
  
  ## less_than_two_scores
  # less than two scores is already covered above I think
  dbDisconnect(db)
})


test_that('adding person and item properties',
{
  # items
  db = verbAggCopy()
  
  expect_error({add_item_properties(db,tibble(a=1,b=2))},'item_id')
  
  # item id "1" does not exist
  expect_output({add_item_properties(db,tibble(item_id=1,b=2))}, '0 items')
  
  expect_output({add_item_properties(db,tibble(item_id='S4DoCurse', BlAme='society', news='olds'))},
                '2 item properties? for 1 items? added or updated', perl=TRUE)

  items = get_items(db)
  
  expect_true(all(items |> 
                    filter(item_id=='S4DoCurse') |> 
                    select(blame,news) |> 
                    unlist() == c("society", "olds" )))
  
  expect_output({add_item_properties(db, default_values = list(a=4L))},
                '1 new item_properties defined')
  
  # can only define once
  expect_output({add_item_properties(db, default_values = list(a=4L))},
                '0 new item_properties defined')
  
  items = get_items(db)
  
  expect_true(all(items$a == 4) && class(items$a) == 'integer')
  
  #persons
  
  expect_error({add_person_properties(db,tibble(a=1,b=2))},'person_id')
  
  # person "1" does not exist
  expect_output({add_person_properties(db,tibble(person_id=1,b=2))}, '0 persons')
  
  
  
  expect_output({add_person_properties(db,tibble(person_id=get_persons(db)$person_id[1], gender='unchecked', news='olds'))},
                '2 person properties? for 1 persons? added or updated', perl=TRUE)
  
  persons = get_persons(db)
  
  expect_true(all(persons |> 
                    filter(person_id=='dxP1') |> 
                    select(gender,news) |> 
                    unlist() == c("unchecked", "olds" )))
  
  expect_output({add_person_properties(db, default_values = list(x=4L))},
                '1 new person_properties defined')
  
  # can only define once
  expect_output({add_person_properties(db, default_values = list(x=4L))},
                '0 new person_properties defined')
  
  persons = get_persons(db)
  
  expect_true(all(persons$x == 4) && class(persons$x) == 'integer')
  dbDisconnect(db)
})


test_that('add_response_data',{
  db = verbAggCopy()
  
  dat = tibble(person_id='x',item_id="S1DoScold",response='2',booklet_id='agg')

  expect_message(expect_error({add_response_data(db, dat, missing_value='_missing_')},
                              'Unknown responses'),
                 '_missing_.+rules')
  
  i=get_items(db)
  i$response='_missing_'
  i$item_score=0
  touch_rules(db,i)
  
  expect_output({add_response_data(db,dat,missing_value='_missing_')},'23 missing responses')
  
  dbDisconnect(db)
})

test_that('add_booklet',{
  db = verbAggCopy()
  
  prs = get_persons(db)$person_id[1:5]
  dat = verbAggrData[1:5,]
  dat$person_id = prs
  
  expect_error({add_booklet(db, dat, "agg")},'overlap.+already imported',label='correct error for already existing booklet-persons')   
  
  dat$person_id='we are all the same person'
  
  expect_error({add_booklet(db, dat, "agg")},'person.+ duplicate',label='correct error for duplicate persons')   
  
  
  dbDisconnect(db)
})

RcppArmadillo::armadillo_reset_cores()