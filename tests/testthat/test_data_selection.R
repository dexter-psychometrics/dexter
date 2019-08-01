context('test data selection')

library(dplyr)
library(DBI)




expect_no_error = function(object, info=NULL) expect_error(object, regexp=NA, info=info)



# equivalent
expect_equal_respData = function(a,b, info='equal respData', ignore_booklet_levels = TRUE)
{
  a_ = a
  prep = function(rd)
  {
    if(ignore_booklet_levels)
    {
      rd$x$booklet_id = as.integer(rd$x$booklet_id)
      rd$design$booklet_id = as.integer(rd$design$booklet_id)
    }
    rd$x$person_id = as.character(rd$x$person_id)
    rd$design = as.data.frame(mutate_if(rd$design, is.factor,as.character))
    rd$x = rd$x %>%
      mutate_if(is.factor,as.character) %>%
      arrange(person_id, booklet_id, item_id) %>%
      as.data.frame()

    rd
  }
  a = prep(a)
  b = prep(b)

  expect_equal(a$summarised, b$summarised, info=info)
  
  expect_equal(a$design %>% arrange(booklet_id,item_id), 
               a$design %>% arrange(booklet_id,item_id), 
               info=info)
  
  expect_true(setequal(colnames(a$x), colnames(b$x)), info=info)
  
  expect_equal(a$x, b$x[,colnames(a$x)], info=info)
  
  invisible(a_)
}

# to do: check no grouping etc
expect_valid_respData = function(respData, msg='respData')
{
  expect_true(is.integer(respData$x$person_id) || is.factor(respData$x$person_id),
              info = sprintf("%s - x$person_id is not a factor but '%s'",
                             info, typeof(respData$x$person_id)))
  
  expect_true(is.factor(respData$x$booklet_id),
              info = sprintf("%s - x$booklet_id is not a factor but '%s'",
                             info, typeof(respData$x$booklet_id)))
  
  expect_true(is.factor(respData$design$booklet_id),
              info = sprintf("%s - design$booklet_id is not a factor but '%s'",
                             info, typeof(respData$design$booklet_id)))
  
  expect_true(is.factor(respData$design$item_id),
              info = sprintf("%s - design$item_id is not a factor but '%s'",
                             info, typeof(respData$design$item_id)))
  
  expect_true(is.integer(respData$x$booklet_score),
              info = sprintf("%s - x$isumSscore is not an integer but '%s'",
                             info, typeof(respData$x$item_id)))
  
  # to do: check factor levels
  
  if(!respData$summarised)
  {
    expect_true(is.factor(respData$x$item_id),
                info = sprintf("%s - x$item_id is not a factor but '%s'",
                               info, typeof(respData$x$item_id)))
    
    expect_true(is.integer(respData$x$item_score),
                info = sprintf("%s - x$item_score is not an integer but '%s'",
                               info, typeof(respData$x$item_id)))
    
    expect_false(is.unsorted(as.integer(respData$person_id)), info=sprintf("%s - person_id is unsorted", info))
    
    split(as.integer(respData$x$booklet_id), respData$x$person_id) %>%
      lapply(is.unsorted) %>%
      unlist() %>%
      any() %>%
      expect_false(info=sprintf("%s - (person_id, booklet_id) is unsorted", info))
    
    respData$x %>%
      group_by(person_id, booklet_id) %>%
      mutate(booklet_score2 = sum(item_score)) %>%
      ungroup() %>%
      summarise(res = all(booklet_score == booklet_score2)) %>%
      pull(res) %>%
      expect_true(info=sprintf("%s - booklet_score incorrect", info))
    
  }
  
  
  invisible(respData)
}

test_that('merging works',
{
  # a set connected over persons only
  rsp = tibble(person_id = rep(rep(1:50,each=20),2), 
               booklet_id = rep(1:2, each=1000),
               item_id = c(rep(1:20, 50),rep(21:40, 50)),
               item_score=sample(0:3,2000,replace=TRUE))
  
  # also make a database
  rules = distinct(rsp, item_id, item_score) %>%
    mutate(response=item_score)
  
  db = start_new_project(rules, ':memory:')
  add_response_data(db, rename(rsp, response=item_score))
  
  get_resp_data(db) %>% 
    expect_valid_respData() %>%
    expect_equal_respData(
      expect_valid_respData(get_resp_data(rsp)))

  
  expect_error({f=fit_enorm(db)},'not connected')
  expect_error({f=fit_enorm(rsp)},'not connected')
  
  # non booklet safe merge, should still not be connected
  expect_error({f=fit_enorm(db, item_id!='3')},'not connected')
  expect_error({f=fit_enorm(rsp, item_id!='3')},'not connected')
  
  # merge over booklets (not fit because data is random)
  
  a = get_resp_data(db,merge_within_person=TRUE) %>% 
    expect_valid_respData() %>%
    expect_equal_respData(
      expect_valid_respData(get_resp_data(rsp, merge_within_person=TRUE)))
  
  
  expect_length(levels(a$design$booklet_id),1)
  close_project(db)
  
  
  # a set that should not be mergable
  rsp = tibble(person_id = rep(rep(1:50,each=20),2), 
               booklet_id = rep(1:2, each=1000),
               item_id = c(rep(1:20, 50),rep(11:30, 50)),
               item_score=sample(0:3,2000,replace=TRUE))
  
  # also make a database
  rules = distinct(rsp, item_id, item_score) %>%
    mutate(response=item_score)
  
  db = start_new_project(rules, ':memory:')
  add_response_data(db, rename(rsp, response=item_score))
  
  
  expect_no_error(get_resp_data(rsp, merge_within_person=FALSE))
  expect_no_error(get_resp_data(db, merge_within_person=FALSE))
  
  expect_error(get_resp_data(rsp, merge_within_person=TRUE),'common items')
  expect_error(get_resp_data(db, merge_within_person=TRUE),'common items')
  
  expect_error(get_resp_data(rsp, merge_within_person=TRUE,summarised=TRUE),'common items')
  expect_error(get_resp_data(db, merge_within_person=TRUE,summarised=TRUE),'common items')
  
  close_project(db)
})




# to also do: check parms and profiles
            





test_that('get responses works correctly with predicates',
{
  
  db = open_project('../verbAggression.db')
  
  #two ways to do the same
  
  r1 = get_responses(db, item_id %like% 'S1%')

  r2 = get_responses(db, grepl('S1', item_id))
 
  expect_true(dexter:::df_identical(r1, r2))
  close_project(db)

})
   