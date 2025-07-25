
library(dplyr)
library(DBI)

RcppArmadillo::armadillo_throttle_cores(1)


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
    rd$x = rd$x |>
      mutate_if(is.factor,as.character) |>
      arrange(person_id, booklet_id, item_id) |>
      as.data.frame()

    rd
  }
  a = prep(a)
  b = prep(b)

  expect_equal(a$summarised, b$summarised, info=info)
  
  expect_equal(a$design |> arrange(booklet_id,item_id), 
               a$design |> arrange(booklet_id,item_id), 
               info=info)
  
  expect_true(setequal(colnames(a$x), colnames(b$x)), info=info)
  
  expect_equal(a$x, b$x[,colnames(a$x)], info=info)
  
  invisible(a_)
}

expect_valid_respData = function(respData, msg='respData')
{
  expect_true(is.integer(respData$x$person_id) || is.factor(respData$x$person_id),
              info = sprintf("%s - x$person_id is not a factor but '%s'",
                             msg, typeof(respData$x$person_id)))
  
  expect_true(is.factor(respData$x$booklet_id),
              info = sprintf("%s - x$booklet_id is not a factor but '%s'",
                             msg, typeof(respData$x$booklet_id)))
  
  expect_true(is.factor(respData$design$booklet_id),
              info = sprintf("%s - design$booklet_id is not a factor but '%s'",
                             msg, typeof(respData$design$booklet_id)))
  
  expect_true(is.factor(respData$design$item_id),
              info = sprintf("%s - design$item_id is not a factor but '%s'",
                             msg, typeof(respData$design$item_id)))
  
  expect_true(is.integer(respData$x$booklet_score),
              info = sprintf("%s - x$isumSscore is not an integer but '%s'",
                             msg, typeof(respData$x$item_id)))
  
  # check factor levels
  expect_true(n_distinct(respData$design$booklet_id) == nlevels(respData$design$booklet_id),
              info="empty levels in booklet_id")
  expect_true(n_distinct(respData$design$item_id) == nlevels(respData$design$item_id),
              info='empty levels in item_id')
  
  
  
  if(!respData$summarised)
  {
    expect_true(is.factor(respData$x$item_id),
                info = sprintf("%s - x$item_id is not a factor but '%s'",
                               msg, typeof(respData$x$item_id)))
    
    expect_true(is.integer(respData$x$item_score),
                info = sprintf("%s - x$item_score is not an integer but '%s'",
                               msg, typeof(respData$x$item_id)))
    
    expect_false(is.unsorted(as.integer(respData$person_id)), info=sprintf("%s - person_id is unsorted", msg))
    
    split(as.integer(respData$x$booklet_id), respData$x$person_id) |>
      lapply(is.unsorted) |>
      unlist() |>
      any() |>
      expect_false(info=sprintf("%s - (person_id, booklet_id) is unsorted", msg))
    
    respData$x |>
      group_by(person_id, booklet_id) |>
      mutate(booklet_score2 = sum(item_score)) |>
      ungroup() |>
      summarise(res = all(booklet_score == booklet_score2)) |>
      pull(res) |>
      expect_true(info=sprintf("%s - booklet_score incorrect", msg))
    
  } 
  
  expect_false(is_grouped_df(respData$x), info = sprintf("%s - x is grouped", msg))
  expect_false(is_grouped_df(respData$design), info = sprintf("%s - design is grouped", msg))
  
  
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
  rules = distinct(rsp, item_id, item_score) |>
    mutate(response=item_score)
  
  db = start_new_project(rules, ':memory:')
  add_response_data(db, data = rename(rsp, response=item_score), design=distinct(rsp,booklet_id,item_id))
  
  
  get_resp_data(db) |> 
    expect_valid_respData() |>
    expect_equal_respData(
      expect_valid_respData(get_resp_data(rsp)))

  
  expect_error({f=fit_enorm(db)},'not connected')
  expect_error({f=fit_enorm(rsp)},'not connected')
  
  # non booklet safe merge, should still not be connected
  expect_error({f=fit_enorm(db, item_id!='3')},'not connected')
  expect_error({f=fit_enorm(rsp, item_id!='3')},'not connected')
  
  # merge over booklets (not fit because data is random)
  
  a = get_resp_data(db,merge_within_person=TRUE) |> 
    expect_valid_respData() |>
    expect_equal_respData(
      expect_valid_respData(get_resp_data(rsp, merge_within_person=TRUE)))
  
  
  expect_length(levels(a$design$booklet_id),1)
  
  expect_equal(
    rsp |>
      group_by(person_id) |>
      summarise(booklet_score=sum(item_score)) |>
      ungroup() |>
      inner_join(get_resp_data(rsp,merge_within_person=TRUE,summarised=TRUE)$x, by=c('person_id','booklet_score')) |>
      NROW(),
    50)
    
  
  close_project(db)
  
  
  # a set that should not be mergable
  rsp = tibble(person_id = rep(rep(1:50,each=20),2), 
               booklet_id = rep(1:2, each=1000),
               item_id = c(rep(1:20, 50),rep(11:30, 50)),
               item_score=sample(0:3,2000,replace=TRUE))
  
  # also make a database
  rules = distinct(rsp, item_id, item_score) |>
    mutate(response=item_score)
  
  db = start_new_project(rules, ':memory:')
  add_response_data(db, data=rename(rsp, response=item_score), design = distinct(rsp,booklet_id, item_id))
  
  
  expect_no_error(get_resp_data(rsp, merge_within_person=FALSE))
  expect_no_error(get_resp_data(db, merge_within_person=FALSE))
  
  expect_error(get_resp_data(rsp, merge_within_person=TRUE),'more than once')
  expect_error(get_resp_data(db, merge_within_person=TRUE),'more than once')
  
  expect_error(get_resp_data(rsp, merge_within_person=TRUE,summarised=TRUE),'more than once')
  expect_error(get_resp_data(db, merge_within_person=TRUE,summarised=TRUE),'more than once')
  
  
  
  close_project(db)
})

test_that('integers and factors',{
  x = tibble(item_id=rep(1:10,10),person_id=rep(c(1,4:12),each=10),item_score=sample(0:3,100,replace=TRUE))
  
  r = get_resp_data(x)
  
  tst = inner_join(mutate_if(x,is.numeric,as.integer),
                    mutate(r$x,across(c(person_id,item_id),\(i) as.integer(as.character(i)))),
                   by=c('person_id','item_id'), suffix=c('.in','.out'))
  
  expect_true(nrow(tst) == nrow(x) && all(tst$item_score.in == tst$item_score.out),label='integer match respdata')
  
  x = mutate(x,
             person_id = factor(as.character(person_id),levels=as.character(sample(1:100,100))),
             item_id = factor(as.character(item_id),levels=as.character(sample(1:100,100))))
  
  r = get_resp_data(x)
  
  tst = inner_join(x,r$x,by=c('person_id','item_id'), suffix=c('.in','.out'))
  
  expect_true(nrow(tst) == nrow(x) && all(tst$item_score.in == tst$item_score.out),label='factor match respdata')
  
})


test_that('input data.frames survives',  {

  # do new project, guarantees nice ordering
  db = start_new_project(verbAggrRules, ":memory:")
  add_booklet(db, verbAggrData, "agg")        
  
  r = get_responses(db)
  r2 = rlang::duplicate(r)
  
  v=get_resp_data(r,summarised=TRUE)
  v=get_resp_data(r,summarised=FALSE)
  
  expect_identical(r,r2, label="get_resp_data should not mutilate input")  
  
  v=get_resp_data(r, summarised=TRUE, protect_x=FALSE)
  
  expect(!all(r$item_score==r2$item_score), 'when protect_x is false we would like some input mutilation')
  
  
  close_project(db)        
})



test_that('get responses works correctly with predicates',
{

  db = open_project(test_path('testdata/verbAggression.db'))
  
  #two ways to do the same
  
  r1 = get_responses(db, item_id %like% 'S1%')

  r2 = get_responses(db, grepl('S1', item_id))
 
  expect_true(dexter:::df_identical(r1, r2))
  close_project(db)

})

test_that('empty levels are resolved',{
  dat = matrix(sample(0:2,100,TRUE),10,10)
  dat[1,] = NA
  dat[,3] = NA
  dat[5,5] = NA
  r1 = get_resp_data(dat)
  expect_valid_respData(r1)
  
  x=r1$x
  levels(x$item_id) = c(levels(x$item_id),'aap','bla')
  r2=get_resp_data(x)
  expect_valid_respData(r2)
  expect_equal_respData(r1,r2)
})
   
test_that('reconstruct output vartypes',{
  
  test = tibble(
    a = factor(letters[1:10], levels=rev(letters)),
    b = rnorm(10),
    c = 2L * 15:6,
    d = rev(letters[6:15]),
    e = as.double(21:30))
  
  info = dexter:::get_datatype_info(test, columns=colnames(test)) 
  
  out_char = dexter:::df_format(mutate(test, across(everything(),factor))) |>
    as_tibble()
  

  expect_identical(mutate(test, across(everything(), as.character)),out_char,label='out format defaults to character')
  
  
  out = dexter:::df_format(mutate(test, across(everything(),factor)), datatype_info=info) |>
    as_tibble()
  
  # doubles are never correctly preserved, unless they are integer like
  # # to do: should we even allow doubles as id's?
  expect_equal(test,out, label='vartypes are correctly preserved')
  
  
  
})


test_that('sql translation',
{

  trans = function(x, vars=NULL, variant='sqlite')
  {
    env = rlang::caller_env()
    p = eval(substitute(quote(x)))
    dexter:::translate_sql(dexter:::partial_eval(p, env=env, vars=vars),variant=variant) 
  }
  
  a=3
  vrs = c('a','b')
  expect_equal(trans(!!a==b, vrs), '3 = "b"')
  expect_equal(trans(local(a)==b, vrs), '3 = "b"')
  expect_equal(trans(a==b,'b'), '3 = "b"')
  expect_equal(trans(a==b, vrs), '"a" = "b"')
  
  expect_equal(trans(a == paste(b,'c'), vrs,'sqlite'), "\"a\" = \"b\"||' '||'c'")
  expect_equal(trans(a == paste(b,'c'), vrs, 'ansi'), "\"a\" = CONCAT_WS(' ',\"b\",'c')")
  
  # get
  v = 'gender'
  expect_equal(trans(get(v)=='bla','gender'),
               '"gender" = \'bla\'')
  
  
  # named and unnamed arguments
  expect_equal(trans(b==substr(a,4,7),vrs),
               trans(b==substr(a,stop=7,4),vrs))
  
  #missing arguments
  expect_error(trans(between(a,b),vrs))
  
  # named vector
  b = c(blaat=1,geit=2)
  
  expect_equal(trans(2 %in% b,'a','ansi'),'TRUE')
  expect_equal(trans(a %in% b,'a','ansi'),  '"a" in (1,2)')
  
  # ranges
  expect_equal(tolower(trans(a %in% 5:10,vrs)), "\"a\" in (5,6,7,8,9,10)")
  
  b=7
  expect_equal(tolower(trans(a %in% 5:.env$b,vrs)),"\"a\" in (5,6,7)")
  
  #casting
  expect_equal(trans(as.character(a),'a'), "CAST( \"a\" AS character )")
  expect_equal(trans(as.character(a)), "'3'")
  
  # indexing
  a = list(x=5,y=6)
  
  expect_equal(trans(x == a$x, 'x'),'"x" = 5')
  expect_equal(trans(x == a[['x']], 'x'),'"x" = 5')
  
  # combined c
  y=123
  expect_equal(trans(x %in% c(y,4,c(5,6)),'x'), trans(x %in% c(y,4,5,6),'x'))
  
  # substr
  expect_equal(trans(x == substr(d,5,6),c('x','d')), '"x" = substr( "d" , 5 , 2 )')
  expect_equal(trans(x == substr(d,5,y),c('x','d','y')), '"x" = substr( "d" , 5 , (1+("y")-(5)) )')
  expect_equal(trans(endsWith(mode,'o'),'mode'),"substr(\"mode\",1+length(\"mode\")-length('o'))='o'") 
  
  #unsure if we want to automatically unpack lists of length 1
  #expect_equal(trans(x == a['x'], 'x'),'"x" = 5')
  
  # .data and .env, quasi quotation
  a = 5L
  expect_equal(trans(a==.env$a),'5 = 5')
  expect_equal(trans(a==.env$a, vars='a'),'"a" = 5')
  expect_equal(trans(.data$a==.env$a, vars='a'),'"a" = 5')
  
  expect_equal(trans(.data$a==!!a,vars='a'),'"a" = 5')
  
  expect_error(trans(.data$a==.env$a),'not found')
})    
  

test_that('variable names cross sql',
{

  # variable names are lowercase in sql and do not support special characters such as a dot
  # We make no effort to support dots and such but we do make an effort to support case mismatch
  
  # If a variable does not exist in the db and does not exist in the environment 
  # but it does exists in the db with another case, it should work.
  
  db = start_new_project(verbAggrRules, ":memory:", person_properties=list(Gender='<NA>'))
  add_booklet(db, verbAggrData, "agg")
  
  
  expect_message({rsp = get_responses(db,Gender=='Male')},
                 'Gender.*gender')
  
  rsp1 = get_responses(db,gender=='Male')
  
  # force non sql evaluation by using grepl, use capital G->Gender
  expect_message({rsp2 = get_responses(db, grepl('^male',Gender,ignore.case = TRUE))},
                 'Gender.*gender')
  
  expect_identical(table(rsp$item_score), table(rsp1$item_score),
                   label='sql capital versus non capital var names, expect equal results. ')
  
  expect_identical(table(rsp1$item_score), table(rsp2$item_score),
                   label='case mismatch sql non sql should not cause a difference, expect equal results. ')
  
  
  # test if unknown names fail
  a = 1
  
  expect_error({get_responses(db,item_id==a | gndr=='Male')},
               "'gndr' not found")
  
  close_project(db)

})

RcppArmadillo::armadillo_reset_cores()
