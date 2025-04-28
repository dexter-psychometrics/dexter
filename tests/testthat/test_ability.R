

RcppArmadillo::armadillo_throttle_cores(1)


test_that('inconsistencies between data and parms are handled correctly',{

  db = open_project(test_path('testdata/verbAggression.db'))
  
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


test_that('verbAgg abilities', {

  db = open_project(test_path('testdata/verbAggression.db'))
  f = fit_enorm(db)
  
  
  # check ability mle is inverse of expected_score
  es = expected_score(f)
  expect_lt(
    ability_tables(f) |>
      filter(is.finite(theta)) |>
      mutate(error = abs(booklet_score - es(theta))) |>
      pull(error) |>
      mean(),
    0.00001,
    label = "ability_tables mle on average estimated to within .00001 of test_score")
  
  
  nscores = get_rules(db) |>
    group_by(item_id) |>
    summarize(m=max(item_score)) |>
    ungroup() |>
    pull(m) |>
    sum() + 1
  
  test_cases = list(MLE = c('MLE','normal'), WLE = c('WLE','normal'), EAP_normal = c('EAP','normal'), EAP_J = c('EAP','Jeffreys'))
  
  res = lapply(test_cases, function(s){ ability_tables(f, method = s[1], prior = s[2])})
  
  expect_false(any(sapply(lapply(res,'[[','theta'), is.unsorted)), info='abilities not increasing verbAgg')
  
  theta = do.call(cbind,lapply(res,'[[','theta'))
  expect_true(sum(!apply(theta,1,is.finite)) == 2 && !any(is.finite(theta[c(1,nscores),1])), info='inifinity only in MLE')
  theta[!is.finite(theta)] = NA
  r = cor(theta,use='pairwise')
  expect_true(all(r >= .99), info='high correlation ability estimates one booklet')
  expect_true(all(r[upper.tri(r)] < 1), info='different ability methods are different')
  
  
  dbDisconnect(db)
  
})

test_that('ability WLE compared to Norman theta', {
  
  load(test_path('testdata/theta.RData'))
  # test against theta-unweighted

    a = ability_tables(item_param, design=design[sample(nrow(design)),], method='WLE')
  
  tst = inner_join(a,os_theta, by=c('booklet_id','booklet_score'))
  
  expect_true(nrow(tst) == nrow(a), info='proper booklet scores in ability')
  
  expect_lt(max(abs(tst$theta.x-tst$theta.y)), 1e-10, label='wle equal to theta norman, difference')
  # small differences multiply in the se, so < 1e6 is enough
  expect_lt(max(abs(tst$se.x-tst$se.y)), 1e-6, label='se wle equal to theta norman, difference')
  
})

test_that('designs are handled correctly', {
  set.seed(123)
  items = tibble(item_id=1:20, ncat=sample(1:3,20,replace=TRUE)) |>
    group_by(item_id) |>
    do({
      tibble(item_score = cumsum(sample(1:3,.$ncat[1], replace=TRUE)),
             beta = sort(rnorm(.$ncat)))
    }) |>
    ungroup()
  
  nit = sample(10:20,3)
  design = tibble(booklet_id=factor(rep(letters[1:3],nit),levels=letters), item_id=unlist(sapply(nit,sample,x=1:20)))
  
  a1 = ability_tables(items,design=design, method='WLE')
  
  #make sure one category does not occur to see if it is correctly handled
  items2 = items |>
    add_count(item_id) |>
    arrange(desc(n), desc(item_score)) |>
    slice(-1)

  dat = r_score(items2)(rnorm(300)) 
  
  f = fit_enorm(dat,fixed_param=items)
  
  a2 = ability_tables(f,design=design, method='WLE')
  
  a3 = ability_tables(coef(f),design=design, method='WLE')
  
  expect_true(df_join_equal(a1,a2,a3,join_by=c('booklet_id','booklet_score')))
  

  
})



RcppArmadillo::armadillo_reset_cores()


