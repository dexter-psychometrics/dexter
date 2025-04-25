
library(dplyr)


RcppArmadillo::armadillo_throttle_cores(1)

max_rel_diff = function(a,b)
{
  dv = abs(as.double(b))
  dv = pmax(dv,1)
  max(abs(a-b)/dv)
}

# there is a greater risk of errors with 1 core
MAX_CORES = 1L

test_that('dichotomous elsym_i matches elsym',{

  nit = 60
  
  b = exp(rnorm(nit))
  a = sample(1:8,nit,TRUE)
  fl = 0:(nit-1)
  
  ns = sum(a) + 1L
  
  e_vnl = sapply(fl,function(i) dexter:::elsymC(b,a,fl,fl,i))
  
  e_i = dexter:::list_elsymiC(b,a,fl,fl)
  
  expect_lt(max_rel_diff(e_i$gi[1:ns,],e_vnl), 1e-14, label='dichotomous elsym_i equal to vanilla')
  
  g = dexter:::elsymC(rev(b),rev(a),fl,fl,-1)
  
  expect_lt(max_rel_diff(e_i$gfull[1:length(g)],g), 1e-14, label='dichotomous elsym_i gfull equal to vanilla')


  
  lbinom = dexter:::lbnm(ns+3)
  
  e_vnl = sapply(fl,function(i) dexter:::elsym_binomC(lbinom,b,a,fl,fl,i))
  
  e_i = dexter:::list_elsymi_binomC(lbinom,b,a,fl,fl)
  
  expect_lt(max_rel_diff(e_i$gi[1:ns,],e_vnl), 1e-14, label='dichotomous binomial elsym_i equal to vanilla')
  
  
  gbnm = dexter:::elsym_binomC(lbinom,b,a,fl,fl,-1)
  
  expect_lt(max_rel_diff(e_i$gfull[1:length(gbnm)],gbnm), 1e-14, label='dichotomous binomial elsym_i gfull equal to vanilla')
  
  
  expect_lt(max_rel_diff(g,gbnm*choose(ns-1,0:(ns-1))), 1e-13, label='binomial elsym approx equal to elsym')
  
  
  g = dexter:::elsym_binomC(lbinom,rep(1,length(b)),rep(1L,length(a)),fl,fl,-1)
  
  expect_true(all(abs(g-1)<1e-14),label='elsym binomial=1 for b=1')
  
  
})  

test_that('polytomous elsym_i matches elsym',{

  nit = 60
  
  ab = tibble(item_id=sprintf("i%03i",1:nit),ncat=sample(1:4,nit,TRUE)) |>
    rowwise() |>
    do({
        tibble(item_id=.$item_id,item_score=sample(1:(2*.$ncat),.$ncat),b=exp(rnorm(.$ncat)))
    }) |>
    ungroup() |>
    arrange(item_id,item_score)
  
  items = ab |>
    mutate(fl=row_number()-1) |>
    group_by(item_id) |>
    summarise(first=min(fl),last=max(fl),ms=max(item_score))
  
  ns = sum(items$ms)+1L
  
  e_vnl = sapply(1:nrow(items),function(i) dexter:::elsymC(ab$b,ab$item_score,items$first,items$last,i-1L))
  
  e_i = dexter:::list_elsymiC(ab$b,ab$item_score,items$first,items$last)$gi

  expect_lt(max_rel_diff(e_i[1:ns,],e_vnl), 1e-14, label='polytomous elsym_i equal to vanilla')
  


  lbinom = dexter:::lbnm(ns+3)
  
  e_vnl = sapply(1:nrow(items),function(i) dexter:::elsym_binomC(lbinom,ab$b,ab$item_score,items$first,items$last,i-1L))
  
  e_i = dexter:::list_elsymi_binomC(lbinom,ab$b,ab$item_score,items$first,items$last)$gi
  
  expect_lt(max_rel_diff(e_i[1:ns,],e_vnl), 1e-14, label='polytomous binomial elsym_i equal to vanilla')

  

})




test_that('dichotomous regular and binom hessian are equal',{
  
  # test data, no regenerate because there are some
  # data conditions that may (correctly) cause an error
  
  # dichotomous
  # nit  = 60
  # np = 3000
  # items = tibble(item_id=sprintf("i%03i",1:nit),item_score=sample(1:9,nit,TRUE),beta=rnorm(nit))
  # b = exp(-items$beta*items$item_score)
  # 
  # dat = dexter::r_score(items)(rnorm(np))
  # 
  # #unequal length booklets with first longer
  # dat[1:as.integer(np/2), as.integer(2*nit/3):nit] = NA_integer_
  # dat[(1L+as.integer(np/2)):np, 1:as.integer(nit/2)] = NA_integer_
  # ss = elsym:::sufStats_nrm(dexter::get_resp_data(dat))
  # 
  # save(ss,b,file=test_path('testdata/sufstats_dich.RData'))
  load(test_path('testdata/sufstats_dich.RData'))
  
  lbinom = dexter:::lbnm(max(ss$booklet$max_score)+3)
  
  bkfirst = as.integer(ss$design$first - 1L)
  bklast = as.integer(ss$design$last - 1L)
  
  # test expect
  
  ev = dexter:::Expect(b, ss$ssIS$item_score, bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit)
  ebnm = dexter:::Expect_binom(lbinom,b, ss$ssIS$item_score, bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit)
  
  expect_lt(max_rel_diff(ev,ebnm), 5e-15, label='binomial expect equal to regular in dichotomous case')
  
  
  nit=nrow(ss$ssI)
  
  ebnm = double(nit)
  ev = double(nit)
  hbnm = matrix(0,nit,nit)
  hv = matrix(0,nit,nit )
  
  
  dexter:::Hess_binom(lbinom,b, ss$ssIS$item_score, 
                  bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit, MAX_CORES, ebnm, hbnm)
  
  dexter:::Hess(b, ss$ssIS$item_score, 
            bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit, MAX_CORES, ev,hv)
  
  expect_true(max_rel_diff(ev,ebnm)< 5e-15, label='binomial h expect equal to regular in dichotomous case')
  
  # we have to accept some compounded machine error for second derivatives
  expect_true(max_rel_diff(hv,hbnm)< 1e-12, label='binomial hessian approx equal to regular in dichotomous case')
  
  

})  


test_that('polytomous regular and binom hessian are equal',{

  # test data, no regenerate because there are some
  # data conditions that may (correctly) cause an error
  
  # nit  = 40
  # np = 3000
  # 
  # items = tibble(item_id=sprintf("i%03i",1:nit),ncat=sample(1:4,nit,TRUE)) |>
  #   rowwise() |>
  #   do({
  #     tibble(item_id=.$item_id,item_score=sample(1:(2*.$ncat),.$ncat),beta=rnorm(.$ncat))
  #   }) |>
  #   ungroup() |>
  #   arrange(item_id,item_score)
  # 
  # b = exp(-items$beta*items$item_score)
  # 
  # 
  # dat = dexter::r_score(items)(rnorm(np))
  # 
  # dat[1:as.integer(np/2), 1:as.integer(nit/2)] = NA_integer_
  # dat[(1L+as.integer(np/2)):np, as.integer(2*nit/3):nit] = NA_integer_
  # 
  # ss = elsym:::sufStats_nrm(dexter::get_resp_data(dat))
  # 
  # save(ss,b,file=test_path('testdata/sufstats_poly.RData'))
  
  load(test_path('testdata/sufstats_poly.RData'))
  
  
  lbinom = dexter:::lbnm(max(ss$booklet$max_score)+3)
  
  bkfirst = as.integer(ss$design$first - 1L)
  bklast = as.integer(ss$design$last - 1L)
  

  ev = dexter:::Expect(b, ss$ssIS$item_score, bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit)
  ebnm = dexter:::Expect_binom(lbinom,b, ss$ssIS$item_score, bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit)

  expect_lt(max_rel_diff(ev,ebnm), 5e-15, label='binomial expect equal to regular in polytomous case')
  
  npar = n_distinct(ss$ssIS$item_id,ss$ssIS$item_score)
  
  ebnm = double(npar)
  ev = double(npar)
  hbnm = matrix(0,npar,npar)
  hv = matrix(0,npar,npar )
  
  
  dexter:::Hess_binom(lbinom,b, ss$ssIS$item_score, 
                     bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit, MAX_CORES, ebnm, hbnm)
  
  dexter:::Hess(b, ss$ssIS$item_score, 
               bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit,MAX_CORES, ev,hv)
  
  expect_true(max_rel_diff(ev,ebnm)< 5e-15, label='binomial h expect equal to regular in polytomous case')
  expect_true(max_rel_diff(hv,hbnm)< 1e-12, label='binomial hessian approx equal to regular in polytomous case')
  
})

test_that('vanilla sstable with normal elsym approx equal to long double version in c',{
  
  a=sample(1:2,16,replace=TRUE)
  b=rnorm(16)
  
  tst = dexter:::sstable_nrmC(a,b,0:7,0:7,8:15,8:15)
  
  g1 = dexter:::elsymC(b, a, 0:7, 0:7)
  g2 = dexter:::elsymC(b, a, 8:15, 8:15)
  g_all = dexter:::elsymC(b, a, 0:15, 0:15)
  
  
  m1 = sum(a[1:8])
  m2 = sum(a[9:16])
  
  m = matrix(0,m1+1,m2+1)
  
  for(s1 in 0:m1) 
    for(s2 in 0:m2)
      m[s1+1,s2+1] = g1[s1+1]*g2[s2+1]/g_all[s1+s2+1]
  
  
  expect_lt(max(abs(tst-m)/pmax(1,abs(tst))),1e-12, label='sstable two ways')
  
})
