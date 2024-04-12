context('check c++ level elsym')

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

  
  # # against dexter
  # dx = elsym:::add_zero(rev(a),rev(b),fl,fl)
  # 
  # dx_e_i = sapply(nit:1, \(i) dexter:::elsym(dx$b,dx$a,dx$bkfirst+1L,dx$bklast+1L,i))
  # 
  # # this works
  # max_rel_diff(dx_e_i,e_i[1:ns,])
  # max_rel_diff(dx_e_i,e_vnl)
  
  
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
  
  # # against dexter
  # 
  # dx = elsym:::add_zero(ab$item_score,ab$b,items$first,items$last)
  # 
  # dx_e_i = sapply(1:nit, \(i) dexter:::elsym(dx$b,dx$a,dx$bkfirst+1L,dx$bklast+1L,i))
  # 
  # #also works
  # max_rel_diff(dx_e_i,e_vnl)
  

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
  # save(ss,b,file=test_path('sufstats_dich.RData'))
  load(test_path('sufstats_dich.RData'))
  
  lbinom = dexter:::lbnm(max(ss$booklet$max_score)+3)
  
  bkfirst = as.integer(ss$design$first - 1L)
  bklast = as.integer(ss$design$last - 1L)
  
  # test expect
  
  ev = dexter:::Expect(b, ss$ssIS$item_score, bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit)
  ebnm = dexter:::Expect_binom(lbinom,b, ss$ssIS$item_score, bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit)
  
  expect_lt(max_rel_diff(ev,ebnm), 5e-15, label='binomial expect equal to regular in dichotomous case')
  
  
  # dx = elsym:::add_zero(ss$ssIS$item_score,b,bkfirst, bklast)
  # 
  # dx_ev = dexter:::E_bkl(dx$b, dx$a, dx$bkfirst, dx$bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit, use_mean=FALSE)
  # dx_ebnm = dexter:::E_bkl(dx$b, dx$a, dx$bkfirst, dx$bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit, use_mean=TRUE)
  # 
  # # goed
  # max_rel_diff(ev,dx_ev[-(sort(unique(dx$bkfirst)+1))])
  # max_rel_diff(ebnm,dx_ebnm[-(sort(unique(dx$bkfirst)+1))])
  
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
  
  
  # dx_h = matrix(0,2*nit,2*nit)
  # dx_e=double(2*nit)
  # 
  # dexter:::NR_bkl(dx$b, dx$a, dx$bkfirst, dx$bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit, 
  #        max(2*ss$booklet$nit), dx_e, dx_h, use_mean=FALSE)
  # 
  # 
  # # 1.5e-13, long double etc. it's believable
  # max_rel_diff(hv,dx_h[-(sort(unique(dx$bkfirst)+1)),-(sort(unique(dx$bkfirst)+1)) ])
  # 
  # dx_h = matrix(0,2*nit,2*nit)
  # dx_e=double(2*nit)
  # 
  # dexter:::NR_bkl(dx$b, dx$a, dx$bkfirst, dx$bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit, 
  #                 max(2*ss$booklet$nit), dx_e, dx_h, use_mean=TRUE)
  # 
  # # 5e-13 
  # max_rel_diff(hbnm,dx_h[-(sort(unique(dx$bkfirst)+1)),-(sort(unique(dx$bkfirst)+1)) ])
  
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
  # save(ss,b,file=test_path('sufstats_poly.RData'))
  
  load(test_path('sufstats_poly.RData'))
  
  
  lbinom = dexter:::lbnm(max(ss$booklet$max_score)+3)
  
  bkfirst = as.integer(ss$design$first - 1L)
  bklast = as.integer(ss$design$last - 1L)
  

  ev = dexter:::Expect(b, ss$ssIS$item_score, bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit)
  ebnm = dexter:::Expect_binom(lbinom,b, ss$ssIS$item_score, bkfirst, bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit)

  expect_lt(max_rel_diff(ev,ebnm), 5e-15, label='binomial expect equal to regular in polytomous case')
  

  # dx = elsym:::add_zero(ss$ssIS$item_score,b,bkfirst, bklast)
  # 
  # dx_ev = dexter:::E_bkl(dx$b, dx$a, dx$bkfirst, dx$bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit, use_mean=FALSE)
  # dx_ebnm = dexter:::E_bkl(dx$b, dx$a, dx$bkfirst, dx$bklast, ss$scoretab$N, ss$booklet$n_scores, ss$booklet$nit, use_mean=TRUE)
  # 
  # # ok
  # max_rel_diff(dx_ev,dx_ebnm)
  # max_rel_diff(ev,ebnm)
  

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