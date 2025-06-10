

test_that('enorm plot',{
  db = open_project(test_path('testdata/verbAggression.db'))
  

  ff = fit_enorm(db)
  
  tht = ability(db,ff,method='MLE') |>
    filter(is.finite(theta)) |>
    arrange(theta)
  
  res = plot(ff,items='S3DoScold',plot=FALSE,nbins=6)
  
  tht$abgroup=rep(1:6,res$n)
  
  # binning should not put same scores in different bins
  expect_equal(n_distinct(tht$booklet_score,tht$abgroup),n_distinct(tht$booklet_score))
  
  tst = tht |>
    group_by(abgroup) |>
    summarise(gr_theta = mean(theta),n=n()) |>
    inner_join(res,by='abgroup', suffix=c('.abl','.plt'))
  
  expect_lt(max(abs(tst$gr_theta.abl-tst$gr_theta.plt)), 1e-14)
  
  
})

test_that('fstr',{
  expect_equal(dexter:::fstr('$bla, $b',list(bla='some string')),"some string, $b")
  expect_equal(dexter:::fstr('$bla:.1f, $b',list(bla=3.2423, b='item')),"3.2, item")
})
