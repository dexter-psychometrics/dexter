
library(dplyr)

RcppArmadillo::armadillo_throttle_cores(1)

tia_equal = function(tia1,tia2, tolerance = 1e-8)
{
  if(!setequal(colnames(tia1),colnames(tia2)))
    return(FALSE)
  
  cmp = inner_join(tia1,tia2,by=intersect(c('booklet_id','item_id'), colnames(tia1)))
  
  if(nrow(cmp) != nrow(tia1) || nrow(tia1) != nrow(tia2))
    return(FALSE)
  
  for(cn in setdiff(colnames(tia1), c('booklet_id','item_id')))
  {
    if(is.integer(tia1[[cn]]))
    {
      if(!all(cmp[[paste0(cn,'.x')]] == cmp[[paste0(cn,'.x')]]))
        return(FALSE)
    }
    if(!all(abs(cmp[[paste0(cn,'.x')]] - cmp[[paste0(cn,'.x')]]) <= tolerance))
      return(FALSE)
  }
  TRUE
}

test_that('tia computations are correct',{
  set.seed(123)
  
  #--- prepare test data
  nit_dich = 70L
  nit_poly = 10L
  nit = nit_poly + nit_dich
  N = 5000L
  
  
  items = tibble(item_id = sprintf('itm%03i',1:nit_dich), 
                 item_score = sample(1:3,nit_dich,replace=TRUE),
                 beta = rnorm(nit_dich))
  
  # a couple polytomous items
  
  items = rbind(items,
                tibble(item_id = rep(sprintf('itm_poly%03i',1:nit_poly), each=2),
                       item_score = rep(1:2,nit_poly),
                       beta = rnorm(2L*nit_poly)))
  
  
  dat = r_score(items)(rnorm(N))
  # remove disco's dich items
  dat[,!grepl('poly',colnames(dat))][dat[,!grepl('poly',colnames(dat))] > 1L] = 1L
  
  # incomplete design, 2 booklets with overlap
  dat[seq(1,N,2), seq(1,nit,3)] = NA_integer_
  dat[seq(2,N,2), seq(2,nit,3)] = NA_integer_
  
  # items without variation give warnings in the dplyr code below so we remove them if necessary
  dat = dat[,apply(dat,2,sd,na.rm=TRUE)>0]
  nit = ncol(dat)
  
  # to lf data.frame for dplyr
  rsp = get_responses(dat, columns = c("person_id", "item_id", "item_score")) |>
    mutate(booklet_id = as.character(1L+(as.integer(person_id) %% 2L == 0)))

  #--- test
  # manual tia separate per booklet
  tia_dpl = rsp |>
    group_by(person_id, booklet_id) |>
    mutate(booklet_score = sum(item_score)) |>
    ungroup() |>
    group_by(booklet_id,item_id) |>
    summarise(mean_score=mean(item_score),
              sd_score=sd(item_score),
              max_score=max(item_score),
              rit = cor(item_score,booklet_score),
              rir = cor(item_score,booklet_score-item_score),
              n_persons=n()) |>
    group_by(item_id) |>
    mutate(pvalue = mean_score/max(max_score))
  
  tia_dx = tia_tables(rsp)
  
  # statistics by booklet item are the basis for all other computations and are heavily optimized in cpp
  # if there is an error anywhere it will be here
  expect_true(tia_equal(tia_dpl,tia_dx$items), label='tia raw equal to dplyr')
  
  # combining sd's is a little tricky and was done wrong in the past
  
  tia_dx = tia_tables(rsp,type='averaged')
  
  sd_item = rsp |>
    group_by(item_id) |>
    summarise(sd_score=sd(item_score))
  
  tia_dpl = tia_dpl |>
    group_by(item_id, max_score) |>
    summarise(n_booklets = n(),
              mean_score = weighted.mean(mean_score, n_persons),
              rit = weighted.mean(rit, n_persons),
              rir = weighted.mean(rir, n_persons),
              pvalue = weighted.mean(pvalue, n_persons),
              nps = sum(n_persons)) |>
    ungroup() |>
    rename(n_persons='nps') |>
    inner_join(sd_item,by='item_id')
  
  
  expect_true(tia_equal(tia_dx$items,tia_dpl), label='tia averaged equal to dplyr')
  
  # this is just a pivot, see if it does not raise a trivial error
  tia_dx = tia_tables(rsp,type='compared')
  
  # see if lack of variation is handled correctly
  
  item1 = rsp$item_id[1]
  
  rsp = mutate(rsp,item_score = if_else(item_id==item1,1L,item_score) )
  
  expect_warning({tia_dx = tia_tables(rsp)},
                 regexp = '^.+without score variation.+$')
  
  expect_true(all(filter(tia_dx$items,item_id==item1)$sd_score == 0))
})

RcppArmadillo::armadillo_reset_cores()
