
# to do: dexterities links

library(stringr)
library(RCurl)
library(dplyr)

# precompile blog function to prevent missing packages on git
# > pre_compile('ExponentialFamily_Blog.Rmd.orig')
# static images should be included in blog/img and be referred to with 'img/...'
# to do: this does not work anymore since moving to dexter-psychometrics

# I read all datasets from /Rdatasets/... 
# since not all of them can be made public it can be any dir on your computer outside of 
# the dexter directory but it's easier if we all keep to /Rdatasets
# dont use webshot or knitr::include_app( but an iframe
# internal function call in calculating score distribution, no contract with the user
pre_compile = function(blog_file_name)
{
  options(dplyr.summarise.inform=FALSE)
  if(!endsWith(blog_file_name,".Rmd.orig"))
    stop("file name should end in .Rmd.orig")
  
  if(dir.exists('figure'))
    stop('found figure dir in home folder')

  if(!check_sanity(file.path('vignettes','blog',blog_file_name)))
    stop("sanity check found problems")
  
  newf = file.path('vignettes','blog',gsub('\\.orig$','',blog_file_name))
  # select = dplyr::select
  # knitr::knit(file.path('vignettes','blog',blog_file_name), output = newf,
  #             envir=new.env())
  knit_separately(file.path('vignettes','blog',blog_file_name), output = newf)
  
  if(dir.exists('figure'))
  {
    fign = list.files('figure','png$')
    fig = lapply(fign, function(fn){
      fn = file.path('figure',fn)
      base64Encode(readBin(fn, "raw", file.info(fn)[1, "size"]), "txt")
    })
    names(fig) = fign
     
    suppressWarnings({txt = paste(readLines(newf),collapse='\n')})
    
    # for some reason it is arbitrary if images are exported in markdown or html format
    # prbl to do with fig margins being set in knitr
    img = str_extract_all(txt,'<img src="[^>]+>')
    for(i in unlist(img))
    {
      fn = str_extract(i,'(?<=src=")[^"]+')
      if(!is.na(fn) && !startsWith(fn,'img'))
      {
        j = gsub('src="[^"]+"',sprintf('src="data:image/png;base64,%s"', fig[[basename(fn)]]),i)
        txt = gsub(i,j,txt,fixed=TRUE)
      }
    }
    img = str_extract_all(txt,'!\\[[^]]+\\]\\([^)]+\\)')
    for(i in unlist(img))
    {
      fn = str_extract(i,'[^()]+\\.png')
      if(!is.na(fn) && !startsWith(fn,'img'))
      {
        txt = gsub(i,sprintf('<img src="data:image/png;base64,%s" />', fig[[basename(fn)]]),txt,fixed=TRUE)
      }
    }

    writeLines(txt,newf)
    unlink('./figure',force=TRUE,recursive=TRUE)
  }
}

knit_separately = function(...)
{
  callr::r(function(...) knitr::knit(..., envir = globalenv()), args = list(...), show = TRUE)
}

check_sanity = function(fn)
{
  sane = TRUE
  suppressWarnings({txt = paste(readLines(fn),collapse='\n')})
  
  ddd = str_extract_all(txt,'\\w+:::\\w+',simplify=TRUE) %>%
    as.character() %>%
    unique()
  if(length(ddd)>0)
  {
    message("Found the following internal function calls:")
    print(ddd)
    sane = FALSE
  }
  # will add more whenever dexter has a (soft) deprecation or similar
  sane
}


# f = list.files('vignettes/blog',pattern='Rmd\\.orig$',full.names=TRUE)
# find_install_libs(f)
find_install_libs = function(file_names)
{
  libs = sapply(file_names, function(fn) {
    txt = suppressWarnings({paste(readLines(fn),collapse='\n')})
    str_match_all(txt,'library\\((\\w+)\\)')[[1]][,2]
  }) %>%
    unlist() %>%
    unique()
  
  install.packages(setdiff(libs, rownames(installed.packages())))  
}


compile_all = function()
{
  f = list.files('vignettes/blog',pattern='Rmd\\.orig$',full.names=FALSE)
  invisible(sapply(f,pre_compile))
}

compile_all_new = function()
{
  f = list.files('vignettes/blog',pattern='Rmd\\.orig$',full.names=FALSE)
  x = list.files('vignettes/blog',pattern='Rmd$',full.names=FALSE) %>% 
    gsub('Rmd$','Rmd.orig',.)
  
  f=setdiff(f,x)
  
  invisible(sapply(f,pre_compile))
}

blog_yaml = function()
{
  s = function(n=6) paste0(rep(' ',n),collapse='')

  
  fls = tibble(fn=list.files('vignettes/blog',pattern='.Rmd$',full.names=FALSE)) %>%
    mutate(fn=gsub('\\.Rmd$','',fn)) %>%
    mutate(bn = gsub('^\\d{4}-\\d{2}-\\d{2}-','',fn)) %>%
    mutate(bn=gsub('(?<=\\d)-(?=\\d)','.',bn,perl=TRUE)) %>%
    mutate(bn=gsub('-',' ',bn,fixed=TRUE)) %>%
    mutate(yml = paste0(s(6),'- text: ',bn,'\n',s(8),'href: articles/blog/',fn)) %>%
    arrange(desc(fn))
  
  cat(paste(fls$yml,collapse='\n'))
}



