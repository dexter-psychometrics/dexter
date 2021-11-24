
# to do: dexterities links

library(stringr)
library(RCurl)
library(dplyr)

# precompile blog function to prevent missing packages on git
# > pre_compile('ExponentialFamily_Blog.Rmd.orig')
# static images should be included in blog/img and be referred to with 'img/...'
# to do: this does not work anymore since moving to dexter-psychometrics
# refer to blog form within blog [](name)
# refer to blog from vignette: full link [](https://dexter-psychometrics.github.io/dexter/articles/blog/...)
# refer to vignette from blog [](../name)
# refer to function `function_name`
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
  
  find_install_libs(file.path('vignettes','blog',blog_file_name))
  
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
  
  # check for files that are not in Rdatasets
  # will add more whenever dexter has a (soft) deprecation or similar
  sane
}


# f = list.files('vignettes/blog',pattern='Rmd\\.orig$',full.names=TRUE)
# find_install_libs(f)
find_install_libs = function(file_names)
{
  libs = sapply(file_names, function(fn) {
    txt = suppressWarnings({paste(readLines(fn),collapse='\n')})
    str_match_all(txt,'library\\(([^\\)]+)\\)')[[1]][,2]
  }) %>%
    unlist() %>%
    unique()
  
  libs = setdiff(libs, rownames(installed.packages()))
  if(length(libs)>0)
    install.packages(libs)  
  invisible(NULL)
}


compile_all = function()
{
  f = list.files('vignettes/blog',pattern='Rmd\\.orig$',full.names=FALSE)
  invisible(sapply(f,pre_compile))
}

compile_all_new = function()
{
  orig = list.files('vignettes/blog',pattern='Rmd\\.orig$',full.names=FALSE)
  rmd = list.files('vignettes/blog',pattern='Rmd$',full.names=FALSE)
  names(rmd) = gsub('Rmd$','Rmd.orig',rmd)
  
  p = function(f) file.path('vignettes','blog',f)
  
  for(f in orig)
  {
    if(!f %in% names(rmd) || file.mtime(p(f)) > file.mtime(p(rmd[f])))
    {
      print(f)
      pre_compile(f)
    } 
  }
}

blog_yaml = function()
{
  s = function(n=6) paste0(rep(' ',n),collapse='')

  
  fls = tibble(fn=list.files('vignettes/blog',pattern='.Rmd$',full.names=FALSE)) %>%
    filter(grepl('^\\d',fn)) %>%
    mutate(fn=gsub('\\.Rmd$','',fn)) %>%
    mutate(bn = gsub('^\\d{4}-\\d{2}-\\d{2}-','',fn)) %>%
    mutate(bn=gsub('(?<=\\d)-(?=\\d)','.',bn,perl=TRUE)) %>%
    mutate(bn=gsub('-',' ',bn,fixed=TRUE)) %>%
    mutate(yml = paste0(s(6),'- text: ',bn,'\n',s(8),'href: articles/blog/',fn)) %>%
    arrange(desc(fn))
  
  cat(paste(fls$yml,collapse='\n'))
}


blurb_info = function(fns, n_lines=10)
{
  lapply(fns, function(fn)
  {
    lines = readLines(file.path('vignettes','blog',fn))
    lines = lines[nchar(trimws(lines))>0]
    start = max(which(grepl('^---',lines))) + 1
    if(startsWith(lines[start],'#')) start = start+1
    end = min(start + n_lines - 1,length(lines))
    if(end != length(lines))
    {
      if(sum(lines[start:end]=='$$') %% 2 == 1)
      {
        end = end + min(which(lines[(end+1):length(lines)] == '$$'))
      }
      if(sum(startsWith(lines[start:end],'```')) %% 2 == 1)
      {
        end = end + min(which(lines[(end+1):length(lines)] == '```'))
      }
    }

    tibble(
      filename = fn,
      blurb = paste(lines[start:end],collapse='\n'),
      abbreviated = end != length(lines),
      author = gsub(' * ',', ', trimws(str_extract(lines[startsWith(lines,'author:')][1],'(?<=:).+')) , fixed=TRUE),
      title = trimws(str_extract(lines[startsWith(lines,'title:')][1],'(?<=:).+')),
      href = gsub('.Rmd','',fn,fixed=TRUE)
    )
  }) %>%
    bind_rows()
}



make_blog_index = function()
{
  fn = list.files('vignettes/blog',pattern='.Rmd$',full.names=FALSE)
  fn = fn[grepl('^\\d',fn)]
  
  content = blurb_info(fn, 4) %>%
    mutate(rdmore = if_else(abbreviated,sprintf('<a href="%s">read more...</a>',href),'')) %>%
    mutate(txt = sprintf('<h2><a href="%s">%s</a></h2><p class="blog-authors">%s</p><p>\n%s\n</p>%s',
                         href,title,author,blurb,rdmore)) %>%
    arrange(desc(filename)) %>%
    pull(txt)
  
  style = 'img{max-width:8cm;max-height:8cm;display:block;}\n#refs{display:none;}\np.blog-authors{font-style:italic}'
  
  cat(sprintf('---\ntitle: Dexterities\nbibliography: dexter.bib\n---\n\n<style>\n%s\n</style>\n\n%s', style, paste(content,collapse='\n')),
      file='vignettes/blog/index.Rmd')
}



