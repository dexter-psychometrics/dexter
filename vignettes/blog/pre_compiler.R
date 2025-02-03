

library(stringr)
library(RCurl)
library(dplyr)


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
    pcap = ""
    for(i in rev(unlist(img)))
    {
      fn = str_extract(i,'(?<=src=")[^"]+')
      if(!is.na(fn) && !startsWith(fn,'img'))
      {
        cap = str_extract(i, '(?<=title=")[^"]+')
        j = gsub('src="[^"]+"',sprintf('src="data:image/png;base64,%s"', fig[[basename(fn)]]),i)
        # figure captions are lost in the conversion
        if(!is.na(cap) && cap!= pcap && !startsWith(cap,'plot of chunk'))
        {
          j = sprintf('%s<p class="caption">%s</p>',j,cap)
          pcap=cap
        }
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

# is it safe to cut after line ...
blurb_safe_cut = function(lines)
{
  tag_open = c()
  in_math = FALSE
  in_code = FALSE
  safe = logical(length(lines))
  lines = trimws(lines)
  for(i in seq_along(lines))
  {
    l = lines[i]
    if(l=='$$')
    {
      in_math = !in_math
    } else if(startsWith(l,'```'))
    {
      in_code = !in_code
    } else if(!startsWith(l,'<img src="data:') && !startsWith(l,'<!') && !in_code && !in_math)
    {
      tags = str_extract_all(l,'<[^>]+>')[[1]]
      tags_start = grepl('^<\\w', tags)
      tags = str_extract(tags, '\\w+(?= |>)') 
      
      for(j in seq_along(tags))
      {
        if(tags_start[j]) tag_open = c(tags[j],tag_open)
        else if(tags[j]==tag_open[1]) tag_open = tag_open[-1]
        else stop('invalid html')
      }
    }
    safe[i] = !in_code && !in_math && length(tag_open)==0
  }
  safe
}

blurb_info = function(fns, n_lines=10)
{
  lapply(fns, function(fn)
  {
    lines = readLines(file.path('vignettes','blog',fn))
    lines = lines[nchar(trimws(lines))>0]
    start = which(grepl('^---',lines))[2] + 1
    
    if(lines[start] == '<style>')
      start = which(endsWith(lines,'</style>'))[1]+1
    
    pre_amble = lines[1:(start-1)]
    
    if(startsWith(lines[start],'#')) start = start + 1
    
    lines = lines[- (1:(start-1))]
    
    end = min(n_lines, length(lines))
    if(end != length(lines))
    {
      end = which(blurb_safe_cut(lines) & 1:length(lines) >=end)[1]
    }

    tibble(
      filename = fn,
      blurb = paste(lines[1:end],collapse='\n'),
      abbreviated = end != length(lines),
      author = gsub('"','', gsub(' * ',', ', trimws(str_extract(pre_amble[startsWith(pre_amble,'author:')][1],'(?<=:).+')) , fixed=TRUE)),
      title = trimws(gsub('"','',str_extract(pre_amble[startsWith(pre_amble,'title:')][1],'(?<=:).+'))),
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
    mutate(rdmore = if_else(abbreviated,sprintf('[read more...](%s)',href),'')) %>%
    mutate(txt = sprintf('\n## [%s](%s)\n<p class="blog-authors">%s</p><p>\n%s\n</p>%s',
                         title,href,author,blurb,rdmore)) %>%
    arrange(desc(filename)) 
  
  content2 = content |>
    filter(abbreviated) |>
    mutate(txt = sprintf('<li><a href="%s">%s</a></li>',href,title)) |>
    pull(txt) |>
    paste(collapse='\n')

  minimal = sprintf('<div class="roll-minimized"><ul>%s</ul></div>',content2)
  
  style = 'img{max-width:8cm;max-height:8cm;display:block;}
  #refs{display:none;}\nsmall.dont-index{display:none;}
  p.blog-authors{font-style:italic}
  div.roll-minimized{
    position:absolute;
    width:35%;
    left:106%;
  }
  /*div.roll-minimized ul{list-style:none;}*/
  div.roll-minimized li{margin-bottom:6px;}
  @media (max-width: 992px)
  {
    div.roll-minimized{display:none;}
  }
  '
  pream = 'title: Dexterities\nbibliography: dexter.bib'

  
  cat(sprintf('---\n%s\n---\n\n<style>\n%s\n</style>\n%s\n%s', pream, style, minimal, paste(content$txt,collapse='\n')),
      file='vignettes/blog/index.Rmd')
}



