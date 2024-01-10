
global_prog = new.env(parent = emptyenv(), size=1L)
global_prog$bar = NULL


# check if it makes sense to show a progress bar
show_progress = function()
{
  getOption("dexter.progress", TRUE) && interactive() && !getOption("knitr.in.progress", FALSE)
}


# get prog bar with highlander principle
get_prog_bar = function(nsteps=NULL, retrieve_data=FALSE, prog_show = show_progress(), 
                        type = if.else(is.null(nsteps),'indeterminate','determinate'),
                        lock=FALSE)
{
  caller = paste(as.character(sys.call(-1)[[1]]),collapse='_')
  type = match.arg(type)
  if(is.null(global_prog$bar))
  {
    global_prog$bar = progress_bar(nsteps=nsteps, retrieve_data=retrieve_data,
                                   prog_show=prog_show, type=type, caller=caller, lock=lock)

  } else if(!is.null(nsteps))
  {
    global_prog$bar$set_nsteps(nsteps, caller=caller)
  }
  global_prog$bar
}


# base prog bar where the number of iterations is unknowable
progress_bar = setRefClass('prog_bar',
  fields = list(pshow='logical', dplyr_prog='logical', data_step='logical',determinate='logical',
                pcaller='character', areas='list',w='integer',locked='logical'),
  methods = list(
    close = function()
    {
      caller = paste(as.character(sys.call(-1)[[1]]),collapse='_')
      if(caller == areas[[1]]$caller)
      {
        areas <<- areas[1]
        fill()
        options(dplyr.show_progress = dplyr_prog)
        if(pshow) cat('\n')
        global_prog$bar = NULL
      } else if(!locked)
      {
        a = which(sapply(areas,'[[', 'caller') == caller)
        if(length(a)>0)
        {
          areas <<- areas[1:a[1]]
          close_area()
        }
      }
    },
    initialize = function(caller, type=c('determinate','indeterminate'), nsteps=50L, 
                          retrieve_data=FALSE, prog_show = show_progress(), lock=FALSE)
    {
      pshow <<- prog_show
      dplyr_prog <<- as.logical(options(dplyr.show_progress=FALSE))
      data_step <<- as.logical(retrieve_data)
      locked <<- as.logical(lock)
      pcaller <<- as.character(caller)
      
      determinate <<- match.arg(type) == 'determinate'
      areas <<-list()
      areas[[1]] <<- list(nsteps=nsteps, step=0L, caller=caller)

      if(pshow)
      {
        w <<- getOption("width") - nchar('| 100%') - 2L
        if(data_step)
        {
          cat('| retrieving data...')  
        } else
        {
          if(determinate) draw_perc()
          else cat('|')
        }
      } else
      {
        w <<- -1L
      }
    },
    permission = function()
    {
      !locked || pcaller == paste(as.character(sys.call(-2)[[1]]),collapse='_')
    },
    lock = function()
    {
      if(pcaller == paste(as.character(sys.call(-1)[[1]]),collapse='_'))
        locked <<- TRUE
    },
    unlock = function()
    {
      if(pcaller == paste(as.character(sys.call(-1)[[1]]),collapse='_'))
        locked <<- FALSE
    },
    close_area = function()
    {
      if(permission())
      {
        a = length(areas)
        if(a>1)
        {
          areas[[a-1]]$step <<- areas[[a-1]]$step + areas[[a]]$steps_allocated 
          areas <<- areas[1:(a-1)]
          if(determinate) draw_perc()
        }
      }
    },
    new_area = function(steps_allocated, nsteps=10L)
    {
      if(permission())
      {
        caller = paste(as.character(sys.call(-1)[[1]]),collapse='_')
        if(caller == areas[[length(areas)]]$caller)
          close_area()
        areas[[length(areas)+1]] <<- list(nsteps=nsteps, step=0L, caller=caller, steps_allocated = steps_allocated)
      }
    },
    set_nsteps = function(nsteps, caller=NULL)
    {
      if(permission())
      {
        areas[[length(areas)]]$nsteps <<- nsteps
        draw_perc()
      }
    },
    draw_perc = function()
    {
      if(pshow && determinate)
      {
        stw = 1/areas[[1]]$nsteps
        p=0
        for(a in areas)
        {
          
          if(!is.null(a$steps_allocated))
            stw = stw * a$steps_allocated/a$nsteps
          p = p + stw * a$step
        }
        p = p+1e-10
        prc = as.integer(100 * p)
        l =  as.integer(w * p)
        cat(sprintf('\r|%s%s| %3i%%',strrep('=',l),strrep(' ',w-l),prc))
      }
    },
    fill = function()
    {
      a=length(areas)
      areas[[a]]$step <<- areas[[a]]$nsteps
      draw_perc()
    },
    tick = function(nticks=1L)
    {
      if(pshow && permission())
      {
        if(data_step)
        {
          data_step <<- FALSE
          cat(paste('\r',strrep(' ',25L)))
          cat('\r|')
        }
        if(determinate)
        {
          a=length(areas)
          areas[[a]]$step <<- pmin(areas[[a]]$step + nticks, areas[[a]]$nsteps)
          draw_perc()
        }
        else cat(rep('=',nticks))
      }
    },
    cpp_prog_init = function()
    {
      a = length(areas)
      v = if(a==1)
      {
        c(areas[[1]]$step, areas[[1]]$nsteps - areas[[1]]$step, areas[[1]]$nsteps - areas[[1]]$step,
          areas[[1]]$nsteps, w)
      } else
      {
        p=1
        s=0
        for(ar in areas[2:a])
        {
          p = p * ar$steps_allocated/ar$nsteps
          s = s + p * ar$steps
        }
        c(s, areas[[a]]$steps, p*areas[[1]]$nsteps, areas[[1]]$nsteps, w)
      }
      as.integer(v)
    }
  ))

