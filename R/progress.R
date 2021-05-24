
global_prog = new.env(parent = emptyenv(), size=2L)
global_prog$bar = NULL
global_prog$bar_caller = NULL


# check if it makes sense to show a progress bar
show_progress = function()
{
  getOption("dexter.progress", TRUE) && interactive() && !getOption("knitr.in.progress", FALSE)
}


# get prog bar with highlander principle
get_prog_bar = function(nsteps=NULL, retrieve_data=FALSE, prog_show = show_progress())
{
  if(is.null(global_prog$bar))
  {
    pb = if(is.null(nsteps)) prog_bar(retrieve_data, prog_show)
         else perc_bar(nsteps, retrieve_data, prog_show)
    global_prog$bar = pb
    global_prog$bar_caller = as.character(sys.call(-1)[[1]])
    
  } else if(!is.null(nsteps))
  {
    global_prog$bar$set_nsteps(nsteps)
  }
  global_prog$bar
}

# closes current prog bar only if it is called from the same function that opened it
close_prog_bar = function()
{
  if(!is.null(global_prog$bar_caller) && global_prog$bar_caller == as.character(sys.call(-1)[[1]]))
  {
    get_prog_bar()$close()
    global_prog$bar = NULL
    global_prog$bar_caller = NULL
  }
}

# base prog bar where the number of iterations is unknowable
prog_bar = setRefClass('prog_bar',
  fields = list(pshow='logical', dx_prog='logical', dplyr_prog='logical', data_step='logical'),
  methods = list(
    close = function()
    {
      options(ddplyr.progress = dplyr_prog)
      if(pshow) cat('\n')
    },
    initialize = function(retrieve_data=FALSE, prog_show = show_progress())
    {
      pshow <<- prog_show
      dplyr_prog <<- as.logical(options(dplyr.show_progress=FALSE))
      data_step <<- as.logical(retrieve_data)
      if(pshow)
      {
        if(data_step)
        {
          cat('| retrieving data...')  
        } else
        {
          cat('|')
        }
      }
    },
    tick = function(nticks=1L)
    {
      if(pshow)
      {
        if(data_step)
        {
          data_step <<- FALSE
          cat(paste0('\r|',strrep(' ',20L)))
          cat('\r|')
        }
        cat(rep('=',nticks))
      }
    }
  ))


# prog bar where the number of iterations is known in advance
# can also be used as an hierarchical progress bar, max one level of depth
perc_bar = setRefClass('perc_bar', contains="prog_bar",
  fields = list(nsteps='integer', step='integer', w='integer', l='integer',p='integer',
                sub_bar='logical', sub_step='integer', sub_nsteps='integer', sub_length='integer',
                step_before_sub='integer'),
  methods = list(
    initialize = function(nsteps_=50L, retrieve_data=FALSE, prog_show=show_progress())
    {
      callSuper(retrieve_data, prog_show)
      nsteps <<- as.integer(nsteps_)
      p <<- -1L
      l <<- -1L
      step <<- 0L
      sub_bar <<- FALSE
      if(pshow)
      {
        w <<- getOption("width") - nchar('| 100%') - 2L
        if(!data_step)
        {
          draw_perc()
        }
      } else
      {
        w <<- -1L
      }
    },
    draw_perc = function()
    {
      old = l+p
      if(step > nsteps) step <<- nsteps
      p <<- as.integer(100L*step/nsteps)
      l <<- as.integer(w * step/nsteps)
      if(pshow && old != l+p)
      {
        cat(sprintf('\r|%s%s| %3i%%',strrep('=',l),strrep(' ',w-l),p))
      }
    },
    close_sub_bar = function()
    {
      fill()
      sub_bar <<- FALSE
    },
    open_sub_bar = function(length_=10L, nsteps_=50L)
    {
      if(sub_bar)
        close_sub_bar()
      sub_bar <<- TRUE
      step_before_sub <<- step
      sub_length <<- as.integer(length_)
      sub_nsteps <<- as.integer(nsteps_)
      sub_step <<- 0L
    },
    tick = function(nticks=1L)
    {
      if(sub_bar)
      {
        sub_step <<- sub_step + as.integer(nticks)
        sub_prog = pmin(sub_step/sub_nsteps,1)
        step <<- step_before_sub + as.integer(1e-6 +  sub_prog * sub_length)
      } else
      {
        step <<- step + as.integer(nticks)
      }
      data_step <<- FALSE
      draw_perc()
    }, 
    set_nsteps = function(nsteps_)
    {
      if(sub_bar)
      {
        sub_nsteps <<- as.integer(nsteps_)
      } else
      {
        nsteps <<- as.integer(nsteps_)
      }
      draw_perc()
    },
    fill = function()
    {
      if(sub_bar)
      {
        sub_step <<- sub_nsteps
        step <<- step_before_sub + sub_length
      } else
      {
        step <<- nsteps
      }
      draw_perc()
    },
    cpp_prog_init = function()
    {
      if(sub_bar)
        c(step_before_sub, sub_nsteps, sub_length, nsteps, w)
      else
        c(0L, nsteps, nsteps, nsteps, w)
    }
  ))
