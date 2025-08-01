# common utility functions


.onLoad = function(libname, pkgname) 
{
  dexter_options = list(dexter.use_tibble=FALSE, dexter.progress=TRUE)
  
  to_set = setdiff(names(dexter_options), names(options()))
  if(length(to_set)>0)
    options(dexter_options[to_set])

  invisible()
}

get_ncores = function(desired=256L, maintain_free=0L)
{
  # relevant discussion: https://github.com/Rdatatable/data.table/issues/5658
  
  min_not_NA = function(...)
  {
    x = c(...)
    x = x[!is.na(x)]
    if(length(x)==0) NA_integer_
    else min(x)
  }
  
  if(Sys.getenv("_R_CHECK_LIMIT_CORES_") != "" || 
     Sys.getenv("_R_CHECK_EXAMPLE_TIMING_CPU_TO_ELAPSED_THRESHOLD_") != "")
  {
    # we are on cran
    as.integer(min(desired, 2L))
  } else
  {
    n_cores = omp_ncores()
    
    user_maximum = if(!is.null(getOption('dexter.max_cores'))){
      getOption('dexter.max_cores')
    } else
    {
      min_not_NA(as.integer(Sys.getenv("OMP_THREAD_LIMIT")), getOption("Ncpus"))
    }
    
    if(is.na(user_maximum))
      user_maximum = n_cores - maintain_free
    
    available = min(c(n_cores, user_maximum), 
                    na.rm=TRUE)
    
    as.integer(min(desired, max(1L, available )))
  }
}



explicit_NA = function(x, replace_NA_with = c('<NA>','.NA.','<<NA>>','__NA__'))
{
  if(!is.character(x) || !anyNA(x))
    return(x)
  
  for(v in replace_NA_with)
  {
    if(!v %in% x)
    {
      x[is.na(x)] = v
      return(x)
    }
  }
  stop('could not resolve NA values')
  
}


get_datatype_info = function(dataSrc, columns)
{
  if(is.null(dataSrc)) return(NULL)
  
  res = if(is_db(dataSrc) || is.matrix(dataSrc))
  {
    lapply(columns, function(colname) list(type='character'))
  } else
  {
    columns = intersect(columns,colnames(dataSrc))
    
    lapply(columns,function(colname)
    {
      if(is.factor(dataSrc[[colname]]))
        list(type='factor',levels=levels(dataSrc[[colname]]))
      else
        list(type=class(dataSrc[[colname]]))
    })
  }
  names(res) = columns
  res
}


df_format = function(df, datatype_info=NULL)
{
  if(getOption('dexter.use_tibble', FALSE))
  {
    df = as_tibble(ungroup(df))
  } else
  {
    df = as.data.frame(df, stringsAsFactors = FALSE)
  }
  
  compatible_factor_levels = function(vec, levels)
  {
    # factoring should not lead to NA values
    if(is.factor(vec))
      x = levels(vec)
    else
      x = unique(as.character(vec))
    length(setdiff(x, levels)) == 0
  }
  
  vec_equal = function(a,b) length(a)==length(b) && all(a==b)
  
  for(name in names(datatype_info))
  {
    if(!name %in% colnames(df))
      next
    
    if(datatype_info[[name]]$type == 'factor')
    {
      if(is.factor(df[[name]]) && vec_equal(levels(df[[name]]), datatype_info[[name]]$levels))
        next
      
      if(compatible_factor_levels(df[[name]], datatype_info[[name]]$levels))
        df[[name]] = factor(df[[name]], levels = datatype_info[[name]]$levels)
      else 
        df[[name]] = as.character(df[[name]])
    } else
    {  
      df[[name]] = as.character(df[[name]])
      if(datatype_info[[name]]$type == 'integer')
      {
        df[[name]] = as.integer(df[[name]])
      } else if(datatype_info[[name]]$type == 'numeric')
      {
        
        df[[name]] = as.numeric(df[[name]])
      } 
    }
  }
  
  for(name in setdiff(colnames(df),names(datatype_info)))
  {
    if(is.factor(df[[name]]))
      df[[name]] = as.character(df[[name]])
  }
  
  df
}


is.date = function(x) inherits(x, "Date")
is.time = function(x) inherits(x,'POSIXt')



# x is a vector, str should include one %s to place the collapsed vector
# example: str = "%s [is a/are] reserved variable name[s] in a dexter project"
# first (optional) element in [] is singular, second element is plural
format_plural = function(str, x, sep=', ', last_sep=' and ', qt="'")
{
  len = length(x)
  str = gsub(x = str, pattern='\\[(?:([^\\]/]+)/)?([^\\]/]+)\\]', 
    replacement = paste0('\\', 1+(len>1)),  perl=TRUE)
  
  x = paste0(qt,x,qt)
  
  if(len>=2)
    x = c(head(x,len-2), paste(tail(x,2), collapse=last_sep))
  
  sprintf(str, paste(x,collapse=sep))
}


# format string with named arguments

# txt: format string with named arguments prefixed by a dollar sign, formatting can be done with postfixing with : 
# arglist: list with named arguments to interpolate in the format string. Use only alphanumerical characters in the names

# return: 
# formatted string

# examples

# fstr('$bla, $b',list(bla='some string'))
# fstr('$bla:.1f, $b',list(bla=3.2423))
#
fstr = function(txt, arglist)
{
  if(length(txt) != 1 || length(arglist)==0 || is.null(txt) || !grepl('$',txt,fixed=TRUE)) 
    return(txt)

  txt = gsub('%','%%',txt,fixed=TRUE)
  
  spr_args = list()
  txt_copy = txt
  matches = gregexpr('\\$[a-z]\\w*(\\:[\\.\\d]*[ifscega])?', txt, ignore.case = TRUE, perl=TRUE)[[1]]
  matches_end = attr(matches, 'match.length') + matches - 1L
  for(i in seq_along(matches))
  {
    matched_str = substring(txt, matches[i], matches_end[i])
    matched_tpl = unlist(strsplit(matched_str,':',fixed=TRUE))
    matched_name = substring(matched_tpl[1],2)
    # text can come after a match if match not ends in :
    if(length(matched_tpl) == 1 && !(matched_name %in% names(arglist)))
    {
      l = -1
      m = NULL
      for(nm in names(arglist))
      {
        if(startsWith(matched_name, nm) && nchar(nm) > l)
        {
          l = nchar(nm)
          m = nm
        }
      }
      if(!is.null(m))
      {
        matched_name = m
        matched_str = paste0('$',m)
      }
    }
    
    if(matched_name %in% names(arglist))
    {
      arg = arglist[[matched_name]]
      if(length(matched_tpl) == 2)
      {
        txt_copy = sub(matched_str, paste0('%',matched_tpl[2]) ,txt_copy, fixed=TRUE)
      } else
      {
        arg = as.character(arg)
        txt_copy = sub(matched_str, '%s' ,txt_copy, fixed=TRUE)
      }
      spr_args = append(spr_args, arg) 
    }
  }
  
  do.call(sprintf, append(txt_copy, spr_args))
}






# non vectorized version of ifelse
if.else = function(test, yes, no)
{
  if(isTRUE(test)) return(yes)
  no
}


#  basic argument type and attribute checks with error messages
stop_ = function(...) stop(..., call. = FALSE)

check_dataSrc = function(x)
{
  force(x)
  if(inherits(x, 'dx_resp_data'))
     return(NULL)
  
  if(is.matrix(x))
  {
    if(!mode(x) %in% c('numeric','integer'))
      stop_('dataSrc must be a matrix of positive numbers, a data.frame or a database connection')
    return(NULL)
  } 
  
  if(inherits(x, 'data.frame'))
  {
    if(length(setdiff(c('person_id','item_id','item_score'), colnames(x)))>0)
      stop_("dataSrc must contain the columns: 'person_id','item_id','item_score'")
    return(NULL)
  }
  if(is_db(x))
  {
    if(dbIsValid(x)) return(NULL)
    stop_('your database connection is no longer valid, you need to reconnect. See: ?open_project for details')
  }
    
     
  if(length(x)== 1 && is.character(x) && file.exists(x))
      stop_('dataSrc is a string but must be a database connection, data.frame or matrix. ',
           'Did you forget to do: `db = open_project("',x,'")`?')
 
  stop_("dataSrc must be of type 'DBIconnection', 'data.frame' or 'matrix'")
}

check_file = function(x, name = deparse(substitute(x)))
{
  if(!length(x)== 1 && is.character(x))
    stop_(paste(name, 'must be a string'))
  
  if(!file.exists(x))
    stop_(sprintf('file "%s" does not exist', x))
}

check_db = function(x)
{
  if(length(x)== 1 && is.character(x) && file.exists(x))
  {
      stop_('db is a string but must be a database connection. ',
            'Did you forget to do: `db = open_project("',x,'")`?')
  } else if(!is_db(x))
  {
    stop_("Argument 'db' is not a database connection.")
  } else if(!dbIsValid(x))
  {
    stop_('your database connection is no longer valid, you need to reconnect. see: ?open_project for details')
  }
  
  NULL
}


check_vector = function(x, type=c('character','numeric','integer'), name = deparse(substitute(x)), 
                        nullable = FALSE, .length = NULL, .min=NULL, .max=NULL, min_length=NULL )
{
  if(nullable && is.null(x))
    return(NULL)
  
  type = match.arg(type)
  
  if(!is.vector(x) || is.list(x))
    stop("Argument '",name, "' must be a vector of type ", type, call.=FALSE)
  
  if(type != 'character')
  {
    if(!is.numeric(x))
      stop_("Argument '",name, "' must be ", type)
  
    if(type=='integer' && !(is.integer(x) || all(x %% 1 == 0)))
       stop_("Argument '",name, "' must be an integer")
    
    if(!is.null(.min) && any(x<.min))
      stop_("Argument '",name, "' must be >= ", .min)
    
    if(!is.null(.max) && any(x>.max))
      stop_("Argument '",name, "' must be <= ", .max)
  }
  

  if(!is.null(.length) && length(x) != .length )  
    stop_("Argument '",name, "' must have length ", .length)
  
  if(is.null(length) && !is.null(min_length) && length(x)<min_length)
    stop_("Argument '",name, "' must have minimum length ", min_length)
 
  NULL
}
  
  
check_num = function(x, type=c('numeric','integer'), name = deparse(substitute(x)), 
                     nullable = FALSE, .length = NULL, .min=NULL, .max=NULL, min_length=1 )
{
  type=match.arg(type)
  name=force(name)
  check_vector(x,type=type,name=name,nullable=nullable,.length=.length,.min=.min,.max=.max,min_length=min_length)
}


check_character = function(x, name = deparse(substitute(x)), nullable = FALSE, .length = NULL,min_length=1)
{
  name=force(name)
  check_vector(x,type='character',name=name,nullable=nullable,.length=.length,min_length=min_length)
}

check_string = function(x, name = deparse(substitute(x)), nullable = FALSE)
{
  if(nullable && is.null(x))
    return(NULL)
  
  if(!is.character(x) && length(x)!=1)
    stop_("Argument '",name, "' must be a string of length 1")
}


check_df = function(x, columns=NULL, n_rows=NULL, name = deparse(substitute(x)), nullable=FALSE)
{
  if(nullable && is.null(x))
    return(NULL)
  
  if(!inherits(x, 'data.frame'))
    stop_("Argument'",name, "' must be a data.frame")
  
  missing_col = setdiff(columns, colnames(x))
  
  if(!is.null(columns) && length(missing_col>0))
    stop_(format_plural(paste('Column[s]: %s must be present in', name), missing_col, qt='`'))
  
  if(!is.null(n_rows) && NROW(x)!=n_rows)
    stop_('argument`', name, '` must have ', n_rows, ' rows')
  
}

check_list = function(x, name = deparse(substitute(x)), nullable=FALSE)
{
  if(nullable && is.null(x))
    return(NULL)
  if(!inherits(x,'list'))
    stop_(name,' must be a list')
}


# start, stop
# rng_fl(1:6) => c(1,6)
# rng_fl(c(1,6)) => c(1,6)
# rng_fl(c(6,1)) => error, etc.
rng_fl = function(x, name = deparse(substitute(x)))
{
  if(length(x)==2)
  {
    if(x[1]>x[2])
      stop_('first element of ',name,' must be smaller or equal than second element')
    x
  } else if(length(x)<2)
  {
    stop_(name, ' must be a vector of length 2')
  } else
  {
    test = x[1]:x[length(x)] 
    if(length(x) != length(test) || !all(test == x) || x[1]>x[length(x)])
      stop_(name, ' is not a valid range')
    c(x[1],x[length(x)])
  }
}




dropNulls = function (x) 
{
  x[!vapply(x, is.null, FUN.VALUE = logical(1))]
}

# use for forwarding arguments to e.g. plot function
merge_arglists = function(args, default = NULL, override = NULL)
{
  if(!is.null(default))
    args = modifyList(default, args)
  
  if(!is.null(override))
    args = modifyList(args, override)

  args
}


df_identical = function(a, b)
{
  # check all values in dataframe equal, disregard column order
  
  if(!all(dim(a)==dim(b))) return(FALSE)
  if(!length(intersect(colnames(a),colnames(b))) == ncol(a)) return(FALSE)
  
  a = a |> mutate_if(is.factor, as.character) 
  b = b |> mutate_if(is.factor, as.character)
  
  return(all(a == b[,colnames(a)]))
}

is_connected = function(design)
{
  design = droplevels(design)
  # usually best via items but with extreme predicates and large data this might lead to 
  # adj larger than working memory
  if(nlevels(design$item_id) >= nlevels(design$booklet_id))
  {
    items = as.matrix(table(design$item_id, design$booklet_id))
    adj = crossprod(items, items)
  } else
  {
    booklets = as.matrix(table(design$booklet_id, design$item_id))
    adj = crossprod(booklets, booklets)
  }
    
  mode(adj) = 'integer'
  
  all(ds_connected_groups(adj)==1)
}
