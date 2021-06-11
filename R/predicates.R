
# attempt to translate a quoted predicate and an environment to an SQL 'WHERE' statement
# also returns vars intersected with db and all_vars
qtpredicate_to_sql = function(qtpredicate, db, env)
{
  if(is.null(qtpredicate))
    return(list(success=TRUE, where='', db_vars = character(), all_vars=character()))

  dbvars = unique(c(dbListFields(db,'dxitems'), dbListFields(db,'dxbooklets'), 
                    dbListFields(db,'dxpersons'), dbListFields(db,'dxbooklet_design'),
                    dbListFields(db,'dxresponses'), dbListFields(db,'dxscoring_rules'), 
                    dbListFields(db,'dxadministrations')))
  
  variant = switch(class(db), SQLiteConnection = 'sqlite', PostgreSQLConnection = 'postgresql', 
                   PqConnection = 'postgresql','ansi')
  
  out = list(success=TRUE)
  
  if(is.sql(qtpredicate))
  {
    out$all_vars = sql_vars(qtpredicate)
    
    qtpredicate = trimws(qtpredicate)
    
    if(qtpredicate == '')
    {
      out$where = ''
    } else
    {
      if(!startsWith(qtpredicate,'WHERE '))
        qtpredicate = paste('WHERE', qtpredicate)
      out$where = qtpredicate
    }
  } else
  {
    pred = try(partial_eval(qtpredicate, vars=dbvars, env=env), silent=TRUE)
    # error is extremely unlikely
    if(is.null(pred) || inherits(pred,'try-error'))
    {
      out$success = FALSE
      out$where = NULL
      out$all_vars = all.vars(qtpredicate)
    } else if(is.logical(pred) && isTRUE(pred))
    {
      out$where = ''
      out$all_vars = character()
    } else
    {
      out$all_vars = all.vars(pred)
      out$where = try(paste(' WHERE', translate_sql(pred, variant = variant)), silent=TRUE)
      if(inherits(out$where, 'try-error'))
        out$success = FALSE
    }
  }
  out$db_vars = intersect(dbvars, out$all_vars)

  for( v in setdiff(out$all_vars, out$db_vars))
    if(!env_has(env, v, inherit=TRUE))
      stop("object '", v, "' not found")
  
  out
}


# evaluates quoted expression as used in get_responses
# to see if it is safe to trust the booklet_id column from the database
#  
# We can be sure the booklets are not mutilated if
# no item or respons level columns are used in the expression.
# This can err on the safe side but never on the fast but unsafe side.
is_bkl_safe = function(dataSrc, qtpredicate, env)
{
  if(!is_db(dataSrc)) 
    return(FALSE)
  
  if(is.null(qtpredicate)) 
    return(TRUE)
 
  pred_vars = qtpredicate_to_sql(qtpredicate, dataSrc, env)$db_vars
  
  blacklist = unique(c( dbListFields(dataSrc,'dxitems'),
                        dbListFields(dataSrc,'dxscoring_rules'),
                        dbListFields(dataSrc,'dxbooklet_design'),
                        dbListFields(dataSrc,'dxresponses'))) 
  
  blacklist = setdiff(blacklist, c('person_id','booklet_id'))

  return(length(intersect(pred_vars, blacklist)) == 0 )
}

sql = function(txt, vars = character())
{
  stopifnot(is.character(txt) && length(txt) == 1)
  if(!is.sql(txt))
    class(txt) = append(class(txt),'sql')
  attr(txt,'vars') = vars
  txt
}

is.sql = function(obj) ('sql' %in% class(obj)) 

# if dbvars is null then only works if vars are sql-quoted
sql_vars = function(sql,dbvars=NULL)
{
  if(!is.null(attr(sql,'vars')))
  {
    attr(sql,'vars')
  } else if(!is.null(dbvars))
  {
    dbvars[sapply(dbvars, function(v){
      grepl(paste0('(^|[^\\w\\d])',v,'($|[^\\w\\d])'),sql,perl=TRUE,ignore.case=TRUE)})]
  } else
  {
    qt = '"'
    if(grepl('`',sql,fixed=TRUE))
      qt = '`'
    
    xpr = paste0('\\',qt,'[^\\',qt,']+','\\', qt)
    gsub(qt, '', unlist(regmatches(sql, gregexpr(xpr, sql, perl=TRUE))))
  }
}


eval_symbol = function(sbl, vars, env)
{
  name = as.character(sbl)
  
  if(name %in% vars)
    return(sbl)
  # functions should possibly/probably be excluded here
  if (env_has(env, name, inherit = TRUE)) 
    return(eval(sbl, env))
  
  if(tolower(name) %in% vars)
  {
    message("All variable names in a project are lowercase, changed '", name, "' -> '", tolower(name),"'")
    return(as.symbol(tolower(name)))
  }
  
  sbl
}



non.nse.vars = function(e)
{
  type = typeof(e)
  
  if(type == 'symbol')
    return(as.character(e))
  
  if(type == 'language')
  {
    name = as.character(e[[1]])
    if(name %in% c('filter','arrange','mutate','group_by','summarize',
                   'summarise','subset','select','rename','pull'))
      return(character())
    
    if(length(e) == 1)
      return(character())
    
    return(unlist(sapply(e[-1], non.nse.vars, simplify = TRUE)))
  }
  
  if(type == 'list')
    return(unlist(sapply(e, non.nse.vars, simplify = TRUE)))
  
  character()
}

eval_lang = function(call, vars, env)
{
  if(is.function(call[[1]]) || is.call(call[[1]]))
  {
    return(eval(call,env))
  }

  if(length(intersect(c('%like%','get'), all.names(call))) == 0 && 
     length(intersect(vars, non.nse.vars(call))) == 0)
  {
    out = try(eval(call, env), silent=TRUE)
    if(!inherits(out,'try-error'))
      return(out)
  }
  
  name = as.character(call[[1]])
  
  if(name=='local')
    return(eval(call[[2]], env))
  
  if(name == 'get')
  {
    a = partial_eval(call[[2]],vars=vars,env=env)
    if(is.character(a) && length(a)==1 && a %in% vars)
      return(as.name(a))
  }
  
  # support dplyrs nasty/quasi quotation, I much prefer local but ok
  if(name=='!' && startsWith(paste(as.character(call),collapse=''),'!!'))
  {
    if(typeof(call[[2]][[2]]) == 'symbol')
    {
      return(eval(call[[2]][[2]], env))
    } else if(typeof(call[[2]][[2]]) == 'language')
    {
      call[[2]][[2]][[2]] = eval(call[[2]][[2]][[2]], env)
      return(eval_lang(call[[2]][[2]], vars=vars, env=env))
    }
  }
    
  if(name %in% c("$", "[[", "["))
    return(eval(call, env))
  
  if(name == "remote")
    return(call[[2]])
  
  call[-1] = lapply(call[-1], partial_eval, vars = vars, env = env)
  

  if(substring(name,1,3)=='as.')
  {
    if(all(all.names(call) %in% c('c',name)) && length(all.vars(call)) == 0)
      return(eval(call, env))
  }
  
  call
}

# partial_eval is meant for pre-evaluating statements to be turned into sql predicates
# do not use it for other purposes
partial_eval = function (e, vars = character(), env) 
{
  type = typeof(e)
  
  if(type == 'symbol')
    return(eval_symbol(e, vars, env))

  if(type == 'language')
    return(eval_lang(e, vars, env))
  
  if(type == 'list')
    return(lapply(e, partial_eval, vars=vars, env=env))
  
  e
}

check_function_call = function(call)
{
  f = sapply(formals(as.character(call[[1]])), is.symbol)
  
  if('...' %in% names(f))
    return(NULL)
    
  arg_names = names(call)[2:length(call)]
  stopifnot(length(call)-1<=length(f))
  
  if(is.null(arg_names))
  {
    if(length(call)-1<length(f))
      stopifnot(all(!f[length(call):length(f)]))
  } else
  {
    stopifnot(length(setdiff(arg_names,c(names(f),"")))==0) 
    stopifnot(length(setdiff(names(f)[f], arg_names)) <= sum(arg_names==""))
  }
}

get_arg = function(name, call)
{
  f = formals(as.character(call[[1]]))
  
  # support positional
  if(is.numeric(name))
    name = names(f)[name]
  
  pos = which(names(f) == name)
  stopifnot(pos>=1)
  
  arg_names = names(call)[2:length(call)]
  
  if(is.null(arg_names))
    return(call[[pos+1]])
  
  if(name %in% arg_names)
    return(call[[name]])
  
 
  i=j=1
  for(a in arg_names)
  {
    if(a=="" && i==pos)
      return(call[[j+1]])
    
    if(a=="" || which(names(f)==a) < pos)
      i = i+1
    j = j+1
  }
  if(is.symbol(f[[name]]))
    stop('argument ', name, ' is missing with no default')
  return(f[[name]])
}


sql_quote = function (x, quote) 
{
  if (length(x) == 0) 
    return(x)
  
  y = gsub(quote, paste0(quote, quote), x, fixed = TRUE)
  y = paste0(quote, y, quote)
  y[is.na(x)] = "NULL"
  y
}

# only for length 1 strings
trim_brackets = function(text)
{
  sapply(text, function(x)
  {
    if(startsWith(x,'(') && endsWith(x,')'))
      return(substr(x,2,nchar(x)-1))
    x
  })
}



sql_infix = function(e, op, variant)
{
  paste(translate_sql(e[[2]], variant), op, translate_sql(e[[3]],variant))
}

sql_cast = function(e, type, variant)
{
  if(typeof(e) == 'language' && as.character(e[[1]]) == 'c')
  {
    if(length(e)>2)
    {
      return(paste0('(',
                    paste('CAST(',sapply(e[2:length(e)], translate_sql, variant=variant),
                          gsub('as.','AS ', type, fixed=TRUE),')', collapse=','),')'))
    }
    e = e[[2]]
  }
  return(paste('CAST(', translate_sql(e, variant=variant), gsub('as.','AS ', type, fixed=TRUE),')'))
}


translate_sql_lang = function(call, variant)
{
  name = as.character(call[[1]])
  
  if(name == '=')
    stop('assignment')
  
  if(name == '==')
    name = '='
  
  if(grepl('^[\\<\\>\\=]\\=?$',name,perl=TRUE))
    return(sql_infix(call, name, variant))
  
  if(name == '!=')
    return(sql_infix(call, '<>', variant))
  
  if(name == '-' && length(call) == 2)
    return(paste0('-',translate_sql(call[[2]], variant)))
  
  if(name %in% c('+','-','/','*'))
    return(sql_infix(call, name, variant))
  
  if(name == '%in%' && typeof(call[[3]]) == 'language' && as.character(call[[3]][[1]]) == ':')
    return(paste('CAST(', translate_sql(call[[2]], variant), ' AS INTEGER) BETWEEN', 
                   translate_sql(call[[3]][[2]],variant), 'AND', translate_sql(call[[3]][[3]],variant)))

  if(name %in% c('(','{') )
    return(paste0('(',translate_sql(call[[2]], variant),')'))
  

  if(name %in% c('&&','&'))
    return(paste0('(',sql_infix(call, 'AND', variant),')'))
  
  if(name %in% c('||','|'))
    return(paste0('(',sql_infix(call, 'OR', variant),')'))
  
  if(name == '!')
  {
    if(as.character(call[[2]][[1]]) == '!')
      stop('possible quasiquotation, untranslatable') 
    return(paste(' NOT (', translate_sql(call[[2]], variant),')'))
  }
  
  if(is.function(call[[1]]))
    check_function_call(call)
  
  if(name == 'xor')
  {
    a = translate_sql(call[[2]], variant)
    b = translate_sql(call[[3]], variant)
    return(paste(a,'OR',b,'AND NOT (',a,'AND',b,')'))
  }
  
  if(startsWith(name,'as.'))
    return(sql_cast(call[[2]], name, variant))
  
  if(name == 'c')
  {
    return(paste0('(',
              paste0(trim_brackets(sapply(call[-1], translate_sql, variant=variant)), collapse=','),
            ')'))
  }
  if(name %in% c('is.na','is.null'))
    return(paste(translate_sql(call[[2]], variant),'IS NULL'))
  
  if(name == 'between')
    return(paste(translate_sql(get_arg(1,call), variant), 'BETWEEN',
                 translate_sql(get_arg(2,call), variant), 'AND',
                 translate_sql(get_arg(3,call), variant)))
  
  if(name %in% c('toupper','tolower'))
    name = toupper(gsub('to','',name, fixed=TRUE))
  
  # simple %in% : already done 
  if(name == ':')
    stop('untranslatable') 
  
  if(name %in% c('nchar','str_length'))
    return(paste0("LENGTH(",translate_sql(call[[2]], variant),")"))
  
  if(name == '%in%')
  {
    type3 = typeof(call[[3]])
    type3_is_vec = type3 %in% c('integer','double','character','logical','Date','POSIXct','POSIXlt')
    if(type3 == 'symbol' ||
       (type3_is_vec && length(call[[3]]) == 1))
    {
      return(paste(translate_sql(call[[2]], variant), 'IN(', translate_sql(call[[3]],variant),')'))
    }
  }
  
  if(startsWith(name,'%') && endsWith(name,'%'))
    return(sql_infix(call, gsub('%','',name,fixed=TRUE), variant))
  
  
  if(name == 'substr')
  {
    start = translate_sql(get_arg('start',call), variant)
    stop = translate_sql(get_arg('stop',call), variant)
    n = if.else(grepl('^\\d+$', start) && grepl('^\\d+$', stop),
                as.integer(stop) - as.integer(start) + 1L,
                paste0('(1+(', stop, ')-(', start, '))'))
    
    x = translate_sql(get_arg('x',call), variant)
    if(variant=='sqlite')
        return(paste("substr(",x, ',', start, ',', n, ')'))
    
    return(paste("substring(",x, 'from', start, 'for', n, ')'))
  }
  
  
  if(variant %in% c('postgresql','ansi'))
  {
    if(name == 'pmax')
      name = 'GREATEST'
    if(name == 'pmin')
      name = 'LEAST'
    if(name == '%%')
      name = 'MOD'
    if(name == '^')
      name = 'POW'
    if(name == 'paste0')
      name = 'CONCAT'
    
    if(name=="paste")
    {
      sep = sql_quote(if.else(is.null(call$sep), " ",call$sep),"'")
      return(paste0('CONCAT_WS(',sep,',',paste(sapply(call[2:length(call)], translate_sql, variant=variant), collapse=','),')'))
    }
    
    if(name == 'startsWith')
      x = translate_sql(get_arg('x',call), variant)
      prefix = translate_sql(get_arg('prefix',call), variant)
      return(paste0("substring(",x," from 1 for char_length(",prefix,"))=",prefix))
      
    if(name == 'endsWith')
      x = translate_sql(get_arg('x',call), variant)
      prefix = translate_sql(get_arg('prefix',call), variant)
      return(paste0("substring(",x,
                    " from (char_length(",x,")-char_length(",prefix,")))=",prefix))
    

    
    # let grepl be sorted out in R, implementation may differ
  }
  if(variant == 'sqlite')
  {
    if(name == 'startsWith')
      return(paste0("substr(",translate_sql(call[[2]], variant),",1,length(",translate_sql(call[[3]], variant),
                    "))=",translate_sql(call[[3]], variant)))
    if(name == 'endsWith')
      return(paste0("substr(",translate_sql(call[[2]], variant),
                    ",length(",translate_sql(call[[2]], variant),")-length(",translate_sql(call[[3]], variant),
                    "))=",translate_sql(call[[3]], variant)))
    
    if(name=="paste")
    {
      sep = sql_quote(" ","'")
      if(!is.null(call$sep))
        sep = translate_sql(call$sep, variant)
      return(paste(sapply(call[2:length(call)], translate_sql, variant=variant), collapse=paste0('||',sep,'||')))
    }
    if(name=="paste0")
    {
      return(paste(sapply(call[2:length(call)], translate_sql, variant=variant), collapse='||'))
    }
  }
  
  if(name=='~')
    stop("no sql translation for formula's")
  
  # rtrim,ltrim, abs, coalesce are official sql, possibly more
  return(paste0(name,'(',paste(sapply(call[2:length(call)], translate_sql, variant=variant), collapse=','),')'))
  
}


translate_sql = function(e, variant) # variant = c('ansi','sqlite','postgresql')
{
  type = typeof(e)
  
  # symbol quote dependent on db seems not necessary, although sqlite supports `
  if(type == 'symbol')
    return(paste0('"',as.character(e),'"'))  
  
  if(type == 'language')
    return(translate_sql_lang(e, variant))
  
  if(type == 'list')
    stop('lists cannot be translated')
  
  if(class(e) == 'matrix')
  {
    e = drop(e)
    stopifnot(is.null(dim(e)))
  }
  if(class(e) == 'array')
  {
    stopifnot(length(dim(e))==1)
    e = as.vector(e)
  }
  
  if(type == 'NULL')
    return('NULL')

  if(is.factor(e))
  {
    e = as.character(e)
    type = "character"
  }
  
  if(type == "character")
    e = sql_quote(e,"'")
  
  if(type == 'Date')
  {
    e = sql_quote(format(e, "%Y-%m-%d"),"'")
    if(variant == 'postgres')
      e = paste('date', e)
  }
  
  if(type %in% c('POSIXct','POSIXlt'))
  {
    e = sql_quote(format(e),"'")
    if(variant == 'postgres')
      e = paste('timestamp', e)
  }
  if(type == 'logical' && variant == 'sqlite')
    e = as.integer(e)
  
  if(type %in% c('integer','double','character','logical','Date','POSIXct','POSIXlt'))
  {
    if(length(e) == 0)
      stop('vector of length 0')
    
    out = as.character(e)
    out[is.na(e)] = 'NULL'
    
    if(type == 'double' && any(is.infinite(e)))
    {
      if(variant == 'sqlite')
        stop('sqlite and inf do not go well together')
      out[is.infinite(e) & e<0] = "'-Infinity'"
      out[is.infinite(e) & e>0] = "'Infinity'"
    }
    
    if(length(out) == 1)
      return(out)
    

    return(paste0('(',paste(out,collapse=','),')'))
  }
  
  # e.g. complex not natively supported anywhere, will give errors
  stop(type)
}





correct_symbol_case = function (e, vars = character(), env) 
{
  type = typeof(e)
  
  if(type == 'symbol')
  {
    name = as.character(e)
    if(name %in% vars || env_has(env, name, inherit = TRUE))
      return(e)
    
    if(tolower(name) %in% vars)
      return(as.symbol(tolower(name)))
    
    return(e)
  }

  if(length(e)>1)
    e[] = lapply(e, correct_symbol_case, vars=vars, env=env)
  
  e
}


