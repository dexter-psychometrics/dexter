

# attempt to translate a quoted predicate and an environment to an SQL 'WHERE' statement
#
qtpredicate2where = function(qtpredicate, db, env)
{
  if(is.sql(qtpredicate))
  {
    qtpredicate = trimws(qtpredicate)
    
    if(qtpredicate == '')
      return(qtpredicate)
    
    if(startsWith(qtpredicate,'WHERE '))
      return(qtpredicate)
    
    return(paste('WHERE', qtpredicate))
  }
  
  vars = unique(c(dbListFields(db,'dxItems'), dbListFields(db,'dxBooklets'), dbListFields(db,'dxPersons'), 
                  dbListFields(db,'dxBooklet_design'),dbListFields(db,'dxResponses'),
                  dbListFields(db,'dxScoring_rules'), dbListFields(db,'dxAdministrations')))
  
  variant = switch(class(db), SQLiteConnection = 'sqlite', PostgreSQLConnection = 'postgresql', 'ansi')
  
  pred = partial_eval(qtpredicate, vars=vars, env=env)
  
  if(is.logical(pred) && isTRUE(pred))
    return('')
  
  if(is.null(pred))
    stop('cannot be translated')
  
  sql = translate_sql(pred, variant = variant)

  return(paste(' WHERE ', sql))
}


# evaluates quoted expression as used in _get_responses
# to see if it is safe to trust the booklet_id column from the database
#  
# We can be sure the booklets are not mutilated if
# no item or respons level columns are used in the expression.
# This can err on the safe side but never on the fast but unsafe side.
is_bkl_safe = function(dataSrc, qtpredicate)
{
  if(inherits(dataSrc,'data.frame')) 
    return(FALSE)
  
  if(is.null(qtpredicate)) 
    return(TRUE)
 
  db = dataSrc
  
  blacklist = unique(c( dbListFields(db,'dxItems'),
                        dbListFields(db,'dxScoring_rules'),
                        dbListFields(db,'dxBooklet_design'),
                        dbListFields(db,'dxResponses'))) 
  
  blacklist = setdiff(blacklist, c('person_id','booklet_id'))

  if(is.sql(qtpredicate))
  {
    # to do:dplyr has an sql class, support it
    pred_vars = attr(qtpredicate,'vars')
  } else
  {
    pred_vars = all.vars(qtpredicate)
  }
  
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

## replacement of dbplyr, testing

# local, !!, !!!
# list/df
# functions with vars as arguments
# evaluation
# user defined functions
# assignments
# type sql
# dplyr shit, is_namespaced_dplyr_call
# factors
# not allow formula's I think, check (anonymous function syntax with filter?)
# %>%
# in lang eval:
#if (is_namespaced_dplyr_call(fun)) {
#  call[[1]] <- fun[[3]]
# nested lists
# in as.character
# translate sql 1:4, etc
# .data
# named function arguments
# function calls without arguments


# to do: make unit tests

# testf = function()
# {
#   library(rlang)
#   library(devtools)
#   load_all()
#   
#   df = data.frame(a=rnorm(10),b=rep('a',10))
#   
#   vars=c('a','b','c','d','e')
#   env = current_env()
#   tst = quote(a!=df$a & !(b %in% df$b) & r==TRUE) 
#   partial_eval(tst, vars = vars, env)
#   
#   t2=quote(a %in% x)
#   x=df$b # factor
#   
#   partial_eval(t2)
#   
#   t3 = quote(a==3 && b %in% as.character(c('1','2','3')) )
#   
#   partial_eval(t3,'',env)
#   
#   t3 = quote(a==3 && b %in% as.integer(c('1','2','3')) )
#   
#   partial_eval(t3,c('a','b'),env)
#   
#   t4 = quote(a==3 && b %in% as.character(c(1,3,a)) )
#   
#   pev = partial_eval(t4,c('a','b'),env)
#
#   cat(translate_sql(pev,'sqlite'))
#   
#   t5 = quote(grepl(a,b, !(b %in% df$b)))
#   
#   partial_eval(t5)
#   
#   t6 = quote(y %in% (df %>% mutate(a=a-1) %>% pull(a)))
#   
#   partial_eval(t6,vars='y')
#   
#   
#   t7 = quote(a == 1L)
#   
#   
#   translate_sql(quote(!(b %in% c(1L, 1L, 1L, 1L))))
#   
#   a=list(1,2,b=list(a=1))
#   partial_eval(quote(p < a$b$a),vars='p',env)
#   partial_eval(quote(p < a$b$a),vars=c('p','a'),env)
#
# }

eval_symbol = function(sbl, vars, env)
{
  name = as.character(sbl)
  
  if(name %in% vars)
    return(sbl)
  # functions should possibly/probably be excluded here
  if (env_has(env, name, inherit = TRUE)) 
    return(eval(sbl, env))
  
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
  if(is.function(call[[1]]))
  {
    return(eval(call,env))
  }
  if(is.call(call[[1]]))
  {
    return(eval(call,env))
  }
  
  # to do: local has to be tested
  if(length(intersect(c('%like%'), all.names(call))) == 0 && 
     length(intersect(vars, non.nse.vars(call))) == 0)
  {
    out = try(eval(call, env), silent=TRUE)
    if(!inherits(out,'try-error'))
      return(out)
  }
  
  name = as.character(call[[1]])
  
  if(name=='local')
    return(eval(call[[2]], env))
  
  if (name %in% c("$", "[[", "["))
    return(eval(call, env))
  
  if (name == "remote")
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

sql_quote = function (x, quote) 
{
  if (length(x) == 0) 
    return(x)
  
  y = gsub(quote, paste0(quote, quote), x, fixed = TRUE)
  y = paste0(quote, y, quote)
  y[is.na(x)] = "NULL"
  y
}

# to do: maybe not necessary, check how db engines react to in((1,2,3))
sql_in_bracket = function(txt)
{
  if(startsWith(txt,'(') && endsWith(txt,')'))
  {
    if( !grepl('\\)|\\(', substring(txt,2,nchar(txt)-1),perl=TRUE))
      return(txt)
    if(regexpr('(',substring(txt,2),fixed=TRUE)[1] < regexpr(')',substring(txt,2),fixed=TRUE)[1])
      return(txt)
  }  
  
  paste0('(',txt,')')
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
# to do: {}
# to do, research db varname quotes
translate_sql_lang = function(call, variant)
{
  name = tolower(as.character(call[[1]]))
  
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
  
  if(name == '%in%')
  {
    # for floats this might not be entirely correct
    if(as.character(call[[3]][[1]]) == ':')
      return(paste(translate_sql(call[[2]], variant), 'BETWEEN', 
                   translate_sql(call[[3]][[2]],variant), 'AND', translate_sql(call[[3]][[3]],variant)))

    return(paste(translate_sql(call[[2]], variant), 'IN', sql_in_bracket(translate_sql(call[[3]], variant))))
  }
  # to do: pass on infixes other than %>%, %*%, etc, e.g. %\w{2,}%  
  
  # to do: quote if not quoted
  if(name == '%like%')
    return(sql_infix(call, 'LIKE', variant))
  
  if(name %in% c('(','{') )
    return(paste0('(',translate_sql(call[[2]], variant),')'))
  
  #to do: xor
  
  if(name %in% c('&&','&'))
    return(paste0('(',sql_infix(call, 'AND', variant),')'))
  
  if(name %in% c('||','|'))
    return(paste0('(',sql_infix(call, 'OR', variant),')'))
  
  if(name == '!')
  {
    if(as.character(call[[2]][[1]]) == '!')
      stop('possible quasiquotation, untranslatable') 
    return(paste(' NOT', translate_sql(call[[2]], variant)))
  }
  
  if(startsWith(name,'as.'))
    return(sql_cast(call[[2]], name, variant))
  
  if(name == 'c')
    return(paste0('(',paste(sapply(call[-1], translate_sql, variant=variant), collapse=','),')'))
  
  if(name %in% c('is.na','is.null'))
    return(paste(translate_sql(call[[2]], variant),'IS NULL'))
  
  # to do: extra brackets?
  if(name == 'between')
    return(paste(translate_sql(call[[2]], variant), 'BETWEEN',
                 translate_sql(call[[3]], variant), 'AND',
                 translate_sql(call[[4]], variant)))
  
  if(name %in% c('toupper','tolower'))
    name = toupper(gsub('to','',name, fixed=TRUE))
  
  if(name == ':')
    stop('untranslatable') # to do: maybe in postgress
  
  if(variant == 'postgresql')
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
  }
  
  # to do: grepl, substring and some other functions have sql variants in some db's
  #rtrim,ltrim, abs,coalesce are official sql
  #print(name)
  #in sqlite 
  # pmax=max,pmin=min (multi-arg), min/max should not be allowed unless contain vector
  # instr
  # like()
  # nullif, round, replace, soundex,substr,typeof
  # switch? between?
  
  #I guess a function
  return(paste0(name,'(',paste(sapply(call[2:length(call)], translate_sql, variant=variant), collapse=','),')'))
  
}




# also vars?

translate_sql = function(e, variant) # variant = c('ansi','sqlite','postgresql')
{
  type = typeof(e)
  
  # to do: test what happens with backticks
  # symbol quote dependent on db
  if(type == 'symbol')
    return(paste0('"',as.character(e),'"'))  
  
  if(type == 'language')
    return(translate_sql_lang(e, variant))
  
  if(type == 'list')
    stop('no lists')
  
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
  
  if(type == 'NULL' || is.na(e))
    return('NULL')

  if(is.factor(e))
  {
    e = as.character(e)
    type = "character"
  }
  
  if(type == "character")
    e = sql_quote(e,"'")
  
  # to do: all, any?
  
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
    
    return(paste('(',paste(out,collapse=','),')'))
  }
  
  # e.g. complex not natively supported anywhere, will give errors
  stop(type)
}


