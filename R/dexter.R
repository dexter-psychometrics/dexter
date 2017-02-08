################################
make_ReST = function(db) {
  q = "create view if not exists ReST as
  select s.booklet, s.item, s.person, s.response, coalesce(r.score,0) as score
  from responses as s left outer join rules as r using(item,response)"
  junk = dbSendQuery(db,q)
  dbClearResult(junk)
  return(NULL)
}

################################
testlev = function(x) {
  x = x[complete.cases(x), ]
  N = max(x$n)
  k = nrow(x)
  if (k>1) {
    alpha = k/(k-1) * ( 1 - sum(x$sdScore^2) / sum(x$rit * x$sdScore)^2 )
  } else alpha = NA
  meanP = mean(x$pvalue)
  meanRit = mean(x$rit)
  meanRir = mean(x$rir)
  maxScore = sum(x$maxScore)
  res = c(k,alpha,meanP,meanRit,meanRir,maxScore,N)
  names(res) = c("nItems","alpha","avgPval","avgRit","avgRir","maxScore","N")
  res
}

############################
#void ElSym(double *b, int *a, int *first, int *last, int *item1, int *item2, int *nI, int *mS, double *g, int *idx);
elsym = function(b,a,first,last)
{
  n=length(first)
  ms=sum(a[last])
  g=double((ms+1))
  tmp<-.C("ElSym",
          as.double(b),
          as.integer(a),
          as.integer(first-1),
          as.integer(last-1),
          as.integer(-1),
          as.integer(-1),
          as.integer(n),
          as.integer(ms),
          as.double(g),
          as.integer(0))
  return(tmp[[9]])
}


###############################
ittot = function(b,c,a,first,last,it)
{
  pi=matrix(0,last[it]-first[it]+1,sum(a[last])+1)
  logb=log(b); logc=log(c)
  for (s in 0:(sum(a[last])))
  {
    eta=exp(logb+(a*s)*logc)
    g=elsym(eta,a,first,last)
    gi=elsym(eta,a,first[-it],last[-it])
    k=1
    for (j in first[it]:last[it])
    {
      idx=s+1-a[j]
      if ((idx>0)&(idx<=length(gi))) { pi[k,s+1]=exp(log(eta[j])+log(gi[idx])-log(g[s+1])) }
      k=k+1
    }
  }
  return(pi)
}


###################################
EstIM  = function(ss) {
  first = ss$il$first
  last = ss$il$last
  a = ss$sl$score
  sufI = ss$sl$sufI
  sufC = ss$il$sufC
  C = rep(1:nrow(ss$il), ss$il$nCat)
  scoretab = ss$tl$N
  
  nI=length(last)
  b=rep(0,length(sufI))
  ic=rep(1,length(sufC))
  se.ic=vector("numeric", nI)
  HRM=matrix(0,length(b),length(b))
  
  # Identification
  b[sufI>0]=1
  local_upd_set=vector("list",nI)
  global_upd_set=vector("list",nI)
  for (i in 1:nI)
  {
    local_upd_set[[i]]=which(sufI[first[i]:last[i]]>0)
    local_upd_set[[i]]=local_upd_set[[i]][-1]
    global_upd_set[[i]]=(first[i]:last[i])[local_upd_set[[i]]]
  }
  
  converged=2
  while(converged>0.001)
  {
    converged=-1
    for (i in 1:nI)
    {
      if (length(local_upd_set[[i]])>0)
      {
        pi=ittot(b,ic[C],a,first,last,i)
        pi[is.na(pi)]=0
        pi=pi[local_upd_set[[i]],,drop=FALSE]
        E=sufI[global_upd_set[[i]]]-pi%*%scoretab
        H=-pi%*%diag(scoretab)%*%t(pi)
        diag(H)=pi%*%scoretab+diag(H)
        
        # NR update for parameters of item i
        update=solve(H,E)
        b[global_upd_set[[i]]]=b[global_upd_set[[i]]]*exp(update)
        converged=max(converged,max(abs(E)))
        HRM[global_upd_set[[i]],global_upd_set[[i]]]=H
      }
    }
  }
  
  bRM=b
  cRM=ic
  ## IM
  
  b[sufI>0]=1
  
  converged=2
  while(converged>0.001)
  {
    converged=-1
    for (i in 1:nI)
    {
      # gradient and hessian for thresholds of item i
      if (length(local_upd_set[[i]])>0)
      {
        pi=ittot(b,ic[C],a,first,last,i)
        pi[is.na(pi)]=0; pi=pi[local_upd_set[[i]],,drop=FALSE]
        E=sufI[global_upd_set[[i]]]-pi%*%scoretab
        H=-pi%*%diag(scoretab)%*%t(pi)
        diag(H)=pi%*%scoretab+diag(H)
        
        # gradient and hessian for interaction parameter
        ncol_pi=ncol(pi)
        E=c(E,sufC[i])
        H=cbind(H,rep(0,nrow(H)))
        H=rbind(H,rep(0,ncol(H)))
        k=1
        e0=0; e1=0
        f=matrix(0,nrow(pi),ncol_pi)
        g=matrix(0,nrow(pi),ncol_pi)
        h=0
        for (j in global_upd_set[[i]])
        {
          E[length(E)]=E[length(E)]-a[j]*sum((0:(ncol_pi-1))*scoretab*pi[k,])
          e0=e0+a[j]*pi[k,]
          e1=e1+a[j]^2*pi[k,]
          f[k,]=a[j]*(0:(ncol_pi-1))*pi[k,]
          g[k,]=pi[k,]
          h=h+a[j]*(0:(ncol_pi-1))*pi[k,]
          k=k+1
        }
        H[nrow(H),nrow(H)]=sum((0:(ncol_pi-1))^2*(e1-e0^2)*scoretab)
        for (k in 1:nrow(f))
        {
          H[k,nrow(H)]=sum((f[k,]-g[k,]*h)*scoretab)
          H[nrow(H),k]=H[k,nrow(H)]
        }
        # NR update for parameters of item i
        update=solve(H,E)
        b[global_upd_set[[i]]]=b[global_upd_set[[i]]]*exp(update[-length(update)])
        ic[i]=ic[i]*exp(update[length(update)])
        se.ic[i]=solve(H)[nrow(H),nrow(H)]
        converged=max(converged,max(abs(E)))
      }
    }
  }
  return(list(group=ss$group,bRM=bRM,cRM=cRM,bIM=b,cIM=ic,se.c=se.ic,HRM=HRM))
}


#################################
bty = function (n, h = c(265, 75), c. = c(61, 66),
                l = c(25, 73), power = c(0.7, 1.742),
                fixup = TRUE, gamma = NULL, alpha = 1, ...)
{
  if (!is.null(gamma))
    warning("'gamma' is deprecated and has no effect")
  if (n < 1L)
    return(character(0L))
  h <- rep(h, length.out = 2L)
  c <- rep(c., length.out = 2L)
  l <- rep(l, length.out = 2L)
  power <- rep(power, length.out = 2L)
  rval <- seq(1, 0, length = n)
  rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L],
                       C = c[2L] - diff(c) * rval^power[1L], H = h[2L] - diff(h) *
                         rval), fixup = fixup, ...)
  if (!missing(alpha)) {
    alpha <- pmax(pmin(alpha, 1), 0)
    alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)),
                    width = 2L, upper.case = TRUE)
    rval <- paste(rval, alpha, sep = "")
  }
  return(rval)
}


##############################
my_layout = function(npic, nr, nc) {
  if(npic==1) nr=nc=1
  nc = min(nc, 3)
  nc = min(nc, npic)
  nw = npic %/% nc + npic %% nc
  nr = min(nr, 3)
  nr = min(nr, nw)
  list(nr=nr, nc=nc)
}

####################################
simpleTIA = function(x) {
  data.frame(item = x$item[1],
             nBooklets = nrow(x),
             Pval = weighted.mean(x$pvalue, w=x$n, na.rm=TRUE),
             Rit = weighted.mean(x$rit, w=x$n, na.rm=TRUE),
             Rir = weighted.mean(x$rir, w=x$n, na.rm=TRUE),
             N = sum(x$n, na.rm=TRUE),
             mnScore = weighted.mean(x$meanScore, w=x$n, na.rm=TRUE),
             sdScore = weighted.mean(x$sdScore, w=x$n, na.rm=TRUE)
  )
}

##################################
## as with estim.. call with c=ic[C]
SSTable <- function(m, AB, model) {
  if (model=="IM") {ic=m$est$cIM; b=m$est$bIM} else {ic=m$est$cRM; b=m$est$bRM}
  first = m$ss$il$first
  last =  m$ss$il$last
  C = rep(1:nrow(m$ss$il), m$ss$il$nCat)
  a = m$ss$sl$score
  ic = ic[C]
  A = AB[[1]]
  B = AB[[2]]
  ### Check
  if (length(intersect(A,B))!=0) stop("sets not disjunct")
  if (length(union(A,B))!=length(b[first])) warning("Sets A and B may not include all items")
  
  ## Bookkeeping
  bA=NULL; bB=NULL
  aA=NULL; aB=NULL
  cA=NULL; cB=NULL
  lastA=NULL; firstA=NULL
  lastB=NULL; firstB=NULL
  telAF=1; telBF=1
  for (i in 1:length(first))
  {
    if (i %in% A)
    {
      bA=c(bA,b[first[i]:last[i]])
      aA=c(aA,a[first[i]:last[i]])
      cA=c(cA,ic[first[i]:last[i]])
      firstA=c(firstA,telAF)
      lastA=c(lastA,telAF+last[i]-first[i])
      telAF=telAF+last[i]-first[i]+1
    }
    if (i %in% B)
    {
      bB=c(bB,b[first[i]:last[i]])
      aB=c(aB,a[first[i]:last[i]])
      cB=c(cB,ic[first[i]:last[i]])
      firstB=c(firstB,telBF)
      lastB=c(lastB,telBF+last[i]-first[i])
      telBF=telBF+last[i]-first[i]+1
    }
  }
  MscA=sum(aA[lastA])
  MscB=sum(aB[lastB])
  Msc=sum(a[last])
  out=matrix(NA,MscA+1,MscB+1)
  
  ### Rasch Model
  if (model=="RM")
  {
    gA = elsym(bA,aA,firstA,lastA)
    gB = elsym(bB,aB,firstB,lastB)
    g =  elsym(b,a,first,last)
    for (s in 0:Msc)
    {
      for (s_a in max(0,s-MscB):min(s,MscA))
      {
        s_b=s-s_a
        out[s_a+1,s_b+1]=log(gA[s_a+1])+log(gB[s_b+1])-log(g[s+1])
      }
    }
  }
  
  ### IM
  if (model=="IM")
  {
    logb=log(b); logc=log(ic)
    logbA=log(bA); logcA=log(cA)
    logbB=log(bB); logcB=log(cB)
    for (s in 0:Msc)
    {
      eta=exp(logb+(a*s)*logc)
      etaA=exp(logbA+(aA*s)*logcA)
      etaB=exp(logbB+(aB*s)*logcB)
      gA = elsym(etaA,aA,firstA,lastA)
      gB = elsym(etaB,aB,firstB,lastB)
      g =  elsym(eta,a,first,last)
      for (s_a in max(0,s-MscB):min(s,MscA))
      {
        s_b=s-s_a
        out[s_a+1,s_b+1]=log(gA[s_a+1])+log(gB[s_b+1])-log(g[s+1])
      }
    }
  }
  return(list(tbl=exp(out),m=m,AB=AB,model=model))
}



############################################
######      exported functions        ######
############################################


#' Start a new project
#'
#' Imports a complete set of scoring rules and starts a new project (data base)
#'
#'
#' @param rules A data frame with columns \code{item}, \code{response}, and \code{score}.
#' The order is not important but spelling is. Any other columns will be ignored.
#' @param db_name A name for the data base that will be created. If the name does not
#' contain a path, the file will be created in the work
#' directory. Any existing file with the same name will be overwritten.
#' @return If the scoring rules pass a sanity check, a handle to the data base.
#' Otherwise, a data frame listing the problems found.
#' @author Ivailo Partchev
#' @details This package only works with closed items: it does not score any open items.
#' The first step to creating a project is to import an exhaustive list of all items and
#' all admissible responses, along with the score that any of the latter will be given.
#' Responses may be numbers or strings, and they must appear exactly as in the data
#' that will be imported. Scores must be integers, and the minimum score for an item
#' must be 0. When inputting data, all responses not on the list will be treated as
#' missing and ultimately scored 0, but it is good style to include the missing
#' responses in the list.
#' importFrom RSQLite dbConnect dbReadTable dbWriteTable dbListFields
#' importFrom RSQLite dbGetQuery dbSendQuery
#' @export
#'
start_new_project = function(rules, db_name="dexter.db") {
  rules = rules[, c("item", "response", "score")]
  sanity = ddply(rules, ~item, function(x) {
    problem = character(0)
    if (length(unique(x$score))<2) problem = paste(problem, "Need at least two distinct scores")
    if (any(duplicated(x$response))) problem = paste(problem, "Cannot have duplicated responses")
    if (min(x$score) != 0) problem = paste(problem, "The minimum score must be 0")
    if (length(problem)) y = data.frame(problem=problem) else y = data.frame(problem=character(0))
    y
  })
  if (nrow(sanity)) {
    cat("There were problems with your scoring rules.\nCheck the output for possible reasons.\n")
    return(sanity)
  } else {
    if (file.exists(db_name)) file.remove(db_name)
    db = dbConnect(SQLite(), db_name)
    dbWriteTable(db, "rules", rules)
    return(db)
  }
}





#' Open an existing project
#'
#' Opens a data base created by function \code{start_new_project}
#'
#'
#' @param db_name The name of the data base to be opened.
#' @return A handle to the data base.
#' @author Ivailo Partchev
#' @export
#'
open_project = function(db_name="dexter.db") {
  if (file.exists(db_name)) {
    db = dbConnect(SQLite(), db_name)
  } else stop("There is no such file")
  return(db)
}





#' Add a booklet to a project
#'
#' Adds item response data for a test form (a.k.a. booklet)
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @param x A data frame containing the responses and, possibly, some additional
#' person characteristics. See details.
#' @param booklet_id A (short) string identifying the test form (booklet)
#' @return A list of: \item{items}{The names of the columns in \code{x} that were
#' treated as items}
#' \item{covariates}{The names of the columns in \code{x} that were
#' treated as person covariates}
#' \item{not_listed}{A data frame of all responses that will be treated as missing}
#' @author Ivailo Partchev
#' @details It is common practice to keep data in rectangular tables: data frames
#' or foreign software like Excel, SPSS, etc. This function is provided to input
#' data in that form, one booklet at a time. The starting point is a data frame,
#' and getting the data frame into R is left to the user. We have found package
#' \code{readxl} to be very good at reading Excel sheets, and \code{haven} quite
#' efficient with SPSS files.
#'
#' This package is not doing any person management. We assume that each person
#' has responded to only one test form, once, and we don't check for any
#' duplicates. If users supply a sound person ID of their own, they will be able to
#' link ability measures with results of the same person on a different test, but
#' this is their own responsibility. Also, data from two-stage tests (a routing test
#' and a follow-up test) should be matched by the user before entered.
#'
#' Any variable whose name has an exact match in the scoring rules input with
#' function \code{start_new_project} will be treated as an item; any other
#' variables will be treated as person covariates. Any responses to an item that
#' do not have an exact match in the scoring rules will be treated as missing, and
#' ultimately given the lowest score of 0. To score missing data differently, or simply
#' abide to good style, the user can include explicit entries for missing value
#' indicators in the scoring rules.
#' @export
#' 
add_booklet = function(db, x, booklet_id) {
  rules = dbReadTable(db, "rules")
  items = unique(rules$item)
  isItem = names(x) %in% items
  itemCols = which(isItem)
  covCols  = which(!isItem)
  if (!dbExistsTable(db, "booklets")) {
    junk=dbSendQuery(db, "create table booklets(bookletName,nItems,nPersons)")
    dbClearResult(junk)
  }
  booklets = dbReadTable(db, "booklets")
  if (booklet_id %in% booklets$bookletName) stop("There is already a booklet with this ID")
  query = paste("insert into booklets VALUES ('",
                booklet_id,"',",length(itemCols),",",nrow(x),")", sep="")
  junk=dbSendQuery(db, query)
  dbClearResult(junk)
  lastBooklet = dbGetQuery(db, "SELECT last_insert_rowid()")[1,1]
  design = data.frame(booklet=lastBooklet,
                      booklet_id=booklet_id,
                      item=names(x)[itemCols],
                      position=1:length(itemCols))
  x$dxBookletID = lastBooklet
  x$dxPersonID  = row.names(x)
  persons = x[,-itemCols]
  responses = melt(x,
                             id.vars = c("dxBookletID","dxPersonID"),
                             measure.vars = itemCols)
  names(responses) = c("booklet", "person", "item", "response")
  scores = merge(responses, rules, by=c("item","response"), all.x=TRUE, all.y=FALSE)
  badCases = scores[is.na(scores$score),]
  badCases$score[is.na(badCases$score)] = 0
  not_on_list = badCases[!duplicated(badCases[,c("item","score")]),]
  persons= melt(persons, id.vars = c("dxBookletID","dxPersonID"))
  dbWriteTable(db, "persons", persons, append=TRUE)
  dbWriteTable(db, "responses",  responses, append=TRUE)
  dbWriteTable(db, "design",  design, append=TRUE)
  return(
    list(
      items = names(x)[itemCols],
      covariates = names(x)[covCols],
      not_on_list = not_on_list[,c("booklet","item","response","score")],
      n_responses_treated_as_NA = nrow(badCases)
    )
  )
}



#' Add item properties to a project
#'
#' Adds item properties to an existing data base
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @param df A data frame containing the item properties. See details.
#' @return A list of: \item{unknown_items}{Item IDs for any items that were provided
#' in the data frame but could not be found in the data base}
#' \item{items_unaccounted_for}{Item IDs for any items that exist in the data base
#' but were not given properties in the data frame}
#' @author Ivailo Partchev
#' @details When entering response data in the form of a rectangular person x item
#' table, it is easy to provide person properties but practically impossible
#' to provide item properties. This function provides a possibility to do so.
#' The order of the rows and columns in the data frame is not important but
#' (i) there must be a column called exactly \code{item} containing the item IDs
#' exactly as entered before, and (ii) all items in the data frame must be known
#' to the data base and all items in the data base must be given properties --
#' otherwise, there will be a warning message, and nothing else will be done.
#' If all is well, the data frame will be added to the data base as table
#' \code{item_properties}, and any variables in it may be used in analyses involving
#' item properties.
#' @export
#'
add_item_properties = function(db, df){
  itemz = dbGetQuery(db, "select distinct item as dxExistingItems,
                              1 as dxMysteriousJunk from design")
  itemz = merge(itemz, df, by.x="dxExistingItems", by.y="item", all.x=TRUE, all.y=TRUE)
  unknown_items = itemz$existing[is.na(itemz$dxMysteriousJunk)]
  wh = which(!(names(itemz) %in% c("dxExistingItems","dxMysteriousJunk")))[1]
  items_unaccounted_for = itemz$dxExistingItems[is.na(itemz[,wh])]
  if (length(unknown_items)+length(items_unaccounted_for)) {
    warning("Found items that do not exist in the data base and/or
            existing items unaccounted for -- please see the output
            and try again")
  } else {
    vn = match(c("dxExistingItems","dxMysteriousJunk"), names(itemz))
    names(itemz)[vn[1]] = "item"
    dbWriteTable(db, "item_properties", itemz[,-vn[2]], append=FALSE, overwrite=TRUE)
  }
  list(unknown_items=unknown_items, items_unaccounted_for=items_unaccounted_for)
}


#' Derive scoring rules from keys
#'
#' For multiple choice items that will be scored as 0/1, derive the
#' scoring rules from the keys to the correct responses
#'
#'
#' @param keys  A data frame containing columns \code{item}, \code{nOptions}, and
#' \code{key} (the spelling is important). See details.
#' @return A data frame that can be used as input to \code{start_new_project}
#' @author Ivailo Partchev
#' @details
#' This function might be useful in setting up the scoring rules when all items
#' are multiple-choice and scored as 0/1. (Hint: Because the order in which the
#' scoring rules is not important, one can use the function to generate rules for
#' many MC items and then append per hand the rules for a few complex items.)
#'
#' The input data frame must contain the exact name of each item, the number
#' of options, and the key. If the keys are all integers, it will be assumed that
#' responses are coded as 1 through nOptions. If they are all uppercase letters,
#' it is assumed that responses are coded as A,B,C,... All other cases result
#' in an error.
#' @export
#'
keys_to_rules = function(keys) {
  if (is.numeric(keys$key)) ABC=FALSE else {
    if (all(keys$key %in% LETTERS)) ABC=TRUE
    else stop("You have inadmissible keys")
  }
  if (ABC) {
    m = match(keys$key, LETTERS)
    if (any(m>keys$nOptions)) stop("You have out-of-range keys")
  } else {
    if (any(keys$key>keys$nOptions)) stop("You have out-of-range keys")
  }
  ddply(keys, ~item, function(x){
    k = x$nOptions[1]
    if (ABC) {
      response=LETTERS[1:k]
      theKey = match(x$key[1], LETTERS)
    } else {
        response=1:k
        theKey = x$key[1]
    }
    y = data.frame(response=response, score=0L)
    y$score[theKey] = 1
    y
  })
}




#' List booklets in a project
#'
#' Show a list of the test forms (booklets) that have been entered in the db
#' so far
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data frame showing the internal booklet number, the booklet name,
#' the number persons and the number of items
#' @author Ivailo Partchev
#' @export
#'
show_booklets = function(db) {
  dbReadTable(db, "booklets")
}



#' List item properties
#'
#' Show a list of the item properties defined in the project (if any)
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return Nothing but messages printed to the screen
#' @author Ivailo Partchev
#' @export
#'
show_item_properties = function(db) {
  if (dbExistsTable(db, "item_properties")) {
      str(dbReadTable(db, "item_properties"))
  } else cat("No item properties defined yet\n")
}


#' List person properties
#'
#' Show a list of the person properties defined in the project (if any)
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return Nothing but messages printed to the screen
#' @author Ivailo Partchev
#' @export
#'
show_person_properties = function(db) {
  dbGetQuery(db, "select distinct variable from persons")
}


#' List items in a project
#'
#' Show a list of all items that have been entered in the db
#' so far by booklet and position in the booklet
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A data frame (in string format) showing the booklet design: rows are items,
#' columns are booklets, and the numbers in the cells show the position of the item in the
#' booklet (blank it the booklet does not include the item)
#' @author Ivailo Partchev
#' @export
#'
show_items = function(db){
  foo = dbReadTable(db, "design")
  bar = format(dcast(foo, item~booklet_id, value.var = "position"))
  bar[bar=="NA"] = ""
  bar
}



#' Estimate the Rasch and the Interaction model
#'
#' Estimate the parameters of the Rasch model and the Interaction model
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @param booklets A vector with the numbers of the test forms (booklets)
#' to be included in the analysis. See Details.
#' @param grpVar Variable name for subsetting
#' @param group Value(s) for subsetting
#' @return An object of class \code{rim} holding results
#' for the Rasch model and the interaction model.
#' @author Gunter Maris, Timo Bechger, Ivailo Partchev
#' @details Unlike the Rasch model, the interaction model cannot be computed
#' concurrently for a whole design of test forms. This function fits the
#' Rasch model and the interaction model on a complete rectangular array of
#' responses, with comparison between the two models playing an important role.
#' The rectangular array may be either one booklet or the intersection (common items)
#' of two or more booklets. If the intersection is empty (no common items), the
#' function will exit with an error message.
#'
#' Please notice that the booklets are specified with their internal numbers
#' in the data base, not with their names. The names and numbers of the booklets
#' can be inspected with function \code{show_booklets}.
#' @useDynLib dexter
#' @export
#'
fit_models = function(db, booklets, grpVar=NULL, group=NULL) {
  make_ReST(db)
  design = dbGetQuery(db, "select distinct booklet, item from design")
  booklets = booklets[booklets %in% design$booklet]
  des = design[design$booklet %in% booklets,]
  if (nrow(des) < 1) stop("No items to analyse")
  ssq = "select * from ReST"
  if (length(booklets)==1) {
    whq = paste("where booklet=", booklets)
    allScores = dbGetQuery(db, paste("select distinct item, score from rules as r
                                              inner join design using(item)",whq,"order by item, score"))
  } else {
    item_list = dlply(des, "booklet", function(x)x$item)
    common = Reduce(intersect, item_list)
    if (length(common)<1) stop("There are no common items")
    itemz = paste0("(",paste(paste("'",common,"'", sep=""), collapse=","),")")
    whq = paste ("where item IN", itemz)
    allScores = dbGetQuery(db, paste("select distinct item, score from rules",
                                              whq,"order by item, score"))
  }
  grpName = "All"
  if (!is.null(grpVar) & !is.null(group)) {
    addGrp = paste0("inner join (select * from persons where variable='",grpVar,
                    "') as g on ReST.booklet=g.dxBookletID and ReST.person=g.dxPersonID")
    ssq = paste(ssq, addGrp)
    if (length(group)==1) {
      whq = paste0(whq, " and value='",group,"' ")
      grpName = group
    } else {
      groupz = paste0("(",paste(paste("'",group,"'", sep=""), collapse=","),")")
      whq = paste0(whq, " and value IN ",groupz)
      grpName = "Subset"
    }
  }
  respData = dbGetQuery(db, paste(ssq, whq))
  respData = ddply(respData, ~booklet+person, function(x){x$sumScore=sum(x$score);x})
  ssScoreLev = ddply(respData, ~item+score, function(x){data.frame(sufI=length(x$score),sufC=sum(x$score*x$sumScore))})
  ssScoreLev = merge(allScores, ssScoreLev, all.x=TRUE)
  ssScoreLev[is.na(ssScoreLev)] = 0
  ssScoreLev=ssScoreLev[order(ssScoreLev$item,ssScoreLev$score),]
  ssItemLev = ddply(ssScoreLev, ~item, function(x){
    data.frame(nCat=length(x$score), N=sum(x$sufI), sufC=sum(x$sufC))
  })
  ssItemLev$last  = cumsum(ssItemLev$nCat)
  ssItemLev$first = ssItemLev$last - ssItemLev$nCat + 1
  plt = ddply(respData, ~item+sumScore, function(x)data.frame(meanScore=mean(x$score),N=length(x$sumScore)))
  maxTotScore = sum(ddply(ssScoreLev, ~item, function(x){data.frame(mx=max(x$score))})$mx)
  totalLev = plt[plt$item==plt$item[1],]
  totalLev = merge(data.frame(sumScore=0:maxTotScore), totalLev, all.x=TRUE)
  totalLev$N[is.na(totalLev$N)] = 0
  ss = list(group=grpName, il=ssItemLev, sl=ssScoreLev, tl=totalLev, plt=plt)
  result = try(EstIM(ss))
  if (inherits(result, "try-error")) result=NULL
  outpt = list(est=result, ss=ss)
  class(outpt) = "rim"
  outpt
}




#' A print method for the interaction model
#'
#' Print the available items for plots of the Rasch and the interaction models
#'
#'
#' @param x An object produced by function \code{fit_models}
#' @param ... Included to stop check from nagging
#' @author Ivailo Partchev
#' @method print rim
#' @export
#'
print.rim = function(x, ...){
  available_items = x$ss$il[,"item"]
  print.default(available_items)
  invisible(x)
}


#' A plot method for the interaction model
#'
#' Plot the item-total regressions fit by the interaction (or Rasch) model
#'
#'
#' @param x An object produced by function \code{fit_models}
#' @param items The items to plot (column numbers). If NULL, all items will be plotted
#' @param summate If FALSE, regressions for polytomous items will be shown for each
#' response option separately; default is TRUE.
#' @param overlay If TRUE and more than one item is specified, there will be two plots,
#' one for the Rasch model and the other for the interaction model, with all items
#' overlayed; otherwise, multiple plots with the two models overlayed. Default is FALSE
#' @param nc An integer between 1 and 3. Number of columns when putting mutiple plots
#' on the same page. Default is 1. May be ignored or adjusted if it does not make sense.
#' @param nr An integer between 1 and 3. Number of rows when putting mutiple plots
#' on the same page. Default is 1. May be ignored or adjusted if it does not make sense.
#' @param curtains 100*the tail probability of the sum scores to be shaded. Default is 10.
#' Set to 0 to have no curtains shown at all.
#' @param show.observed If TRUE, the observed proportion correct at each sum score
#' will be shown as dots. Default is FALSE.
#' @param ... Any additional plotting parameters
#' @author Ivailo Partchev
#' @method plot rim
#' @export
#'
plot.rim = function(x, items=NULL, summate=TRUE, overlay=FALSE,
                    nc=1, nr=1, curtains=10, show.observed=FALSE, ...){
  nit = length(x$est$cIM)
  qua = curtains/200
  if(qua>0 & qua<.5) qnt = quantile(rep(as.integer(x$ss$tl$sumScore),
                                               x$ss$tl$N), c(qua,1-qua)) else qnt=NULL
  maxScore = nrow(x$ss$tl) - 1
  if (is.null(items)) items=1:nit
  npic = length(items)
  if (length(items)==1) nr=nc=1
  if (overlay & !summate) overlay=FALSE
  C = rep(1:nrow(x$ss$il), x$ss$il$nCat)
  
  if (overlay) {
    # only summate possible
    if (nr*nc==2) graphics::layout(matrix(1:2,nr,nc)) else graphics::layout(1)
    # do the Rasch model
    #
    p = matrix(0, npic, maxScore+1)
    r = 0
    for (i in items) {
      itm = x$ss$il[i,]
      opt = x$ss$sl[x$ss$sl$item==itm$item[1],]
      prb = ittot(x$est$bRM, x$est$cRM[C], x$ss$sl$score, x$ss$il$first, x$ss$il$last, i)
      p[r<-r+1,] = crossprod(opt$score, prb)
    }
    graphics::plot(c(0,maxScore),c(0,max(p)),type="n",main="Rasch model",xlab="Test score",
         ylab="Item score")
    if (!is.null(qnt)) {
      tmp = graphics::par('usr')
      graphics::rect(tmp[1], tmp[3], qnt[1], tmp[2], col="#EEEEEE", border=NA)
      graphics::rect(qnt[2], tmp[3], tmp[2], tmp[4], col="#EEEEEE", border=NA)
    }
    for (i in 1:npic) graphics::lines(0:maxScore,p[i,]) # the actual lines
    lx = sample(0:maxScore, npic, replace = FALSE) # label the lines
    for (i in 1:npic) {
      graphics::points(lx[i], p[i,lx[i]+1], co="white", cex=1.6, pch=19)
      graphics::text(lx[i], p[i,lx[i]+1], items[i], co=1, cex=.6)
    }
    # do the Interaction model
    #
    p = matrix(0, npic, maxScore+1)
    r = 0
    for (i in items) {
      itm = x$ss$il[i,]
      opt = x$ss$sl[x$ss$sl$item==itm$item[1],]
      prb = ittot(x$est$bIM, x$est$cIM[C], x$ss$sl$score, x$ss$il$first, x$ss$il$last, i)
      p[r<-r+1,]= crossprod(opt$score, prb)
    }
    graphics::plot(c(0,maxScore),c(0,max(p)),type="n",main="Interaction model",xlab="Test score",
         ylab="Item score")
    if (!is.null(qnt)) {
      tmp = graphics::par('usr')
      graphics::rect(tmp[1], tmp[3], qnt[1], tmp[2], col="#EEEEEE", border=NA)
      graphics::rect(qnt[2], tmp[3], tmp[2], tmp[4], col="#EEEEEE", border=NA)
    }
    for (i in 1:npic) graphics::lines(0:maxScore,p[i,]) # the actual lines
    lx = sample(0:maxScore, npic, replace = FALSE) # label the lines
    for (i in 1:npic) {
      graphics::points(lx[i], p[i,lx[i]+1], co="white", cex=1.6, pch=19)
      graphics::text(lx[i], p[i,lx[i]+1], items[i], co=1, cex=.6)
    }
    graphics::box()
    # end of overlay
  } else {
    # not overlay: do many plots
    ly = my_layout(npic, nr, nc)
    graphics::layout(matrix(1:(ly$nr*ly$nc), byrow=TRUE, ncol=ly$nc))
    if (summate) {
      # for each item in turn, do both models (with summation), and plot
      for (i in items) {
        itm = x$ss$il[i,]
        opt = x$ss$sl[x$ss$sl$item==itm$item[1],]
        prb = ittot(x$est$bIM, x$est$cIM[C], x$ss$sl$score, x$ss$il$first, x$ss$il$last, i)
        prb = crossprod(opt$score, prb)
        graphics::plot(c(0,maxScore), c(0,max(prb)), type="n",
             main=paste("Item", i, ": ", itm$item[1], sep=""),
             xlab="Test score", ylab="Item score")
        if (!is.null(qnt)) {
          tmp = graphics::par('usr')
          graphics::rect(tmp[1], tmp[3], qnt[1], tmp[2], col="#EEEEEE", border=NA)
          graphics::rect(qnt[2], tmp[3], tmp[2], tmp[4], col="#EEEEEE", border=NA)
        }
        if(show.observed) {
          plt = x$ss$plt[x$ss$plt$item==itm$item[1],]
          graphics::points(plt$sumScore,plt$meanScore,col="coral",pch=20)
        }
        graphics::lines(0:maxScore, prb, col="gray80", lwd=3)
        prb = ittot(x$est$bRM, x$est$cRM[C], x$ss$sl$score, x$ss$il$first, x$ss$il$last, i)
        prb = crossprod(opt$score, prb)
        graphics::lines(0:maxScore, prb)
      }
      graphics::box()
    } else {
      # for each item in turn, similar but possibly multiline and coloured
      for (i in items) {
        itm = x$ss$il[i,]
        opt = x$ss$sl[x$ss$sl$item==itm$item[1],]
        prb = ittot(x$est$bIM, x$est$cIM[C], x$ss$sl$score, x$ss$il$first, x$ss$il$last, i)
        # prb = rbind(1-colSums(prb), prb)
        pte = bty(nrow(prb), alpha=.6)
        graphics::plot(c(0,maxScore), 0:1, type="n",
             main=paste("Item", i, ": ", itm$item[1], sep=""),
             xlab="Test score", ylab="Probability")
        if (!is.null(qnt)) {
          tmp = graphics::par('usr')
          graphics::rect(tmp[1], tmp[3], qnt[1], tmp[2], col="#EEEEEE", border=NA)
          graphics::rect(qnt[2], tmp[3], tmp[2], tmp[4], col="#EEEEEE", border=NA)
        }
        for (j in 1:nrow(prb)) {
          graphics::lines(0:maxScore, prb[j,], col=pte[j], lwd=3)
        }
        prb = ittot(x$est$bRM, x$est$cRM[C], x$ss$sl$score, x$ss$il$first, x$ss$il$last, i)
        # prb = rbind(1-colSums(prb), prb)
        pte = bty(nrow(prb))
        for (j in 1:nrow(prb)) {
          graphics::lines(0:maxScore, prb[j,], col=pte[j])
        }
      } # eol items
      graphics::box()
    } # eo not summate
  } # eo not overlay
}







#' Distractor plot
#'
#' Produce a diagnostic distractor plor for an item
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @param item The ID of the item to plot. A separate plot will be produced
#' for each booklet that contains the item, or an error message if the item ID
#' is not known. Each plot contains a non-parametric regression of each possible
#' response on the total score.
#' @param nc An integer between 1 and 3. Number of columns when putting mutiple plots
#' on the same page. Default is 1. May be ignored or adjusted if it does not make sense.
#' @param nr An integer between 1 and 3. Number of rows when putting mutiple plots
#' on the same page. Default is 1. May be ignored or adjusted if it does not make sense.
#' @author Ivailo Partchev
#' @export
#'
distractor_plot = function(db, item, nc=1, nr=1){
  the_item = dbGetQuery(db, paste("select * from design where item='",item,"'",sep=""))
  if(nrow(the_item)<1) stop(paste("Item",item,"unknown"))
  q = paste0("with sc as
    (select s.booklet, s.item, s.person, s.response, coalesce(r.score,0) as score
    from responses as s
    left outer join rules as r using(item,response)),
    ss as (select booklet, person, sum(score) as sumScore from sc 
      group by booklet, person)
    select sc.*, sumScore from sc inner join ss using(booklet, person) 
    where sc.item='",item,"'")
  dataForItem = dbGetQuery(db, q)

  foo = ddply(dataForItem, ~booklet+item+response+sumScore, function(x)data.frame(n=length(x$sumScore)))
  iSt = ddply(dataForItem, ~booklet+item, function(x)data.frame(meanScore=mean(x$score), rit=cor(x$score,x$sumScore)))
  mxSc = dbGetQuery(db, paste0("select max(score) as maxScore from rules where item='",item,"'")) 
  iSt$pvalue = iSt$meanScore / mxSc$maxScore[1]

  q5 = paste("select * from rules where item='",item,"'",sep="")
  scoru = dbGetQuery(db, paste0("select * from rules where item='",item,"'"))
  foo = merge(foo, scoru)
  
  npic = length(unique(foo$booklet))
  ly = my_layout(npic, nr, nc)
  graphics::layout(matrix(1:(ly$nr*ly$nc), byrow=TRUE, ncol=ly$nc))

  d_ply(foo, "booklet", function(y){
    bk = which(the_item$booklet==y$booklet[1])
    tit = sprintf("%s: %d in %s",the_item$item[bk],the_item$position[bk],the_item$booklet_id[bk])
    subtit = sprintf("Pval: %.3f, Rit: %.3f",iSt$pvalue[bk],iSt$rit[bk])
    graphics::plot(c(0,max(y$sumScore)),c(0,1),type="n",main=tit,sub=subtit,
         xlab="Sum score",ylab="Proportion",cex.sub=.7)
    
    bar = ddply(y, ~sumScore, function(x)data.frame(n=sum(x$n)))
    dAll = density(bar$sumScore, n=51, weights=bar$n/sum(bar$n))
    nnn = sum(bar$n)
    k = 1
    lgnd = ddply(y, ~response, function(z) {
      dxi = density(z$sumScore, n=51, weights=z$n/sum(z$n), bw = dAll$bw, 
            from = min(dAll$x), to = max(dAll$x))
      yy = dxi$y / dAll$y * sum(z$n)/nnn
      graphics::lines(dAll$x, yy, co=k<<-k+1, lw=2)
      data.frame(col=k, resp=paste0(z$response[1]," (",z$score[1],")"))
    }) 
    graphics::legend("right",legend=as.character(lgnd$resp),
           lty=1,col=lgnd$col,cex=.7,box.lty=0)
    graphics::box()
  })
}





#' Add test data for an existing booklet
#'
#' Adds item response data for a test form (a.k.a. booklet) that
#' already exists in the data base. This function has been added 
#' as a workaround by user request -- we do not necessarily approve.
#' 
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @param x A data frame containing the responses and, possibly, some additional
#' person characteristics. See details.
#' @param booklet_id A (short) string identifying the test form (booklet)
#' @return A list of: \item{items}{The names of the columns in \code{x} that were
#' treated as items}
#' \item{covariates}{The names of the columns in \code{x} that were
#' treated as person covariates}
#' \item{not_listed}{A data frame of all responses that will be treated as missing}
#' @author Ivailo Partchev
#' @export
#'
add_tests = function(db, x, booklet_id) {
  booklets = dbReadTable(db, "booklets")
  if (!booklet_id %in% booklets$bookletName) stop("There is no such booklet")
  rules = dbReadTable(db, "rules")
  items = unique(rules$item)
  knownItems = dbGetQuery(db, paste0("select item from design where booklet_id='",
                                              booklet_id,"' order by position"))
  if (any(items!=knownItems$item)) {
    print(data.frame(KnownItems=knownItems$item,NewItems=items))
    stop("Not all items match")
  }
  isItem = names(x) %in% items
  itemCols = which(isItem)
  covCols  = which(!isItem)
  
  query = paste("insert into booklets VALUES ('",
                booklet_id,"',",length(itemCols),",",nrow(x),")", sep="")
  junk=dbSendQuery(db,query)
  dbClearResult(junk)
  
  theBooklet = dbGetQuery(db, paste0("select rowid from booklets where bookletname='",
                                              booklet_id,"'"))$rowid[1]
  lastPerson = booklets$nPersons[booklets$bookletName==booklet_id]
  
  x$dxBookletID = theBooklet
  x$dxPersonID  = row.names(x) + lastPerson
  lastPerson = lastPerson + nrow(x)
  junk=dbSendQuery(db, paste0("UPDATE booklets SET nPersons= ", lastPerson,
                                  " WHERE bookletName='",booklet_id,"'"))
  dbClearResult(junk)
  
  persons = x[,-itemCols]
  responses = melt(x,
                             id.vars = c("dxBookletID","dxPersonID"),
                             measure.vars = itemCols)
  names(responses) = c("booklet", "person", "item", "response")
  scores = merge(responses, rules, by=c("item","response"), all.x=TRUE, all.y=FALSE)
  badCases = scores[is.na(scores$score),]
  badCases$score[is.na(badCases$score)] = 0
  not_on_list = badCases[!duplicated(badCases[,c("item","score")]),]
  persons= melt(persons, id.vars = c("dxBookletID","dxPersonID"))
  dbWriteTable(db, "persons", persons, append=TRUE)
  dbWriteTable(db, "responses",  responses, append=TRUE)
  
  return(
    list(
      items = names(x)[itemCols],
      covariates = names(x)[covCols],
      not_on_list = not_on_list[,c("booklet","item","response","score")],
      n_responses_treated_as_NA = nrow(badCases)
    )
  )
}

#' Estimate the Rasch and the Interaction model per domain
#'
#' Estimate the parameters of the Rasch model and the Interaction model
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @param booklets A vector with the numbers of the test forms (booklets)
#' to be included in the analysis. See Details.
#' @param domains The variable name for the item property defining the 
#' domains (subtests)
#' @return An object of class \code{imp} holding results
#' for the Rasch model and the interaction model.
#' @author Gunter Maris, Timo Bechger, Ivailo Partchev
#' @details Unlike the Rasch model, the interaction model cannot be computed
#' concurrently for a whole design of test forms. This function fits the
#' Rasch model and the interaction model on a complete rectagular array of
#' responses, with comparison between the two models playing an important role.
#' The rectangular array may be either one booklet or the intersection (common items)
#' of two or more booklets. If the intersection is empty (no common items), the
#' function will exit with an error message.
#'
#' Please notice that the booklets are specified with their internal numbers
#' in the data base, not with their names. The names and numbers of the booklets
#' can be inspected with functopn \code{show_booklets}.
#' @export
#'
fit_domains = function(db, booklets, domains) {
  make_ReST(db) 
  # do we have such an item property?
  if (!domains %in% dbListFields(db,"item_properties")) stop("Unknown domain variable")
  
  ### exi = dbGetQuery(db, paste0("select exists(select 1 from persons where variable ='",domains,"')"))
  
  design = dbGetQuery(db, "select distinct booklet, item from design")
  booklets = booklets[booklets %in% design$booklet]
  des = design[design$booklet %in% booklets,]
  if (nrow(des) < 1) stop("No items to analyse")
  
  ssq = paste("select booklet, person, sum(score) as score, ", domains, "as domain",
              "from ReST inner join (select item,",domains, "from item_properties) using(item)")
  
  grq = "group by booklet, person, domain order by booklet, person, domain"
  
  if (length(booklets)==1) {
    whq = paste("where booklet=", booklets)
    allScores = dbGetQuery(db, paste0("select distinct item, score, p.",domains," as domain",
                                               " from rules as r",
                                               " inner join design as d using(item) inner join (select item,",domains,
                                               " from item_properties) as p using(item) ", whq," order by item, score"))
  } else {
    item_list = dlply(des, "booklet", function(x)x$item)
    common = Reduce(intersect, item_list)
    if (length(common)<1) stop("There are no common items")
    itemz = paste0("(",paste(paste("'",common,"'", sep=""), collapse=","),")")
    whq = paste ("where item IN", itemz)
    allScores = dbGetQuery(db,
                                    paste("select distinct item, score,",domains,"as domain",
                                          "from rules inner join (select item,",domains,
                                          "from item_properties) using(item)", whq,"order by item, score"))
  }
  
  bizarre = function(x,y){m=outer(x,y,"+");unique(as.vector(m))}
  add = function(x) Reduce("bizarre", x)
  allScores = ddply(allScores, ~domain, function(x){
    data.frame(score=add(dlply(x, ~item, function(y)y$score)))
  })
  domdat = dbGetQuery(db, paste(ssq, whq, grq))
  # need to add sum scores
  domdat = ddply(domdat, ~booklet+person, function(x){x$sumScore=sum(x$score);x})
  ssScoreLev = ddply(domdat, ~domain+score, function(x){
    data.frame(sufI=length(x$score), sufC=sum(x$score*x$sumScore))
  })
  ssScoreLev = merge(allScores, ssScoreLev, all.x=TRUE)
  ssScoreLev[is.na(ssScoreLev)] = 0
  ssItemLev = ddply(ssScoreLev, ~domain, function(x){
    data.frame(nCat=length(x$score), N=sum(x$sufI), sufC=sum(x$sufC))
  })
  ssItemLev$last  = cumsum(ssItemLev$nCat)
  ssItemLev$first = ssItemLev$last - ssItemLev$nCat + 1
  plt = ddply(domdat, ~domain+sumScore, function(x){
    data.frame(meanScore=mean(x$score), N=length(x$score))
  })
  maxTotScore = sum(ddply(ssScoreLev, ~domain, function(x){data.frame(mx=max(x$score))})$mx)
  totalLev = plt[plt$domain==plt$domain[1],]
  totalLev = merge(data.frame(sumScore=0:maxTotScore), totalLev, all.x=TRUE)
  totalLev$N[is.na(totalLev$N)] = 0
  ssItemLev = rename(ssItemLev, c("domain" = "item"))
  ssScoreLev = rename(ssScoreLev, c("domain" = "item"))
  ssScoreLev=ssScoreLev[order(ssScoreLev$item,ssScoreLev$score),]
  totalLev = rename(totalLev, c("domain" = "item"))
  plt = rename(plt, c("domain" = "item"))
  ss = list(group="All", il=ssItemLev, sl=ssScoreLev, tl=totalLev, plt=plt)
  result = try(EstIM(ss))
  if (inherits(result, "try-error")) result=NULL
  outpt = list(est=result, ss=ss)
  class(outpt) = "rim"
  outpt
}


#' Profile plot
#'
#' Profile plot
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @param booklets A vector with the numbers of the test forms (booklets)
#' to be included in the analysis. See Details.
#' @param domains The item property defining the domains
#' @param grpVar Variable name for subsetting
#' @param group Value(s) for subsetting; if NULL, all distinct values
#' will be used.
#' @param model "IM" (default) or "RM"
#' @return Nothing interesting
#' @author Timo Bechger, Gunter Maris, Ivailo Partchev
#' @export
#'
profile_plot = function(db, booklets, domains, grpVar, group=NULL, model="IM") {
  if (model != "IM") model="RM"
  # do we have such an item property?
  if (!domains %in% dbListFields(db,"item_properties")) stop("Unknown domain variable")
  # do we have such a person property?
  exi = dbGetQuery(db, paste0("select exists(select 1 from persons where variable ='",grpVar,"')"))
  if (!exi[1,1]) stop("Unknown grpVar variable")
  # we only support 2 domains
  theDomains = dbGetQuery(db, paste("select distinct",domains," as dm from item_properties"))$dm
  if (length(theDomains) != 2) stop("Sorry, this function only supports 2 domains")
  if (!is.null(group)) {
    theGroups = group
  } else {
    theGroups = dbGetQuery(db, paste0("select distinct value from persons where variable='",grpVar,"'"))$value
  }
  models = lapply(theGroups, function(x)fit_models(db=db, booklets=booklets, grpVar=grpVar, group=x))
  itemz = as.data.frame(models[[1]]$ss$il$item)
  names(itemz) = "item"
  dom = dbGetQuery(db, paste("select item,",domains,"from item_properties"))
  dom = merge(itemz, dom, all.x=T, all.y=F)
  dom$itemN = as.integer(row.names(dom))
  AB = dlply(dom, as.formula(paste0("~",domains)), function(x)x$itemN)
  tt = lapply(models, function(x)SSTable(x, AB=AB, model=model))
  maxA = nrow(tt[[1]]$tbl)-1
  maxB = ncol(tt[[1]]$tbl)-1
  sg = data.frame(k=0:(maxA+maxB))
  graphics::plot(c(0,maxA), c(0,maxB), type="n", asp=1, main="Profile plot", 
       xlab=theDomains[1], ylab=theDomains[2])
  # The timolines
  k = maxA + maxB
  sg$y0 = pmin(maxB,sg$k)
  sg$x0 = sg$k - sg$y0
  sg$x1 = pmin(maxA,sg$k)
  sg$y1 = sg$k - sg$x1
  graphics::segments(sg$x0, sg$y0, sg$x1, sg$y1, col="gray")
  
  graphics::text(0:maxA,0,0:maxA,cex=.6,col="lightgray")
  graphics::text(maxA,1:maxB,(maxA+1:maxB),cex=.6,col="lightgray")
  
  for (i in 1:length(tt)) {
    stp=ddply(melt(tt[[i]]$tbl),~I(Var1+Var2),function(x)x[which.max(x$value),]-1)
    graphics::lines(stp$Var1, stp$Var2, col=i+1, lw=2)
  }
  lgd = sapply(tt, function(x)x$m$ss$group)
  graphics::legend("topleft", legend=as.character(lgd), lty=1,col=1+(1:length(lgd)),cex=.7,box.lty=0)
  graphics::box()
}




#' Simple Test-Item Analysis
#'
#' Produce tables with simple test-item analysis statistics
#'
#'
#' @param db A handle to the database, i.e. the output of \code{start_new_project}
#' or \code{open_project}
#' @return A list of two data frames, itemStats and testStats.
#' @author Ivailo Partchev
#' @details When an item is contained in more than one booklet, the 
#' item statistics are averaged across booklets with the number of 
#' persons in each booklets serving as weights.
#' @export
#'
tia_tables = function(db) {
  make_ReST(db)
  sco = dbGetQuery(db, "select * from ReST")
  sco = ddply(sco, ~booklet+person, function(x){x$sumScore=sum(x$score);x})
  iStats = ddply(sco, ~booklet+item, function(x){
    data.frame(meanScore=mean(x$score),
               sdScore=sd(x$score),
               rit=cor(x$score,x$sumScore),
               rir=cor(x$score,(x$sumScore-x$score)),
               n=length(x$score))
  })
  mxSc = dbGetQuery(db, "select item, max(score) as maxScore from rules group by item")
  iStats = merge(iStats, mxSc, all.x=T, all.y=F)
  iStats$pvalue = iStats$meanScore / iStats$maxScore
  tStats = ddply(iStats, "booklet", testlev)
  sTia = ddply(iStats, "item", simpleTIA)
  list(itemStats=sTia, testStats=tStats)
}


#' Verbal aggression data
#' 
#' A data set of self-reported verbal behaviour in different frustrating
#' situations (Vansteelandt, 2000)
#' 
#' 
#' @name verbAggrData
#' @docType data
#' @format A data set with 316 rows and 25 columns.
#' @keywords datasets
NULL

#' Scoring rules for the verbal aggression data
#' 
#' A set of (trivial) scoring rules for the verbal 
#' aggression data set
#' 
#' 
#' @name verbAggrRules
#' @docType data
#' @format A data set with 72 rows and 3 columns (item, response, score).
#' @keywords datasets
NULL

#' Item properties in the verbal aggression data
#' 
#' A data set of item properties related to the verbal
#' aggression data
#' 
#' 
#' @name verbAggrProperties
#' @docType data
#' @format A data set with 24 rows and 5 columns.
#' @keywords datasets
NULL

#' Item properties in the PISA 2012 example
#' 
#' A data set of item properties in the PISA 2012 example (see the
#' help screen for function add_booklet)
#' 
#' 
#' @name PISA_item_class
#' @docType data
#' @format A data set with 109 rows and 6 columns.
#' @keywords datasets
NULL
