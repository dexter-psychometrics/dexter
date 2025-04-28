

library(dplyr)


verbAggCopy = function(pth = test_path('testdata/verbAggression.db'))
{
  con = DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  db = open_project(pth)
  
  RSQLite::sqliteCopyDatabase(db, con)
  
  close_project(db)
  return(con)
}

df_join_equal = function(..., join_by, tol_double=1e-10)
{
  l=list(...)
  if(!all(sapply(l, inherits, what='data.frame'))){
    message('not a data.frame')
    return(FALSE)
  } 
  if(length(l)<2)
  {
    message('neeed more than one df to compare')
    return(FALSE)
  }
  a = l[[1]]
  cn = setdiff(colnames(a), join_by)
  for(b in l[2:length(l)])
  {
    if(!all(dim(a)==dim(b))) return(FALSE)
    if(!length(intersect(colnames(a),colnames(b))) == ncol(a)) return(FALSE)
    
    tst = inner_join(a,b,by=join_by,relationship="one-to-one")
    if(!nrow(a)==nrow(b) || nrow(a) != nrow(tst))
    {
      message('no matching join')
      return(FALSE)
    }
    for(column in cn)
    {
      if(class(a[[column]]) != class(b[[column]])){
        message('type mismatch', column)
        return(FALSE)
      } 
      if(is.double(a[[column]]))
      {
        if(max(abs(tst[[sprintf('%s.x',column)]] - tst[[sprintf('%s.y',column)]])) > tol_double)
        {
          message('double mismatch',column)
          return(FALSE)
        }
        
      } else if(!all(tst[[sprintf('%s.x',column)]] == tst[[sprintf('%s.y',column)]])){
        message('value mismatch',column)
        return(FALSE)
      } 
    }
  }
  TRUE
}