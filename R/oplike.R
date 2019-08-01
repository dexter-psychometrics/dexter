


#' Start a new project from oplm files
#'
#' Creates a dexter project database and fills it with response data based on a .dat and .scr file
#'
#'
#' @param dbname filename/path of new dexter project database (will be overwritten if already exists)
#' @param scr_path path to the .scr file
#' @param dat_path path to the .dat file
#' @param booklet_position vector of start and end of booklet position in the dat file, e.g. c(1,4), 
#' all positions are counted from 1, start and end are both inclusive. If NULL, this is read from the scr file.
#' @param responses_start start position of responses in the .dat file. If NULL, this is read from the scr file.
#' @param response_length length of individual responses, default=1
#' @param person_id optionally, a vector of start and end position of person_id in the .dat file.
#' If NULL, person id's will be auto-generated.
#' @param missing_character vector of character(s) used to indicate missing in .dat file, 
#' default is to use both a space and a 9 as missing characters.
#' @param use_discrim if TRUE, the scores for the responses will be multiplied by the
#' discrimination parameters of the items
#' @param format not used, at the moment only the compressed format is supported.
#' @return a database connection object.
#' 
#' @details
#' start_new_project_from_oplm builds a complete dexter database from a .dat and .scr file in
#' the proprietary oplm format. Three custom variables are added to the database: 
#' booklet_on_off, item_local_on_off, item_global_on_off. These are taken from the .scr file
#' and can be used in predicates in the various dexter functions.
#' 
#' Booklet_position and responses_start are usually inferred from the scr file but since they
#' are sometimes misspecified in the scr file they can be overridden. Response_length is not 
#' inferred from the scr file since anything other than 1 is most often a mistake.
#' 
#' @examples
#'\donttest{
#' db = start_new_project_from_oplm('test.db',
#'    'path_to_scr_file', 'path_to_dat_file', 
#'    booklet_position=c(1,3), responses_start=101,
#'    person_id=c(50,62))
#'
#' prms = fit_enorm(db, 
#'    item_global_on_off==1 & item_local_on_off==1 & booklet_on_off==1)
#' 
#' 
#'}
start_new_project_from_oplm = function(dbname, scr_path, dat_path, 
                                       booklet_position = NULL, responses_start = NULL, response_length = 1, 
                                       person_id = NULL, missing_character = c(' ','9'), use_discrim=FALSE,
                                       format='compressed')
{
  if(format != 'compressed') stop(paste('Only compressed format is supported at this time', 
                                        'use the inexpand tool from oplm to compress your data'))
  
  scr = readSCR(scr_path)
  scr$itemLabels = trimws(scr$itemLabels)
  
  if(is.null(booklet_position))
    booklet_position = scr$booklet_position
  
  if(is.null(responses_start))
    responses_start = scr$responses_start
  
  if(response_length==1)
  { 
    rsplit = function(s) strsplit(s,"")
  } else if(response_length==2)
  {
    rsplit = function(s) lapply(strsplit(s, ""), function(x){paste0(x[c(TRUE, FALSE)], x[c(FALSE, TRUE)])})
  } else if(response_length==3)
  {
    rsplit = function(s) lapply(strsplit(s, ""), function(x){paste0(x[c(TRUE, FALSE, FALSE)], x[c(FALSE, TRUE, FALSE)], x[c(FALSE, FALSE, TRUE)])})
  } else 
  {
    stop('response length cannot be larger than 3.')
  }
  missing_character = trimws(missing_character)
  
  if(max(nchar(missing_character))>response_length) 
      stop('missing character may not be longer than the length of an individual response')
  
  if(!is.null(person_id))
  {
    if(length(person_id)!=2) stop('person_id should be a vector of length 2')
    if(person_id[2]<person_id[1]) stop('person_id is invalid, end position is smaller than start position')
  }
  if(length(booklet_position)!=2) stop('booklet_position should be a vector of length 2')
  if(booklet_position[2]<booklet_position[1]) stop('booklet_position is invalid, end position is smaller than start position')
  
  
  if(anyDuplicated(scr$itemLabels))
  {
    warning("screenfile contains duplicate itemlabels, labels have been adjusted to make them unique")
    scr$itemLabels = sprintf('%04i_%s',1:length(scr$itemLabels),scr$itemLabels)
  }
  
  rules = tibble(item_id = rep(scr$itemLabels, as.integer(scr$maxScore)+1+length(missing_character)), 
                 response = unlist(lapply(scr$maxScore, function(n) c(0:n,missing_character))),
                 item_score = unlist(lapply(scr$maxScore, function(n) c(0:n,rep(0,length(missing_character))))))
  
  if(use_discrim) rules$item_score = rules$item_score * rep(scr$discrim, as.integer(scr$maxScore)+1+length(missing_character))
  
  db = start_new_project(rules, dbname)
  
  finished = FALSE
  dbTransaction(db,{
    dbExecute(db,'ALTER TABLE dxBooklets ADD COLUMN booklet_on_off INTEGER NOT NULL DEFAULT 1;')
    dbExecute(db,'ALTER TABLE dxBooklet_design ADD COLUMN item_local_on_off INTEGER NOT NULL DEFAULT 1;')
    dbExecute(db,'ALTER TABLE dxItems ADD COLUMN item_global_on_off INTEGER NOT NULL DEFAULT 1;')
    
    dbExecute(db,'UPDATE dxItems SET item_global_on_off=:onoff WHERE item_id=:item;', 
                tibble(onoff=scr$globCal, item=scr$itemLabels))
  
    design = tibble(booklet_id = as.character(unlist(lapply(1:scr$nBook, function(x) rep(x,scr$nitBook[x])))),
                    item_nbr = unlist(scr$itemsBook),
                    item_position = unlist(lapply(scr$itemsBook, function(x) 1:length(x))),
                    onoff = unlist(scr$itemsOn)) %>%
              inner_join(tibble(item_nbr = 1:scr$nit, item_id = scr$itemLabels), by='item_nbr') %>%
              arrange(.data$booklet_id, .data$item_position) 
            
    dbExecute(db,'INSERT INTO dxBooklets(booklet_id, booklet_on_off) VALUES(:b,:onoff);',
                      tibble(b=as.character(1:scr$nBook),onoff=scr$bookOn))
    dbExecute(db,'INSERT INTO dxBooklet_design(booklet_id, item_id, item_position, item_local_on_off) 
                            VALUES(:booklet_id,:item_id,:item_position, :onoff);', 
					select(design, .data$booklet_id,.data$item_id,.data$item_position, .data$onoff))                
    
    vp = 1
    con = file(dat_path, "r")
    short_line = NULL
    while ( TRUE ) {
      lines = readLines(con, n = 5000, encoding='ascii')
      if ( length(lines) == 0 ) {
        break
      }
      cat('\rreading lines:', vp, '-', as.integer(vp + length(lines)-1))
      bkl = as.integer(substr(lines, booklet_position[1], booklet_position[2]))
      
      if(anyNA(bkl))
      {
        cat('\n')
        stop(paste0("empty booklet id's at position (", paste0(booklet_position, collapse=', ' ),")"))
      } else if(min(bkl) < 1 || max(bkl) > scr$nBook)
      {
        cat('\n')
        stop(paste0("The following booklet_id's in data are not present in the .scr file:\n",
                    paste(setdiff(bkl, 1:scr$nBook), collapse = ' ')))
      }
      
      
      rsp = rsplit(substring(lines, 
                        responses_start, 
                        responses_start + scr$nitBook[bkl] * response_length-1))
      
      if(any(sapply(rsp,length) < scr$nitBook[bkl]))
      {
        lnbr = which(sapply(rsp,length) < scr$nitBook[bkl])
        short_line = c(short_line, lnbr + vp - 1)
        sapply(scr$nitBook[bkl[lnbr]] - sapply(rsp[lnbr], length), rep, x = ' ')
        
        rsp[lnbr] = mapply(c, 
                           rsp[lnbr],  
                           lapply(scr$nitBook[bkl[lnbr]] - sapply(rsp[lnbr], length), rep, x = ' '),
                           SIMPLIFY = FALSE)
        
      }
      

      bkl=as.character(bkl)
      
      if(is.null(person_id))
      {
        prs = as.character(vp:(vp+length(lines)-1) )
      } else 
      {
        prs = trimws(substr(lines,person_id[1],person_id[2]))
      }
      vp = vp+length(lines)

      data = tibble(booklet_id = rep(bkl, sapply(rsp,length)), 
                    person_id = rep(prs, sapply(rsp,length)),
                    response=trimws(unlist(rsp)), 
                    item_position = unlist(sapply(rsp, function(n) 1:length(n),simplify=FALSE))) %>%
        inner_join(design, by=c('booklet_id','item_position')) %>%
        select(.data$person_id, .data$booklet_id, .data$item_id, .data$response)
      
      dbExecute(db,'INSERT INTO dxPersons(person_id) VALUES(:person_id);',tibble(person_id=prs) )
      dbExecute(db,'INSERT INTO dxAdministrations(person_id,booklet_id) VALUES(:person_id,:booklet_id);', 
                      tibble(person_id=prs,booklet_id=bkl))

      dbExecute(db,
                'INSERT INTO dxResponses(booklet_id,person_id,item_id,response) 
                                  VALUES(:booklet_id,:person_id,:item_id,:response);', 
                select(data, .data$booklet_id,.data$person_id,.data$item_id,.data$response))
    
    }
    close(con)
    finished=TRUE
    if(!is.null(short_line))
    {
      message('Too short lines in file:')
      print(short_line)
      warning('Some lines are too short for the number of responses that are supposed to be present ',
              'according to the design', call.=FALSE)
      
    }
    cat('\nvalidating constraints...\n')
 }, on_error = function(e)
    {
      if(!finished) close(con)
   
      if(finished & grepl('foreign key', e$message,ignore.case=TRUE))
      {

        # there are three deferred foreign key constraints
        # that are not caught until the end of the transaction (for efficiency reasons)
        # we can test them one by one
        unknown_responses = dbGetQuery(db,'SELECT item_id, response FROM dxResponses 
                                            EXCEPT SELECT item_id, response FROM dxScoring_rules;')
        if(nrow(unknown_responses) > 0)
        {
          cat('\nThe following responses were found in the data but they are not defined in the .scr file or coded as missing responses.')
          cat('Possible causes are that not all missing characters are correctly specified, your screen and dat files do not match ')
          cat('or responses_start is incorrect.\n')
          print(unknown_responses)
          e$message = 'Invalid responses, see output'
        } 
        # FOREIGN KEY (booklet_id, item_id) is highly unlikely to be violated. 
        # because this will most likely manifest in the previous constraint
        # FOREIGN KEY (booklet_id, person_id) can also not be wrong in this function because they are entered together

      }
      
      dbRollback(db)  
      dbDisconnect(db)
      stop(e, .call=FALSE)
    }, on_error_rollback=FALSE)
  return(db)
}



#' Read item parameters from oplm PAR or CML files
#'
#' @param par_path path to a file in the (binary) OPLM PAR format or the human readable CML format
#' @return depends on the input. For .PAR files a tibble with columns: item_id, item_score, beta, nbr,
#' for .CML files also several statistics columns that are outputted by OPLM as part of the calibration.
#' 
#' @details
#' It is occasionally useful to calibrate new items on an existing scale. This
#' function offers the possibility to read parameters from the proprietary oplm format 
#' so that they can be used to fix a new calibration in Dexter on an existing scale of items
#' that were calibrated in oplm.
#' 
#' @examples
#' \donttest{
#' par = read_oplm_par('/parameters.PAR')
#' f = fit_enorm(db, fixed_params=par)
#' }
read_oplm_par = function(par_path){
  if(grepl('\\.par$',par_path,perl=TRUE, ignore.case=TRUE))
  {
    return(
      lapply(readPAR(par_path), function(p)
      {
        tibble(item_id=trimws(p$label),
               item_score=c(1:p$ncat)*p$discr,
               beta=p$delta,
               nbr=p$ilabel)
      }) %>% 
        bind_rows() %>%
        arrange(.data$item_id, .data$item_score) %>%
        as.data.frame()
    )
  } else
  {
    return(as.data.frame(readCML(par_path)))
  }
}


### occasionally useful mainly for testing purposes to get pars from cml files  ###

readCML = function(cml_path)
{
  con = file(cml_path, "r")
  l = readLines(con, skipNul = TRUE)
  close(con)
  
  i = min(which(grepl('^\\s+RESULTS CALIBRATION\\.',l,perl=TRUE)))
  colnames = strsplit(trimws(l[i + 2]),'\\s+',perl=TRUE)[[1]]
  
  l = l[(i + 5):length(l)]
  l = l[1:(min(which(grepl('^-+$',l,perl=TRUE) == TRUE))-1)]
  
  pars = do.call(rbind, 
                 lapply(strsplit(trimws(l),'\\s+',perl=TRUE),
                        function(x)
                        {
                          if(length(x) == 0) return(NULL)
                          if(substr(x[1], 1, 1) == '[') x = c(NA, NA, gsub('\\D+','',x[1]), x[2:length(x)])
                          x[x=='---'] = NA
                          x
                        })
                 ) %>%
          as.data.frame(stringsAsFactors = FALSE) 
  
  pars[,c(1,3:ncol(pars))] = lapply(pars[,c(1,3:ncol(pars))],as.numeric)
  colnames(pars) = gsub('\\.+$','',make.names(tolower(colnames)))
  colnames(pars)[2:4] = c('item_id','item_score','beta') 
  
  # fill out the NA's in the first two columns
  pars %>% 
    fill(.data$nr, .data$item_id) %>% 
    rename(nbr = 'nr')
}


### functions borrowed from oplike package ###

readSCR = function (scrfile) 
{
  z = file(scrfile, "rb")
  n = readBin(z, 'int', size = 2, 3)
  nit = n[3]
  itemLabels = sapply(1:nit, function(x) {
    sl = readBin(z, 'int', size = 1, 1)
    rawToChar(readBin(z, 'raw', n = 8)[1:sl])
  })
  globCal = readBin(z, 'int', size = 1, nit)
  discrim = readBin(z, 'int', size = 1, nit)
  maxScore = readBin(z, 'int', size = 1, nit)
  parFixed = readBin(z, 'int', size = 1, nit)
  four = readBin(z, 'int', size = 1, 4)
  sl = readBin(z, 'int', size = 1, 1)
  jobname = rawToChar(readBin(z, 'raw', n = 12)[1:sl])
  five = readBin(z, 'int', size = 1, 5)
  sl = readBin(z, 'int', size = 1, 1)
  title = rawToChar(readBin(z, 'raw', n = 79)[1:sl])
  for (i in 1:20) {
    sl = readBin(z, 'int', size = 1, 1)
    if (sl > 0) 
      someComment = rawToChar(readBin(z, 'raw', n = sl))
  }
  sl = readBin(z, 'int', size = 1, 1)
  dataDir = rawToChar(readBin(z, 'raw', n = 60)[1:sl])
  sl = readBin(z, 'int', size = 1, 1)
  dataFile = rawToChar(readBin(z, 'raw', n = 12)[1:sl])
  expanded = readBin(z, 'int', size = 1, 1)
  expanded = 1 - expanded
  sl = readBin(z, 'int', size = 1, 1)
  fmt = rawToChar(readBin(z, 'raw', n = sl))
  nb = readBin(z, 'int', size = 2, 1)
  inUse = readBin(z, 'int', size = 1, nb)
  nMarg = readBin(z, 'int', size = 1, 1)
  nStat = readBin(z, 'int', size = 1, 1)
  margLabels = sapply(1:nMarg, function(x) {
    sl = readBin(z, 'int', size = 1, 1)
  readBin(z, 'raw', n = 8)[1:sl]
  })
  statLabels = sapply(1:nStat, function(x) {
    sl = readBin(z, 'int', size = 1, 1)
    readBin(z, 'raw', n = 8)[1:sl]
  })
  nitb = readBin(z, 'int', size = 2, nb)
  datum = readBin(z, 'int', size = 2, 5)
  margBook = statBook = rep(NA, nb)
  itemsBook = vector(mode = "list", length = nb)
  itemsOn = vector(mode = "list", length = nb)
  for (i in 1:nb) {
    rubbish = readBin(z, 'int', size = 2, 1)
    margBook[i] = readBin(z, 'int', size = 1, 1)
    statBook[i] = readBin(z, 'int', size = 1, 1)
    itemsBook[[i]] = readBin(z, 'int', size = 2, nitb[i])
    itemsOn[[i]] = readBin(z, 'int', size = 1, nitb[i])
  }
  close(z)

  fmt = as.integer(unlist(regmatches(fmt, gregexpr('\\d+',fmt,perl=TRUE))))

  list(nit = nit, nMarg = nMarg, nStat = nStat, itemLabels = itemLabels, 
       booklet_position = c(fmt[1], fmt[1] + fmt[2] - 1L),
       responses_start = fmt[3], response_length = fmt[5],
       margLabels = margLabels, statLabels = statLabels, 
       globCal = globCal, 
       discrim = discrim, maxScore = maxScore, parFixed = parFixed, 
       nBook = nb, bookOn = inUse, nitBook = nitb, 
       margBook = margBook, statBook = statBook, 
       itemsBook = itemsBook, itemsOn = itemsOn, 
       expanded = expanded)
}


readPAR = function(parfile)
{
  z = file(parfile, "rb")
  tit = rawToChar(readBin(z, 'raw', size = 1, 80))
  n = readBin(z, 'int', size = 2, 4)
  nit = n[3]
  nstat = -1 * n[1]
  statlab = sapply(1:nstat, function(x) rawToChar(readBin(z, 'raw', size = 1, 8)))
  pars = vector(mode = "list", length = nit)
  for (i in 1:nit) 
  {
    iLab = readBin(z, 'int', size = 2, 1)
    label = rawToChar(readBin(z, 'raw', size = 1, 8))
    disc = readBin(z, 'int', size = 2, 1)
    ncat = readBin(z, 'int', size = 2, 1)
    if (nstat > 1) 
    {
      idumm = readBin(z, 'int', size = 2, 1)
      dummy = readBin(z, 'int', size = 2, idumm)
    }
    delta = readBin(z, numeric(), size = 8, ncat)
    pars[[i]] = list(ncat = ncat, discr = disc, delta = delta, 
                    ilabel = iLab, label=label)
  }
  close(z)
  pars
}

