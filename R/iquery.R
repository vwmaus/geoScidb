#' @title Build SciDB query
#' 
#' @description The function builds the SciDB query for time 
#' series analysis. 
#' @param INPUTARRAY
#' @param OUTPUTSCHEMA
#' @param PATTERNNAMES
#' @param XMIN
#' @param XMAX
#' @param YMIN
#' @param YMAX
#' @param TMIN
#' @param TMAX
#' @param COLIDS
#' @param ROWIDS
#' @param CHUNKCOL
#' @param CHUNKROW
#' @param THRESHOLD
#' @param FILL
#' @param NORMALIZE
#' @param REALDAY
#' @param METHOD
#' @param THETA
#' @param ALPHA
#' @param BETA
#' @param DELAY
#' @param NCORES
#' @param LIBRARYPATH
#' @param EXPPRPATH
#' @docType methods
#' @export
buildSciDBDTWQuery = function(INPUTARRAY, OUTPUTSCHEMA, PATTERNNAMES, 
                              XMIN, XMAX,
                              YMIN, YMAX,
                              TMIN, TMAX,
                              COLIDS, ROWIDS,
                              CHUNKCOL, CHUNKROW,
                              THRESHOLD, FILL, NORMALIZE,
                              REALDAY, METHOD,
                              THETA, ALPHA,
                              BETA, DELAY,
                              NCORES,
                              LIBRARYPATH, EXPPRPATH){
  
  rOut = 0:(length(PATTERNNAMES)+2)
  attrRename = gsub(" ", ",",paste(paste("expr_value", rOut, sep="_"), c("col", "row", "year", PATTERNNAMES), collapse = ","))
  
  attrRedimension = paste(c("col_id", "row_id", "year_id", PATTERNNAMES), collapse = ",")
  
  res = paste(" redimension(
              project(
              apply(
              attribute_rename(
              r_exec(
              project(
              apply(
              redimension(
              between(",INPUTARRAY,",",XMIN,",",YMIN,",",TMIN,",",XMAX,",",YMAX,",",TMAX,"),
              <evi:int16>[col_id=",COLIDS,",",CHUNKCOL,",0,row_id=",ROWIDS,",",CHUNKROW,",0,time_id=",TMIN,":",TMAX,",",2*TMAX,",0]
              ),
              devi, double(evi), dcol, double(col_id), drow, double(row_id), dtime, double(time_id)
              ), 
              devi, dcol, drow, dtime
              ),
              'output_attrs=",length(rOut),"',
              'expr=
              output_attrs=",length(rOut),"
              LIBRARYPATH=\"",LIBRARYPATH,"\"
              FROM=\"",FROM,"\"
              TO=\"",TO,"\"
              THRESHOLD=",THRESHOLD,"
              FILL=",FILL,"
              TMIN=",TMIN,"
              TMAX=",TMAX,"
              NORMALIZE=",NORMALIZE,"
              REALDAY=",REALDAY,"
              METHOD=\"",METHOD,"\"
              THETA=",THETA,"
              ALPHA=",ALPHA,"
              BETA=",BETA,"
              DELAY=",DELAY,"
              NCORES=",NCORES,"
              source(\"",EXPPRPATH,"\")
              res'
              ),
              ",attrRename,"
              ),
              col_id,  int64(col),
              row_id,  int64(row),
              year_id, int64(year)
              ),",
                              attrRedimension,"
              ),",
                        OUTPUTSCHEMA,"
              )", sep="")
return(res)

}







#' @title Build SciDB query
#' 
#' @description The function builds the SciDB query for time 
#' series analysis. 
#' @param INPUTARRAY
#' @param INPUTSCHEMA
#' @param OUTPUTSCHEMA
#' @param PATTERNNAMES 
#' @param XMIN
#' @param XMAX
#' @param YMIN
#' @param YMAX
#' @param TMIN
#' @param TMAX
#' @param FILL
#' @param NCORES
#' @param EXPPRPATH
#' @docType methods
#' @export
buildSciDBPostProcQuery = function(INPUTARRAY, 
                                   INPUTSCHEMA,
                                   OUTPUTSCHEMA,
                                   PATTERNNAMES, 
                                   XMIN, XMAX,
                                   YMIN, YMAX,
                                   TMIN, TMAX,
                                   FILL, NCORES, 
                                   EXPPRPATH){
  
  attrRename = gsub(" ", ",",paste(paste("expr_value", 0:3, sep="_"), c("col", "row", "year", "class"), collapse = ","))
  
  res = paste(" 
              redimension(
              apply(
              attribute_rename(
              r_exec(
              redimension(
              apply(
              between(",INPUTARRAY,",",XMIN,",",YMIN,",",TMIN,",",XMAX,",",YMAX,",",TMAX,"),
              dcol, double(col_id), drow, double(row_id), dyear, double(year_id)
              ), 
              ",INPUTSCHEMA,"
              ),
              'output_attrs=",4,"',
              'expr=
              output_attrs=",4,"
              FILL=",FILL,"
              NCORES=",NCORES,"
              PATTERNNAMES=c(\"",paste(PATTERNNAMES, collapse ="\",\"", sep=""), "\")
              source(\"",EXPPRPATH,"\")
              res'
              ),
              ",attrRename,"
              ),
              col_id,  int64(col),
              row_id,  int64(row),
              year_id, int64(year),
              class_id, uint8(class)
              ),",OUTPUTSCHEMA,"
              )", sep="")
return(res)
}


#' @title Build SciDB query
#' 
#' @description The function builds the SciDB query for time 
#' series analysis. 
#' @param INPUTARRAY
#' @param INPUTSCHEMA
#' @param OUTPUTSCHEMA
#' @param PATTERNNAMES 
#' @param XMIN
#' @param XMAX
#' @param YMIN
#' @param YMAX
#' @param TMIN
#' @param TMAX
#' @param FILL
#' @param NCORES
#' @param EXPPRPATH
#' @docType methods
#' @export
buildSciDBPostProcSecVegetationQuery = function(INPUTARRAY, 
                                                INPUTSCHEMA,
                                                OUTPUTSCHEMA,
                                                PATTERNNAMES, 
                                                XMIN, XMAX,
                                                YMIN, YMAX,
                                                TMIN, TMAX,
                                                FILL, NCORES, 
                                                EXPPRPATH){
  
  attrRename = gsub(" ", ",",paste(paste("expr_value", 0:3, sep="_"), c("col", "row", "year", "class"), collapse = ","))
  
  res = paste(" 
              redimension(
              apply(
              attribute_rename(
              r_exec(
              redimension(
              apply(
              between(",INPUTARRAY,",",XMIN,",",YMIN,",",TMIN,",",XMAX,",",YMAX,",",TMAX,"),
              dcol, double(col_id), drow, double(row_id), dyear, double(year_id), dclass, double(class_id) 
              ), 
              ",INPUTSCHEMA,"
              ),
              'output_attrs=",4,"',
              'expr=
              output_attrs=",4,"
              FILL=",FILL,"
              NCORES=",NCORES,"
              PATTERNNAMES=c(\"",paste(PATTERNNAMES, collapse ="\",\"", sep=""), "\")
              source(\"",EXPPRPATH,"\")
              res'
              ),
              ",attrRename,"
              ),
              col_id,  int64(col),
              row_id,  int64(row),
              year_id, int64(year),
              class_id, uint8(class)
              ),",OUTPUTSCHEMA,"
              )", sep="")
return(res)
}


#' @title Build SciDB 2 query (Faster)
#' 
#' @description The function builds the SciDB query for time 
#' series analysis. 
#' @param INPUTARRAY
#' @param OUTPUTSCHEMA
#' @param PATTERNNAMES
#' @param XMIN
#' @param XMAX
#' @param YMIN
#' @param YMAX
#' @param TMIN
#' @param TMAX
#' @param COLIDS
#' @param ROWIDS
#' @param CHUNKCOL
#' @param CHUNKROW
#' @param THRESHOLD
#' @param FILL
#' @param NORMALIZE
#' @param REALDAY
#' @param METHOD
#' @param THETA
#' @param ALPHA
#' @param BETA
#' @param DELAY
#' @param NCORES
#' @param LIBRARYPATH
#' @param EXPPRPATH
#' @docType methods
#' @export
buildSciDBDTWQuery2 = function(INPUTARRAY, OUTPUTSCHEMA, PATTERNNAMES, 
                               XMIN, XMAX,
                               YMIN, YMAX,
                               TMIN, TMAX,
                               COLIDS, ROWIDS,
                               CHUNKCOL, CHUNKROW,
                               THRESHOLD, FILL, NORMALIZE,
                               REALDAY, METHOD,
                               THETA, ALPHA,
                               BETA, DELAY,
                               NCORES,
                               LIBRARYPATH, EXPPRPATH){
  
  rOut = 0:(length(PATTERNNAMES)+2)
  attrRename = gsub(" ", ",",paste(paste("expr_value", rOut, sep="_"), c("col", "row", "year", PATTERNNAMES), collapse = ","))
  
  attrRedimension = paste(c("col_id", "row_id", "year_id", PATTERNNAMES), collapse = ",")
  
  res = paste("redimension(
              project(
              apply(
              attribute_rename(
              r_exec(
              project(
              apply(
              redimension(
              subarray(
              apply(
              between(",INPUTARRAY,",",XMIN,",",YMIN,",",TMIN,",",XMAX,",",YMAX,",",TMAX,"), 
              icol, int64(col_id), irow, int64(row_id), itime, int64(time_id)
              ), 
              ",XMIN,",",YMIN,",",TMIN,",",XMAX,",",YMAX,",",TMAX,"
              ),
              <evi:int16>[icol=",COLIDS,",",CHUNKCOL,",0,irow=",ROWIDS,",",CHUNKROW,",0,itime=",TMIN,":",TMAX,",",2*TMAX,",0]
              ),
              devi, double(evi), dcol, double(icol), drow, double(irow), dtime, double(itime)
              ),
              devi, dcol, drow, dtime
              ),
              'output_attrs=",length(rOut),"',
              'expr=
              output_attrs=",length(rOut),"
              LIBRARYPATH=\"",LIBRARYPATH,"\"
              FROM=\"",FROM,"\"
              TO=\"",TO,"\"
              THRESHOLD=",THRESHOLD,"
              FILL=",FILL,"
              TMIN=",TMIN,"
              TMAX=",TMAX,"
              NORMALIZE=",NORMALIZE,"
              REALDAY=",REALDAY,"
              METHOD=\"",METHOD,"\"
              THETA=",THETA,"
              ALPHA=",ALPHA,"
              BETA=",BETA,"
              DELAY=",DELAY,"
              NCORES=",NCORES,"
              source(\"",EXPPRPATH,"\")
              res'
              ),
              ",attrRename,"
              ),
              col_id,  int64(col),
              row_id,  int64(row),
              year_id, int64(year)
              ),",
              attrRedimension,"
              ),",
              OUTPUTSCHEMA,"
              )", sep="")
  return(res)
  
}
