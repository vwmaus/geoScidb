
#' @title Raster from scidb array. 
#' 
#' @description This function creates and project MODIS grid from slices of 
#' a scidb array.
#' 
#' @param arrayname A carachter. SciDB array name.
#' @param slices A numeric vector. The array index for each desired time slices.
#' @param mincol An integer. The global index for the minimum MODIS column.
#' @param maxcol An integer. The global index for the maximum MODIS column.
#' @param minrow An integer. The global index for the minimum MODIS row.
#' @param maxrow An integer. The global index for the maximum MODIS row.
#' @param pixelsize A number. The MODIS resolution in meters.
#' @param host A carachter.
#' @param port An integer. 
#' @param user A carachter.
#' @param password A carachter.
#' @param filepath A carachter.
#' @param overwrite Logical. If TRUE, raster file will be overwritten if it exists. Default is FALSE.
#' 
#' @docType methods
#' @export
modisrasterFromSciDBArray = function(arrayname, slices, mincol, maxcol, minrow, maxrow, pixelsize, 
                                     hostname, port, user, password, filepath, overwrite=FALSE)
{
  #   arrayname="CLASSIFICATION_MT_2000_2014"
  #   slices=2010 
  #   mincol=59221
  #   maxcol=59300
  #   minrow=48501
  #   maxrow=48520
  #   pixelsize=231.656358263889
  #   hostname="https://gis-bigdata.uni-muenster.de"
  #   port=48922
  #   user="maus"
  #   password="96V14w57m24" 
  #   filepath="/home/maus/workspace/results/small_test"
  #   overwrite=TRUE
  projCRSSinu = CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
  hoststring = unlist(strsplit(hostname, split="://"))
  hostscidb = hoststring[length(hoststring)]
  hostname = paste0("https://", hostscidb)
  scidbconnect(host=hostscidb, port=port, username=user, password=password)
  N = length(slices)
  
  # Remove existing temp array from scidb
  if(any(scidblist()=="tmp_raster1"))
    iquery("remove(tmp_raster1)", return=FALSE)
  if(any(scidblist()=="tmp_raster2"))
    iquery("remove(tmp_raster2)", return=FALSE)
  
  cat("\nInsert array slice 1/",N,sep="")
  iquery(paste0("store(between(slice(",arrayname,", year_id,",slices[N],"),",mincol,",",minrow,",",maxcol,",",maxrow,"),tmp_raster1)"), return = FALSE) 
  if(length(slices)>1)
    lapply((N-1):1, function(i){
      cat("\nInsert array slice ",N-i+1,"/",N,sep="")
      iquery(paste0("store(join(",paste0("between(slice(",arrayname,", year_id, ",slices[i],"),",mincol,",",minrow,",",maxcol,",",maxrow,")", collapse = ","), ",tmp_raster1),tmp_raster2)"), return=FALSE)
      iquery("remove(tmp_raster1)", return=FALSE)
      iquery("rename(tmp_raster2,tmp_raster1)", return=FALSE)
    })
  # Create raster from scidb array
  aux = unlist(strsplit(schema(scidb("tmp_raster1")), split="\\["))
  arraydim = paste(dimensions(scidb("tmp_raster1")), unlist(lapply(unlist(strsplit(aux[2], split="=")), function(i) unlist(strsplit(i, split=","))[1]))[-1], sep="=")
  TMPSCHEMA = paste0(aux[1],paste0("[row_id=",minrow,":",maxrow,",1024,0,col_id=",mincol,":",maxcol,",1024,0]"))
  iquery(paste0("store(redimension(tmp_raster1,",TMPSCHEMA,"),tmp_raster2)"))
  tmpfile = paste(c(unlist(strsplit(filepath, split = "/")), "raster.tmp"), collapse = "/")
  systemstring = paste("gdal_translate -of GTiff \"SCIDB:array=tmp_raster2 host=",hostname," port=",port," user=",user," password=",password,"\" ",tmpfile,sep="")
  system(systemstring)
  
  # Remove temp array from scidb
  if(any(scidblist()=="tmp_raster1"))
    iquery("remove(tmp_raster1)", return=FALSE)
  if(any(scidblist()=="tmp_raster2"))
    iquery("remove(tmp_raster2)", return=FALSE)
  
  # Project raster
  ceter_min_coords = modisColRowToLongLat(mincol, maxrow, h=NULL, v=NULL, pixelsize, projCRS=projCRSSinu)
  min_coords = ceter_min_coords - pixelsize/2
  ceter_max_coords = modisColRowToLongLat(maxcol, minrow, h=NULL, v=NULL, pixelsize, projCRS=projCRSSinu)
  max_coords = ceter_max_coords + pixelsize/2
  r_extent = extent(t(rbind(min_coords, max_coords)))
  r2=brick(tmpfile, xmn=min_coords$longitude, xmx=max_coords$longitude, ymn=min_coords$latitude, ymx=max_coords$latitude, crs=projCRSSinu)
  r2@extent = extent(t(rbind(min_coords, max_coords)))
  names(r2) = paste(arrayname,slices, sep=".")
  filename = paste(c(unlist(strsplit(filepath, split = "/")), paste("slice_",slices[1],"_",slices[length(slices)],"_array_",arrayname,".tif",sep="") ), collapse = "/")
  res = writeRaster(x=r2, filename = filename, format = "GTiff", overwrite = overwrite)
  system(paste("rm -rf ",tmpfile))
  res
}




#' @title Raster from scidb array. 
#' 
#' @description This function creates and project MODIS grid from slices of 
#' a scidb array.
#' 
#' @param arrayname A carachter. SciDB array name.
#' @param slices A numeric vector. The array index for each desired time slices.
#' @param mincol An integer. The global index for the minimum MODIS column.
#' @param maxcol An integer. The global index for the maximum MODIS column.
#' @param minrow An integer. The global index for the minimum MODIS row.
#' @param maxrow An integer. The global index for the maximum MODIS row.
#' @param pixelsize A number. The MODIS resolution in meters.
#' @param host A carachter.
#' @param port An integer. 
#' @param user A carachter.
#' @param password A carachter.
#' @param filepath A carachter.
#' @param RGBcode an R data.frame with the class number and R,G,B values in the columns.
#' @param overwrite Logical. If TRUE, raster file will be overwritten if it exists. Default is FALSE.
#' 
#' @docType methods
#' @export
modisrasterFromSciDBArrayRGB = function(arrayname, slices, mincol, maxcol, minrow, maxrow, pixelsize, 
                                        hostname, port, user, password, filepath, RGBcode, overwrite=FALSE)
{
  projCRSSinu = CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
  hoststring = unlist(strsplit(hostname, split="://"))
  hostscidb = hoststring[length(hoststring)]
  hostname = paste0("https://", hostscidb)
  scidbconnect(host=hostscidb, port=port, username=user, password=password)
  N = length(slices)
  
  # Remove existing temp array from scidb
  if(any(scidblist()=="tmp_raster1"))
    iquery("remove(tmp_raster1)", return=FALSE)
  if(any(scidblist()=="tmp_raster2"))
    iquery("remove(tmp_raster2)", return=FALSE)
  
  
  # Get RGB code for each class
  I = as.numeric(RGBcode$class_id)
  R = RGBcode$R
  G = RGBcode$G
  B = RGBcode$B
  strR = paste0("R,uint8(",paste0(lapply(seq_along(I), function(i){
    paste0("iif(class_id=uint8(",I[i],"), uint8(",R[i],")")
  }), collapse = ","), ",200", paste0(rep(")", length(I)+1), collapse = ""))
  strG = paste0("G,uint8(",paste0(lapply(seq_along(I), function(i){
    paste0("iif(class_id=uint8(",I[i],"), uint8(",G[i],")")
  }), collapse = ","), ",200", paste0(rep(")", length(I)+1), collapse = ""))
  strB = paste0("B,uint8(",paste0(lapply(seq_along(I), function(i){
    paste0("iif(class_id=uint8(",I[i],"), uint8(",B[i],")")
  }), collapse = ","), ",200", paste0(rep(")", length(I)+1), collapse = ""))
  
  res = lapply(1:N, function(i){
    cat("\nInsert array slice ",i,"/",N,sep="")
    iquery(
      paste(
        "store(
        project(
        apply(
        between(slice(",arrayname,", year_id,",slices[i],"), ",mincol,",",minrow,",",maxcol,",",maxrow,"),
        ",strR,",",strG,",",strB,"
        ),
        R, G, B
        ), tmp_raster1)"
            )
      , return = FALSE) 
    
    # Create raster from scidb array
    aux = unlist(strsplit(schema(scidb("tmp_raster1")), split="\\["))
    arraydim = paste(dimensions(scidb("tmp_raster1")), unlist(lapply(unlist(strsplit(aux[2], split="=")), function(i) unlist(strsplit(i, split=","))[1]))[-1], sep="=")
    TMPSCHEMA = paste0(aux[1],paste0("[row_id=",minrow,":",maxrow,",1024,0,col_id=",mincol,":",maxcol,",1024,0]"))
    iquery(paste0("store(redimension(tmp_raster1,",TMPSCHEMA,"),tmp_raster2)"))
    tmpfile = paste(c(unlist(strsplit(filepath, split = "/")), ".raster.tmp"), collapse = "/")
    systemstring = paste("gdal_translate -of GTiff \"SCIDB:array=tmp_raster2 host=",hostname," port=",port," user=",user," password=",password,"\" ",tmpfile,sep="")
    system(systemstring)
    
    # Remove temp array from scidb
    if(any(scidblist()=="tmp_raster1"))
      iquery("remove(tmp_raster1)", return=FALSE)
    if(any(scidblist()=="tmp_raster2"))
      iquery("remove(tmp_raster2)", return=FALSE)
    
    # Project raster
    ceter_min_coords = modisColRowToLongLat(mincol, maxrow, h=NULL, v=NULL, pixelsize, projCRS=projCRSSinu)
    min_coords = ceter_min_coords - pixelsize/2
    ceter_max_coords = modisColRowToLongLat(maxcol, minrow, h=NULL, v=NULL, pixelsize, projCRS=projCRSSinu)
    max_coords = ceter_max_coords + pixelsize/2
    r_extent = extent(t(rbind(min_coords, max_coords)))
    r2=brick(tmpfile, xmn=min_coords$longitude, xmx=max_coords$longitude, ymn=min_coords$latitude, ymx=max_coords$latitude, crs=projCRSSinu)
    r2@extent = extent(t(rbind(min_coords, max_coords)))
    names(r2) = paste(arrayname,c("R","G","B"), sep=".")
    filename = paste(c(unlist(strsplit(filepath, split = "/")), paste("RGB_slice_",slices[i],"_array_",arrayname,".tif",sep="") ), collapse = "/")
    res = writeRaster(x=r2, filename = filename, format = "GTiff", overwrite = overwrite)
    system(paste("rm -rf ",tmpfile))
    res
  })
  res
}


# modisrasterFromSciDBArrayDTWValue = function(arrayname, slices, mincol, maxcol, minrow, maxrow, pixelsize, 
#                                      hostname, port, user, password, filepath, overwrite=FALSE)
# {
#   projCRSSinu = CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
#   TMPARRAY = "TMP_RASTER_ARRAY"
#   hoststring = unlist(strsplit(hostname, split="://"))
#   hostscidb = hoststring[length(hoststring)]
#   hostname = paste0("https://", hostscidb)
#   scidbconnect(host=hostscidb, port=port, username=user, password=password)
#   
#   # Remove existing temp array from scidb
#   if(any(scidblist()==TMPARRAY))
#     scidbremove(TMPARRAY, force=TRUE)
#   
#   patternnames = scidb_attributes(scidb(arrayname))
#   
#   out = lapply(slices, function(year){
#     cat("\nSaving slice ",which(slices==year),"/",length(slices))
#      TMPSCHEMA = paste("<",paste(paste(patternnames, collapse = ":double,"), ":double", sep=""),">",
#                        " [row_id=",minrow,":",maxrow,",256,0,col_id=",mincol,":",maxcol,",256,0]", sep="")
#      
#     iquery(paste("store(
#                  redimension(
#                  slice(",arrayname,", year_id, ",year,"),
#                  ",TMPSCHEMA,"
#                  ),
#                  ",TMPARRAY,"
#     )"), return=FALSE)
# 
#     # Create raster from scidb array
#     tmpfile = paste(c(unlist(strsplit(filepath, split = "/")), "raster.tmp"), collapse = "/")
#     systemstring = paste("gdal_translate -of GTiff \"SCIDB:array=",TMPARRAY," host=",hostname," port=",port," user=",user," password=",password,"\" ",tmpfile,sep="")
#     system(systemstring)
#     
#     # Remove temp array from scidb
#     if(any(scidblist()==TMPARRAY))
#       scidbremove(TMPARRAY, force=TRUE)
#     
#     # Project raster
#     ceter_min_coords = modisColRowToLongLat(mincol, maxrow, h=NULL, v=NULL, pixelsize, projCRS=projCRSSinu)
#     min_coords = ceter_min_coords - pixelsize/2
#     ceter_max_coords = modisColRowToLongLat(maxcol, minrow, h=NULL, v=NULL, pixelsize, projCRS=projCRSSinu)
#     max_coords = ceter_max_coords + pixelsize/2
#     r_extent = extent(t(rbind(min_coords, max_coords)))
#     r2=brick(tmpfile, xmn=min_coords$longitude, xmx=max_coords$longitude, ymn=min_coords$latitude, ymx=max_coords$latitude, crs=projCRSSinu)
#     r2@extent = extent(t(rbind(min_coords, max_coords)))
#     names(r2) = paste("slice_",slices,"_array_",arrayname,".",patternnames, sep="")
#     filename = paste(c(unlist(strsplit(filepath, split = "/")), paste("slice_",year,"_array_",arrayname,".tif",sep="") ), collapse = "/")
#     res = writeRaster(x=r2, filename = filename, format = "GTiff", overwrite = overwrite)
#     system(paste("rm -rf ",tmpfile))
#     res 
# })
# }

