#' @title Computes a list of GridTopology objects
#' 
#' @description This function computes a list of GridTopology 
#' objects with the cells to be processed in each thread.  
#' 
#' @param cellcentre.offset see ?GridTopology.
#' @param cellsize see ?GridTopology.
#' @param cells.dim see ?GridTopology.
#' @param thread.size integer, vector with number of cells in each thread.
#' @param ncores An integer with the number of cores. Default is 1.
#' 
#' @docType methods
#' @export
computeGridTopologyList = function(cellcentre.offset, cellsize,
                                   cells.dim, thread.size, ncores=1)
{
  if(missing(thread.size))
    thread.size = c(trunc(cells.dim[1]/ncores), trunc(cells.dim[2]/ncores))
  
  if(length(thread.size)==1)
    thread.size = rep(thread.size, 2)
  
  threadsSizeX = computeThreadsSize(cells.dim[1], thread.size[1])
  threadsSizeY = computeThreadsSize(cells.dim[2], thread.size[2])
  
  res = list()
  k=1
  y = cellcentre.offset[2]
  for(I in threadsSizeY){
    x = cellcentre.offset[1]
    for(J in threadsSizeX){
      res[[k]] = GridTopology(c(x,y), c(cellsize[1],cellsize[2]), c(J, I))
      x = x + J*cellsize[1]
      k = k + 1
    }
    y = y + I*cellsize[2]
  }
  
  res
}

#' @title longitude and latitude to MODIS column and row. 
#' 
#' @description This function computes the MODIS grid 
#' attributes from a pair of coordinates.
#' 
#' @param x An SpatialPoints object or the longitude.
#' @param pixelsize A number. The MODIS resolution.
#' @param projCRS A CRS object.
#' 
#' @docType methods
#' @export
longLatTomodisColRow = function(x, y=NULL, pixelsize, projCRS=NULL){
  require(sp)
  projCRSSinu = CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
  tilesize = 1111950.51966666709631681442
  nhtiles = 36
  nvtiles = 18
  ncells = trunc(tilesize/pixelsize)
  xmin = - tilesize*nhtiles/2
  ymax = tilesize*nvtiles/2
  
  if( class(x)!="SpatialPoints" & is.null(projCRS))
    stop("Missin either an SpatialPoints object or the CRS projection")
  
  sppoint = x
  if( class(sppoint)!="SpatialPoints" )
    sppoint = SpatialPoints(coords = data.frame(longitude=x, latitude=y), projCRS)
  
  longitude = coordinates(sppoint)[,"longitude"]
  latitude = coordinates(sppoint)[,"latitude"]
  sppoint = spTransform(sppoint, projCRSSinu)
  x = coordinates(sppoint)[,"longitude"]
  y = coordinates(sppoint)[,"latitude"]
  
  dx = x - xmin
  dy = ymax - y 
  
  col = trunc(dx/pixelsize)
  row = trunc(dy/pixelsize)
  
  h = abs(trunc(dx / tilesize))
  v = abs(trunc(dy / tilesize))
  
  j = col - ncells * h
  i = row - ncells * v
  
  res = data.frame(longitude, latitude, col, row, h, v, j, i)
  rownames(res) = NULL
  res
}




#' @title MODIS column and row to longitude and latitude. 
#' 
#' @description This function computes the longitude and latitude for a
#' MODIS column and row index.
#' 
#' @param col An integer. A global index of the MODIS column
#' or a local index of the column within a specific h:v tile
#' @param row An integer. A global index of the MODIS row
#' or a local index of the row within a specific h:v tile
#' @param h An integer. Horizontal MODIS tile 
#' @param v An integer. Vertical MODIS tile 
#' @param pixelsize A number. The MODIS resolution.
#' @param projCRS A CRS object.
#' 
#' @docType methods
#' @export
modisColRowToLongLat = function(col, row, h=NULL, v=NULL, pixelsize, projCRS=NULL){
  require(sp)
  projCRSSinu = CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
  if(class(projCRS)!="CRS")
    projCRS = CRS("+proj=longlat +datum=WGS84 +no_defs")
  
  tilesize = 1111950.51966666709631681442 # meters
  ncells = trunc(tilesize/pixelsize)
  nhtiles = 36
  nvtiles = 18
  
  nhcells = ncells*nhtiles/2
  nvcells = ncells*nvtiles/2
  j = col - nhcells
  i = nvcells - row
  if( !is.null(h) )
    j = ncells*(h - nhtiles/2) + col
  if( !is.null(v) )
    i = ncells*(nvtiles/2 - v) - row 
  
  longitude = pixelsize*j + pixelsize/2
  latitude = pixelsize*i - pixelsize/2
  sppoint = SpatialPoints(coords = data.frame(longitude, latitude), projCRSSinu)
  sppoint = spTransform(sppoint, projCRS)
  res = data.frame(coordinates(sppoint))
  res
}

#' @title Recover MODIS dates
#' 
#' @description The function recovers the MODIS image dates. 
#' 
#' @param years An vector of integers with the years.Default 
#' is form 2000 to the system time year format(Sys.time(), \"\%Y\").
#' @param frequency An integer with the frequency in days.
#' @docType methods
#' @export
recoverMODISDates = function(years, frequency=16){
  if(missing(years))
    years = 2000:format(Sys.time(), "%Y")
  
  dates = unlist(lapply(years, function(y){  
    days = seq(from = 0, to = 365, by = frequency)
    as.Date(days, origin = as.Date(paste(y,"-01-01",sep="")))
  }))
  dates = as.Date(dates, origin="1970-01-01")
  return(dates)
}


#' @title Recover dates from year and julian number
#' 
#' @description The function recovers the date from an year and 
#' an julian  number. 
#' 
#' @param year An integers with the year.
#' @param n An integer with the julian number.
#' @docType methods
#' @export
recoverDatesFromJulian = function(year, n){
  dates = seq(as.Date(paste0(year,"-01-01")), as.Date(paste0(year,"-12-31")), by=1)
  N = as.numeric(format(dates, "%j"))
  i = which(n==N)
  return(dates[i])
}


#' @title Recover MODIS time index from a date
#' 
#' @description The function recovers the MODIS time index 
#' from a date
#' 
#' @param date 
#' @docType methods 
#' @export
recoverMODISTimeIndex = function(d){
  years = 2000:format(Sys.time(), "%Y")
  dates = unlist(lapply(years, function(y){  
    days = seq(from = 0, to = 365, by = 16)
    as.Date(days, origin = as.Date(paste(y,"-01-01",sep="")))
  }))
  dates = as.Date(dates, origin="1970-01-01")
  i = which.min(abs(dates-d))
  return(i)
}
