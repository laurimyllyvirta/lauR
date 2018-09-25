#' Proj4string for standard latitude-longitude projection
#'
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
llproj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

#' Constant to convert molec/cm2 to Dobson Units
#'
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
Dobson.unit = 2.69e16

#' Return true for all values of x that are NOT included in y
#'
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
'%notin%' <- function(x,y)!('%in%'(x,y))

#' Return a vector of the values of x that are included in y
#'
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
'%whichin%' <- function(x,y) x[x %in% y]

#' Return a vector of the values of x that are NOT included in y
#'
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
'%whichnotin%' <- function(x,y) x[x %notin% y]

#' Convert a data.frame to SpatialPointsDataFrame
#'
#' Convert a data.frame to SpatialPointsDataFrame by automatically identifying the columns that contain coordinate information (unless specified in a parameter), and assuming a lat-lon pseudoprojection unless otherwise indicated.
#'
#' Default behavior built into the function means that nine times out of ten, you can create your spatial object with a simple spdf(obj) call, instead of the standard sp:SpatialPointsDataFrame function call which requires three parameters. The function is also easy to use in magrittr pipe workflows. However, one time out of ten things could go horribly wrong if you're not careful!
#' @param data An object coercible into a data.frame
#' @param crs The proj4string argument for SpatialPointsDataFrame function. Defaults to WGS-84 lat-lon
#' @param llcols Names of the longitude and latitude (x and y coordinate) columns, in that order. If missing, looks first for column names starting with latitude and longitude, then for column names starting with lat and lon, then with x and y (case insensitive). If any of the lookups result in one hit or more than two hits, or all lookups result in no hits, the function will fail.
#' @param na.action What to do with missing or invalid values in coordinate columns
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
spdf <- function(data,crs=NULL,llcols=NULL,na.action=na.omit) {
  if(grepl('^Spatial',class(data))) {
    warning('Data is already of type Spatial*')
    return(data)
  }

  if(class(data) != 'data.frame')
    as.data.frame(data) -> data

  if(is.null(llcols)) {
    llcols <- unlist(sapply(c("^longitude","^latitude"),grep,tolower(names(data))))

    if(length(llcols)!=2)
      llcols <- unlist(sapply(c("^lon","^lat"),grep,tolower(names(data))))

    if(length(llcols)!=2)
      llcols <- unlist(sapply(c("x","y"),function(str) { which(str == tolower(names(data))) }))

    if(length(llcols)!=2)
      llcols <- unlist(sapply(c("^x","^y"),grep,tolower(names(data))))

    if(length(llcols)<2)
      stop("could not identify coordinate columns, need to start with lat, lon or x, y")

    if(length(llcols)>2)
      stop("could not identify coordinate columns, too many starting with lat, lon or x, y")
  }

  if(anyNA(data[,llcols])) warning("NAs in coordinates")
  data[row.names(na.action(data[,llcols])),] -> data

  if(is.null(crs)) {
    crs <- llproj
    warning("assuming lat-lon WGS84 coordinate system")
  }

  if(class(crs) == "character")
    crs <- CRS(crs)

  return(sp::SpatialPointsDataFrame(coords = data[,llcols],data=data,proj4string = crs))
}

#' Geocode multiple locations, retrying in case of errors
#'
#' @param locs Locations to geocode; a vector of place names or a data.frame or matrix with column called 'name' to be used as queries
#' @param repeats How many times to try before returning error
#' @param source See documentation for ggmap::geocode
#' @return A data.frame with columns name, lat and lon
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
geocode.loop <- function(locs,repeats=10, source=c("google", "dsk"), ...) {
  require(magrittr)
  if(locs %>% dim %>% length <= 1) {
    data.frame(name=locs,lon=NA,lat=NA,stringsAsFactors = F) -> locs
  } else {
    locs <- as.data.frame(locs)
    if(is.null(locs$lat)) locs$lat <- NA
    if(is.null(locs$lon)) locs$lon <- NA
  }

  counter=1
  while(anyNA(locs$lat) & counter <= repeats * length(source)) {
    paste("geocoding",sum(is.na(locs$lat)),"locations") %>% print
    locs[is.na(locs$lat),c('lon','lat')] <-
      ggmap::geocode(locs[is.na(locs$lat),'name'],
              source=source[ceiling(counter / repeats)], ...)
    counter = counter + 1
  }
  return(locs)
}

#' Order factor levels based on another variable
#'
#' Create a factor from a factor or character string, with levels ordered based on another variable.
#' @param var Factor or character vector with the values of the output factor.
#' @param by A vector that can be ordered with base::order.
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
orderfactor <- function(var, by) {
  var = factor(var, levels = var[rev(order(by))])
}



#' Get a list of dates, with ability to specify year
#'
#' Given a vector of dates or strings convertible to dates with as.Date, return all the dates between the earliest and the latest date in the vector.
#' @param year If specified, convert to corresponding dates of the given year
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
getdatelist <- function(dateinput,year=NULL) {
  dateinput <- as.Date(dateinput)
  seq(min(dateinput),max(dateinput),by='day') -> dateoutput
  if(!is.null(year))
    lubridate::year(dateoutput) <- year
  return(dateoutput)
}

#' Determine the appropriate UTM proj4string to use for a location, or get proj4string for a zone
#'
#' @param zone Zone number; if missing, this is determined from loc
#' @param hem Hemisphere. "S" or "s" for southern. If missing, this is determined from loc.
#' @param loc Location. Either a Spatial* object or a vector with lon and lat values.
#' @param units Either "m" or "km"; defaults to km.
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
getUTMproj <- function(zone=NULL,hem=NULL,loc=NULL,units="km") {

  if(!is.null(loc) & (!is.null(zone) | !is.null(hem)))
    warning("using explicit zone / hemisphere settings to override coordinate input")

  if(!is.null(loc) & grepl("Spatial",class(loc))) {
    require(sp)
    if(proj4string(loc) != llproj) loc <- spTransform(loc,CRS(llproj))
    ll <- colMeans(coordinates(loc))
  }

  if(!is.null(loc) & !grepl("Spatial",class(loc))) {
    warning("numeric location input - assuming lon-lat coordinate system")
    ll <- loc
  }

  if(is.null(zone))
    zone <- floor((ll[1] + 180)/6) %% 60 + 1

  if(is.null(hem)) {
    southhemi <- ll[2] < 0
  } else southhemi <- tolower(substr(hem,1,1)) == "s"



  paste0("+proj=utm +datum=WGS84 +no_defs +zone=",zone,ifelse(southhemi," +south ","")," +units=",units)
}

#' Round a number down to specified number of significant digits
#'
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
sigfloor <- function(x,sigdig=1) {
  mag <- 10^floor(log10(x)-sigdig+1)
  return(floor(x/mag)*mag)
}

#' Clean up the global environment to reduce memory use
#'
#' The function shows total memory use and lists all objects with a size larger than a threshold, giving the user the option to delete each one.
#' @param protect Vector of object names that user is not allowed to delete
#' @param threshold Objects with size below threshold are not offered for deletion. Megabytes. Defaults to 1.
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
cleanup <- function(protect,threshold=1) {
  sort( sapply(ls(envir = .GlobalEnv),function(x){object.size(get(x))}),decreasing = T)/1e6 -> memsizes
  trash <- rep(F,length(memsizes))

  for(o in which(memsizes>threshold)) {
    print(paste0("total memory size: ",round(sum(memsizes[!trash]),1),"Mb"))
    print(paste0(ifelse(o==1,"","next "),"largest object: ",names(memsizes)[o],", ",round(memsizes[o],1),"Mb"))
    readline("remove y/n/q ") -> reply
    if(substr(reply,1,1) %in% c("", "q")) break()

    if(substr(reply,1,1) == "y") trash[o] <- T
  }

  if(sum(trash==T)==0) return()

  rmObjs <- names(memsizes[trash])
  print(paste("removing",paste(rmObjs,collapse=" ")))
  readline("continue y/n ") -> reply2
  #rmObjs <- names(memsizes[memsizes>10e6 & !(names(memsizes) %in% c(protectObjects,"total"))])
  if(substr(reply2,1,1) == "y") rm(list=rmObjs, envir = .GlobalEnv)
}

#' Plot all pch symbols
#'
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
pchcheat <- function() {
  plot(x=rep(1:5,5),y=sapply(1:5,rep,5),pch=1:25)
  text(x=rep(1:5,5),y=sapply(1:5,rep,5),1:25,pos=c(rep(3,5),rep(1,20)))
}

#' Plot all named colors
#'
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
colcheat <- function(cl=colors(), bg = "grey",
                     cex = 0.75, rot = 30) {
  m <- ceiling(sqrt(n <-length(cl)))
  length(cl) <- m*m; cm <- matrix(cl, m)
  require("grid")
  grid.newpage(); vp <- viewport(w = .92, h = .92)
  grid.rect(gp=gpar(fill=bg))
  grid.text(cm, x = col(cm)/m, y = rev(row(cm))/m, rot = rot,
            vp=vp, gp=gpar(cex = cex, col = cm,font=2))
}

#' Normalize data.frame column values to 100 in a specified column
#'
#' Normalizes all values in data.frame columns based on values in one column. Most natural use is to set values to 100 on a basedate for plotting.
#' @param cols Integer or character vector to indicate columns on which normalization should be performed. Defaults to all numeric columns.
#' @param basedate Integer or character to indicate column which is used as the basis for normalization. If several values are given, the second one is used in case value in the first one is not finite.
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
norm100 <- function(df,cols=which(sapply(df, is.numeric)),basedate) {
  if(is.null(df$date)) stop('date column not present')

  df[,cols] <-
    lapply(df[,cols],function(x) {
      x.base = x[which(df$date == basedate[1])]
      if(!is.finite(x.base)) x.base <- x[which(df$date == basedate[2])]
      x / x.base  * 100
    } )
  return(df)
}

#' Normalize data.frame column values to zero in a specified column
#'
#' Normalizes all values in data.frame columns based on values in one column. Most natural use is to set values to 0 on a basedate for plotting.
#' @param cols Integer or character vector to indicate columns on which normalization should be performed. Defaults to all numeric columns.
#' @param basedate Integer or character to indicate column which is used as the basis for normalization. If several values are given, the second one is used in case value in the first one is not finite.
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
norm.0 <- function(df,cols,basedate) {
  if(is.null(df$date)) stop('date column not present')
  df[,cols] <-
    lapply(df[,cols],function(x) {
      x.base = x[which(df$date == basedate[1])]
      if(!is.finite(x.base) & length(basedate)>1)
        x.base <- x[which(df$date == basedate[2])]
      x - x.base
    } )
  return(df)
}

#' Capitalize the first letter of each word
#'
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
capitalize.first <- function(x) {
  sapply(x, function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")
  } ) -> x.clean
  return(x.clean)
}

#' Download multiple files, making a specified number of attempts before giving up
#'
#' @param urls URLs to download
#' @param destfiles Filenames for downloaded files - defaults to filename in the URL
#' @param tries How many tries to make before giving up
#' @param verbose Print an update after each file?
#' @param ... Passed on to download.file
#' @return A data.frame reporting whether each file was already on disk, downloaded successfully, or could not be downloaded
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
try.download <- function(urls,
                         destfiles=gsub('.*/','',urls),
                         tries=10,mode='wb',
                         overwrite=T,
                         quiet=T,
                         verbose=T,...) {

  outputs <- data.frame(url=urls,file=destfiles,status='missing',stringsAsFactors = F)

  if(!overwrite)
    outputs$status <- ifelse(file.exists(destfiles),'already on disk','missing')

  if(file.exists('try.download.temp')) file.remove('try.download.temp')
  for(i in 1:nrow(outputs)) {

    if(outputs[i,'status'] == 'missing') {
      tries <- 0
      while(!file.exists('try.download.temp') & tries < 10) {
        try(download.file(outputs[i,'url'],destfile = 'try.download.temp',
                          quiet=quiet,mode=mode, ...))
        tries = tries+1
      }

      if(file.exists('try.download.temp')) {
        file.rename('try.download.temp',outputs[i,'file'])
        outputs[i,'status'] <- 'downloaded'
      }
    }

    if(verbose)
      print(paste('file',outputs[i,'file'],outputs[i,'status'],
                  paste0('(',i,' of ',nrow(outputs),')')))
  }

  return(outputs)
}

#' Simple outlier detection based on standard deviation
#'
#' @param SDs How many standard deviations from mean signify an outlier. Defaults to 10.
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
is.outlier <- function(x, SDs=10,na.rm=F) {
  abs((x - mean(x, na.rm=na.rm)) / sd(x, na.rm = na.rm)) -> devs
  return(is.na(devs) | devs > SDs)
}


#' Calculate mean with maximum amount of NA values specified
#'
#' A wrapper around the base mean function that returns NA when there are more than maxna NAs in input
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
mean.maxna <- function(x,maxna) {
  if(sum(is.na(x))>maxna) { return(NA)
  } else return(mean(x,na.rm=T))
}

#' Calculate statistics with maximum amount of NA values specified
#'
#' A wrapper statistics functions that returns NA when there are more than maxna NAs in input
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
fun.maxna <- function(x,fun,maxna) {
  if(sum(is.na(x))>maxna) { return(NA)
  } else return(fun(x,na.rm=T))
}


#' Specify a transparent color with name and alpha value
#'
#' @param colorname A color name recognized by col2rgb
#' @param alpha An alpha value between 0 and 1
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
col.a <- function(colorname,alpha) {
  require(magrittr)
  colorname %>% col2rgb %>% unlist %>% divide_by(255) -> cn
  rgb(cn[1],cn[2],cn[3],alpha)
}

#' Cheatsheet for lattice pos argument
#'
#' Print a plot to show text label positions corresponding to different lattice pos arguments (1 to 4)
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
poscheat <- function() {
  require(magrittr)
  c(-1:1) %>% data.frame(lat=., lon=.,z=.) %>% spdf -> t1
  sp::spplot(t1,zcol='z') +
    latticeExtra::layer(sp::sp.text(sp::coordinates(t1[c(2,2,2,2),]),
                  txt=as.character(1:4), pos=1:4))
}

#' A convenience function to use the wesanderson Darjeeling1 palette in a ggplot
#'
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
darjeeling <- function(ordering=T,...) {
  ggplot2::scale_color_manual(values=wesanderson::wes_palette("Darjeeling1")[ordering],...)
}

#' A convenience function to initialize a new png plot with defaults for width, height and res
#'
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
quickpng <- function(file, width=2000, height=1500, res=300, ...) {
  png(filename=file, width=width, height=height, res=res, ...)
}

#' Load administrative area boundaries
#'
#' Requires that the GADM RDS files or shapefiles exist in a path specified in either the environment variable GISpath or in the function parameter path. Tries to load first the RDS file and secondarily the shp file. If RDS file does not exist, it will convert the shp into RDS and write to disk.
#'
#' The required files can be downloaded from https://gadm.org/
#' @param level The level of admin boundaries to load.
#' @param path Path to the file
#' @return SpatialPolygonsDataFrame
#' @author Lauri Myllyvirta \email{lauri.myllyvirta@@greenpeace.org}
#' @export
getadm <- function(level=0, path=Sys.getenv('GISpath')) {
  require(sp)

  if(path=="")
    path='~/GIS/boundaries/'

  inF <- paste0(path, '/gadm28_adm',level,'.RDS')

  if(file.exists(inF)) {
    readRDS(inF)
  } else {
    raster::shapefile(gsub('\\.RDS','.shp',inF),
                      encoding='UTF-8', use_iconv=TRUE) -> admout
    writeRDS(admout, inF)
    return(admout)
  }
}
