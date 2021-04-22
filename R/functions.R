#' Load default libraries
#'
#' @param install.missing Attempt to install missing libraries? Default=T
#' @param custom additional libraries to load
#' @param check.updates check for updates
#' @param defaults all default libraries. Default=T
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
loadlibs <- function(general=defaults,
                     plotting=defaults,
                     geo=defaults,
                     xls=defaults,
                     geocoding=defaults,
                     defaults=T,
                     custom=NULL,
                     install.missing=T,
                     check.updates=T,
                     online=F,
                     crea=T) {
  if(is.null(online))
    online <- (system2("ping", "baidu.com", stderr = FALSE, stdout = FALSE) == 0)

  reqlibs <- custom
  if(general)
    reqlibs <- c(reqlibs, "magrittr","plyr","tidyverse",
                 "reshape2","data.table","lubridate","zoo",
                 "pbapply", "countrycode")

  if(plotting)
    reqlibs <- c(reqlibs,"lattice","RColorBrewer","ggplot2","ggthemes", "ggsci")

  if(geo)
    reqlibs <- c(reqlibs,"sp","raster","rasterVis","rworldmap","rgdal","sf","geosphere","gstat", "rgeos", "ggspatial")

  if(xls) reqlibs <- c(reqlibs,"readxl")

  if(geocoding) reqlibs <- c(reqlibs,"ggmap")

  unique(reqlibs) -> reqlibs

  if(check.updates & online) {
    if(nrow(old.packages()) > 0)
      readline('packages need to be updated, continue y/n?') -> r

    if(tolower(r) == 'y')
      update.packages(ask=F)
  }

  inst.libs <- row.names(installed.packages())
  needed.libs <- reqlibs[!(reqlibs %in% inst.libs)]

  if(install.missing & online & length(needed.libs) > 0) {
    install.packages(needed.libs)
  } else {
    if(length(needed.libs)>0) warning(paste(paste(needed.libs,sep=','),'not installed!'))
    reqlibs <- reqlibs[reqlibs %in% inst.libs]
  }

  lapply(reqlibs, library, character.only = T)

  if(crea) source(boxpath('CREA/CREAtheme.R'))
}

#' Proj4string for standard latitude-longitude projection
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
llproj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

#' Constant to convert molec/cm2 to Dobson Units
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
Dobson.unit = 2.69e16

#' Return true for all values of x that are NOT included in y
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
'%notin%' <- function(x,y)!('%in%'(x,y))

#' Return a vector of the values of x that are included in y
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
'%whichin%' <- function(x,y) x[x %in% y]

#' Return a vector of the values of x that are NOT included in y
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
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
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
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

#' Order factor levels based on another variable
#'
#' Create a factor from a factor or character string, with levels ordered based on another variable.
#' @param var Factor or character vector with the values of the output factor.
#' @param by A vector that can be ordered with base::order.
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
orderfactor <- function(var, by) {
  var = factor(var, levels = var[rev(order(by))])
}

#' Get ISO3 codes for country names
#'
#' Shorthand for the countrycode::countrycode function.
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
get_iso3 = function(x, custom_match=NULL) {
  require(countrycode)
  countrycode(x, 'country.name', 'iso3c', custom_match = c(Kosovo='XKX', custom_match))
}


#' Get a list of dates, with ability to specify year
#'
#' Given a vector of dates or strings convertible to dates with as.Date, return all the dates between the earliest and the latest date in the vector.
#' @param year If specified, convert to corresponding dates of the given year
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
getdatelist <- function(dateinput,year=NULL) {
  dateinput <- as.Date(dateinput)
  seq(min(dateinput),max(dateinput),by='day') -> dateoutput
  if(!is.null(year))
    lubridate::year(dateoutput) <- year
  return(dateoutput)
}

#' Strip non-numeric characters and convert to numeric
#'
#' @param selection Vector of positions of numbers to return. TRUE to return all. Defaults to the first number found.
#' @return If length of selection is 1, the function returns a vector; a list otherwise.
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
force.numeric = function(x, selection=1) {
  x %>% gsub(',', '', .) %>% gsub('[^0-9\\.\\-]', ' ', .) %>% strsplit(" ") -> x.out
  x.out %<>% lapply(function(x){x %>% as.numeric %>% subset(!is.na(.)) %>% '['(selection)}) %>%
    lapply(function(x){if(is.null(x)){return(NA)}else return(x)})
  if(length(selection)==1 & !(is.logical(selection) && selection)) x.out %<>% unlist
  x.out
}

#' Calculate the mode of a vector
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
statmode <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  ux <- unique(x)
  ux[which.max(table(match(x, ux)))]
}

#' Replace NA values in a vector with the corresponding values of another vector
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
na.cover <- function(x, x.new) { ifelse(is.na(x), x.new, x) }

#' Determine the appropriate UTM proj4string to use for a location, or get proj4string for a zone
#'
#' @param zone Zone number; if missing, this is determined from loc
#' @param hem Hemisphere. "S" or "s" for southern. If missing, this is determined from loc.
#' @param loc Location. Either a Spatial* object or a vector with lon and lat values.
#' @param units Either "m" or "km"; defaults to km.
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
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
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
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
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
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
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
pchcheat <- function() {
  plot(x=rep(1:5,5),y=sapply(1:5,rep,5),pch=1:25)
  text(x=rep(1:5,5),y=sapply(1:5,rep,5),1:25,pos=c(rep(3,5),rep(1,20)))
}

#' Plot all named colors
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
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
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
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
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
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
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
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
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
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
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
is.outlier <- function(x, SDs=10,na.rm=F) {
  abs((x - mean(x, na.rm=na.rm)) / sd(x, na.rm = na.rm)) -> devs
  return(is.na(devs) | devs > SDs)
}


#' Calculate mean with maximum amount of NA values specified
#'
#' A wrapper around the base mean function that returns NA when there are more than maxna NAs in input
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
mean.maxna <- function(x,maxna) {
  if(sum(is.na(x))>maxna) { return(NA)
  } else return(mean(x,na.rm=T))
}

#' Calculate statistics with maximum amount of NA values specified
#'
#' A wrapper statistics functions that returns NA when there are more than maxna NAs in input
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
fun.maxna <- function(x,fun,maxna) {
  if(sum(is.na(x))>maxna) { return(NA)
  } else return(fun(x,na.rm=T))
}


#' Specify a transparent color with name and alpha value
#'
#' @param colorname A color name recognized by col2rgb
#' @param alpha An alpha value between 0 and 1
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
col.a <- function(colorname,alpha) {
  require(magrittr)
  colorname %>% col2rgb %>% unlist %>% divide_by(255) -> cn
  rgb(cn[1,],cn[2,],cn[3,],alpha)
}

#' Cheatsheet for lattice pos argument
#'
#' Print a plot to show text label positions corresponding to different lattice pos arguments (1 to 4)
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
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
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
darjeeling <- function(ordering=T,...) {
  ggplot2::scale_color_manual(values=wesanderson::wes_palette("Darjeeling1")[ordering],...)
}

#' A convenience function to initialize a new png plot with defaults for width, height and res
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
quickpng <- function(file, width=2000, height=1500, res=300, ...) {
  png(filename=file, width=width, height=height, res=res, ...)
}

#' Convert coordinates in degrees-minutes[-seconds] to decimal
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
dms_to_dec <- function(x) {
  x %>%
    lapply(function(x) {
      xn = x %>% strsplit('[^0-9.,]+') %>% unlist %>% gsub(',', '.', .) %>% as.numeric
      if(is.na(xn[3])) xn[3] <- 0
      x_dec = xn[1]+xn[2]/60+xn[3]/3600
      if(grepl("S$|W$", x)) x_dec %<>% multiply_by(-1)
      return(x_dec)
    }) %>%
    unlist
}

#' Load administrative area boundaries
#'
#' Requires that the GADM RDS or shp files exist in a path specified in either the environment variable GISpath or in the function parameter path. Tries to load first the RDS file and secondarily the shp file. If RDS file does not exist, it will convert the shp into RDS and write to disk.
#'
#' The required files can be downloaded from https://gadm.org/
#' @param level The level of admin boundaries to load.
#' @param path Path to the file
#' @return SpatialPolygonsDataFrame
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
getadm <- function(level=0, res='full', version='36', path=Sys.getenv('GISpath')) {
  require(sp)
  resext=''
  if(res!='full') resext=paste0('_', res)
  inF <- paste0(GISpath, '/', 'gadm',version,'_',level,resext,'.RDS')

  if(!file.exists(inF) & res=='full') {
    message('RDS not found, trying shapefile')
    raster::shapefile(gsub('\\.RDS','.shp',inF),
                      encoding='UTF-8', use_iconv=TRUE) -> adm
    saveRDS(adm, inF)
  }

  if(!file.exists(inF) & res!='full') {
    message('simplifying the admin features')
    simplify_adm(level, res, version, path)
  }
  readRDS(inF)
}

#' Simplify administrative area boundaries
#'
#' Requires that the GADM RDS files exist in a path specified in either the environment variable GISpath or in the function parameter path. Tries to load first the RDS file and secondarily the shp file. If RDS file does not exist, it will convert the shp into RDS and write to disk.
#'
#' The required files can be downloaded from https://gadm.org/
#' @param level The level of admin boundaries to load.
#' @param path Path to the file
#' @return SpatialPolygonsDataFrame
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
simplify_adm = function(level=0, resname='low', version='36', path) {
  require(raster)
  require(sf)

  if(resname=='low') tol=.005
  if(resname=='coarse') tol=.025

  getadm(level, res='full', version) -> adm
  adm %>% st_as_sf() -> adm_sf

  adm_sf %>% st_simplify(dTolerance = tol) -> adm_sf_coarse

  outF = paste0(GISpath, '/', 'gadm',version,'_',level,'_',resname,'.shp')
  st_write(adm_sf_coarse, outF, overwrite=T, append=F)
  shapefile(outF) -> adm_coarse

  saveRDS(adm_coarse, gsub('\\.shp', '.RDS', outF))
}


#' Repair and make unique names with one command
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
namerep=function(x){x %>% make.names %>% make.unique}

#' Extract parameter values from proj4 strings
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
get_crs_par = function(crsstring, parname) {
  search_string = paste0('.*\\+',parname,'=')
  if(grepl(search_string, crsstring)) {
    crsstring %>% gsub(search_string, '', .) %>% gsub(' .*', '', .)
  } else return(as.character(NA))
}

#' Repair UTM projections with km as units
#'
#' New versions of GDAL convert UTM projections with km as units into generic transverse Mercator, which causes no end of problems.
#'
#' This function takes an object that has a crs and converts it back to UTM/
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
fixproj = function(x) {
  crs0 = crs(x)
  if(grepl("\\+proj=tmerc", crs0)) {
    UTMZ = get_crs_par(crs0, 'lon_0') %>% as.numeric %>% add(177) %>% divide_by(6) %>% add(1)
    UTMH = ifelse(get_crs_par(crs0, 'y_0')=='0', '', ' +south')
    u = get_crs_par(crs0, 'units')
    datum=get_crs_par(crs0, 'datum')
    if(is.na(datum)) datum='WGS84'
    crs(x)=paste0("+proj=utm +datum=",datum,
                  " +no_defs +zone=",UTMZ,UTMH,
                  " +units=",u)
  }
  return(x)
}

#' Load a population density raster
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
getpopdens <- function(path=Sys.getenv('GISpath')) {
  raster(paste0(path, '/','gpw-v4-population-density-adjusted-to-2015-unwpp-country-totals_2015.tif'))
}

#' Load a population count raster
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
getpopcount <- function(path=Sys.getenv('GISpath')) {
  raster(paste0(path, '/','~/GIS/population/gpw-v4-population-count_2015.tif'))
}

#' Load a population raster
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
getpop <- function(what='density') {
  if(what=='density') pop.out = getpopdens()
  if(what=='count') pop.out = getpopcount()
  return(pop.out)
}


#' Alias for dplyr::select
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
sel <- dplyr::select

#' Alternative for raster::projectExtent that works with UTM with +units=km
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
projectExtent2 = function(r, targetcrs, ...) {
  r %>% extent() -> bb
  bb %>% as('SpatialPolygons') -> bbpol
  crs(bbpol) = crs(r)
  bbpol %>% spTransform(crs(targetcrs)) %>% extent %>% raster(...) -> r.out
  crs(r.out) = targetcrs
  return(r.out)
}

#' Alternative for sp::spTransform that works with UTM with +units=km
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
spTransform2 = function(obj, targetcrs, ...) {
  if(grepl('\\+units=km', targetcrs)) obj %<>% spTransform(gsub('\\+units=km', '\\+units=m', targetcrs, ...))
  obj %>% spTransform(targetcrs, ...)
}

#' Reproject a spatial object, cropping to the extent of target raster
#'
#' This function reprojects the extent of the target raster, crops the spatial object to this extent, and then reprojects it to the crs of the target raster.
#'
#' Using standard functions raster::projectRaster or sp::spTransform directly causes an error when the object cannot be entirely represented in the target crs, f.ex. when reprojecting a global dataset to UTM. Cropping first is also faster.
#'
#' @param shapeobj Object to be reprojected
#' @param rasterobj Raster from which to take the crs and extent
#' @param expand Factor by which to expand the extent before cropping
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
cropProj <- function(shapeobj, rasterobj, expand=1.25, ...) {
  bb = rasterobj %>% projectExtent2(crs(shapeobj)) %>%
    extent %>% multiply_by(expand)
  shapeobj %<>% crop(bb)
  if(grepl("Raster", class(shapeobj))) {
    shapeobj %>% projectRaster(rasterobj, ...) %>% return
  } else shapeobj %>% spTransform2(crs(rasterobj)) %>% return
}

#' Group spatial points into clusters
#'
#' Cluster a set of points.
#'
#' @param sp A SpatialPoints* object or data.frame with coordinate information
#' @param distKM Maximum distance between points in the same cluster, in kilometers
#' @return Numeric vector giving the cluster for each row in the input data.frame
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
cluster <- function(sp, distKM) {
  require(geosphere)
  sp <- spdf(sp)
  hc <- sp %>% coordinates %>% distm %>% as.dist %>% hclust
  cutree(hc,h=distKM*1000)
}

#' Apply a function recursively
#'
#' Applies a function recursively until the result is of desired class.
#'
#' @param L a (nested) list to which to apply the function
#' @param f function to apply
#' @param targetclass expected class of the result, defaults to data.frame
#'
#' @author Lauri Myllyvirta \email{lauri@@energyandcleanair.org}
#' @export
recurse <- function (L, f, targetclass="data.frame",...) {
  if(inherits(L, targetclass) | !is.list(L)) f(L, ...)
  else lapply(L, recurse, f, ...)
}
