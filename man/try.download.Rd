% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{try.download}
\alias{try.download}
\title{Download multiple files, making a specified number of attempts before giving up}
\usage{
try.download(urls, destfiles = gsub(".*/", "", urls), tries = 10,
  mode = "wb", overwrite = T, quiet = T, verbose = T, ...)
}
\arguments{
\item{urls}{URLs to download}

\item{destfiles}{Filenames for downloaded files - defaults to filename in the URL}

\item{tries}{How many tries to make before giving up}

\item{verbose}{Print an update after each file?}

\item{...}{Passed on to download.file}
}
\value{
A data.frame reporting whether each file was already on disk, downloaded successfully, or could not be downloaded
}
\description{
Download multiple files, making a specified number of attempts before giving up
}
\author{
Lauri Myllyvirta \email{lauri.myllyvirta@greenpeace.org}
}
