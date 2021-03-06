% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{getdatelist}
\alias{getdatelist}
\title{Get a list of dates, with ability to specify year}
\usage{
getdatelist(dateinput, year = NULL)
}
\arguments{
\item{year}{If specified, convert to corresponding dates of the given year}
}
\description{
Given a vector of dates or strings convertible to dates with as.Date, return all the dates between the earliest and the latest date in the vector.
}
\author{
Lauri Myllyvirta \email{lauri.myllyvirta@greenpeace.org}
}
