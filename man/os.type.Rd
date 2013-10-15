\name{os.type}
\alias{os.type}
\title{returns TRUE if correct operating system is passed as an argument}
\usage{
  os.type(type = c("linux", "mac", "windows", "else"))
}
\arguments{
  \item{type}{character string of operating system type}
}
\value{
  logical
}
\description{
  returns TRUE if correct operating system is passed as an
  argument
}
\note{
  \code{os.type} is intended to be used internally and
  based on the same function in INLA
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

