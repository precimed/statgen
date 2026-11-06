\name{verbosity}
\alias{get_verbosity}
\alias{set_verbosity}
\title{Get or Set Runtime Verbosity}
\usage{
get_verbosity()
set_verbosity(level)
}
\arguments{
\item{level}{Character scalar, either \code{"quiet"} or \code{"info"}.}
}
\value{
\code{get_verbosity()} returns the current verbosity level. \code{set_verbosity()}
returns the selected level invisibly.
}
\description{
Controls the package-wide statgen runtime verbosity setting. Use
\code{"quiet"} to suppress progress messages and loading summaries.
Warnings remain visible.
}
\seealso{\link{allele_matching}}
\examples{
old <- get_verbosity()
set_verbosity("quiet")
get_verbosity()
set_verbosity(old)
}
