\name{termInfo}
\alias{termInfo}
\title{Details of MeSH terms}
\description{
    This function returns a tree contains the given term.
}
\usage{
    termInfo(term, brief,env=NULL)
}
\arguments{
    \item{term}{  a term name.  }
    \item{brief}{whether to retrive breif tree information of MeSH term}
    \item{env}{  the dataset to use.  }
}
\value{
    Hierarchy information of a given term
}
\examples{
    termInfo("Body Regions")
}