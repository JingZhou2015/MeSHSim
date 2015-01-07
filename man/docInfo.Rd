\name{docInfo}
\alias{docInfo}
\title{Get details of documents}
\description{  Fetch information of a given article from PubMed.  }
\usage{
    docInfo(pmid, verbose=FALSE, major=FALSE)  
}
\arguments{
    \item{pmid}{  pmid of the desired article.  }
    \item{verbose}{  whether the title and abstract of the article should be print out.  }
    \item{major}{  whether only major terms should be returned.  }
}
\value{
    Document information of given PMID including titile, abstract, MeSH headings
}
\examples{
    docInfo("1111111")
}
\note{
    Network connection is required for using this function.  
}
\keyword{MeSH}
