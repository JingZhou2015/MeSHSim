\name{headingSetSim}
\alias{headingSetSim}
\title{Similarity between two MeSH heading sets}
\description{
    Calculate similarity between two MeSH heading sets.
}
\usage{
    headingSetSim(headingSet1, headingSet2, method="SP", frame="node", env=NULL)
}
\arguments{
    \item{headingSet1, headingSet2}{ two lists of headings  }
    \item{method}{  similarity measurment method, see Details for available methods.  }
    \item{frame}{  available options are node and heading, decide whether using node-based or heading-based methods.  }
    \item{env}{  the dataset to use.  }
}
\value{
    Semantic similarity between two MeSH heading sets, whose value is between 0 and 1.
}
\examples{
    headingSet1<-c("Lumbosacral Region", "Body Regions")
    headingSet2<-c("Body Regions", "Abdomen", "Abdominal Cavity")
    headingSetSim(headingSet1,headingSet2,'SP','node')
}
\details{
    Available methods:
    SP: Shortest Path method,
    WL: Weighted Link method,
    WP: Wu and Palmer's method,
    LC: Leacock and Chodorow's method,
    Li: Li's method,
    Lord: Lord's method,
    Resnik: Resnik's method,
    Lin: Lin's method,
    JC: Jiang and Conrath's method.  
}
\keyword{MeSH}