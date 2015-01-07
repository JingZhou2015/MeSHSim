\name{mheadingSim}
\alias{mheadingSim}
\title{Similarity between heading lists}
\description{
    Calculate similarity matrix between two heading lists. 
}
\usage{
    mheadingSim(headingList1, headingList2, method="SP", frame="node", env=NULL)
}
\arguments{
    \item{headingList1, headingList2}{  two headings or two lists of headings  }
    \item{method}{  similarity measurment method, see Details for available methods.  }
    \item{frame}{  available options are node and heading, decide whether using node-based or heading-based methods.  }
    \item{env}{  the dataset to use.  }
}
\value{
    Semantic similarity matrix between two MeSH heading lists.
}
\examples{
    headingList1<-c("Body Regions", "Abdomen", "Abdominal Cavity")
    headingList2<-c("Lumbosacral Region", "Body Regions")
    mheadingSim(headingList1,headingList2)
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