\name{headingSim}
\alias{headingSim}
\title{Similarity between headings}
\description{
    Calculate similarity between two headings.
}
\usage{
    headingSim(heading1, heading2, method="SP", frame="node", env=NULL)
}
\arguments{
    \item{heading1, heading2}{  two headings or two lists of headings  }
    \item{method}{  similarity measurment method, see Details for available methods.  }
    \item{frame}{  available options are node and heading, decide whether using node-based or heading-based methods.  }
    \item{env}{  the dataset to use.  }
}
\value{
    Semantic similarity between two MeSH headings, whose value is between 0 and 1.
}
\examples{
    headingSim("Lumbosacral Region", "Body Regions")
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
