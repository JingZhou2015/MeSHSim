\name{docSim}
\alias{docSim}
\title{Similarity between articles}
\description{  Calculate the similarity between two articles.  }
\usage{
    docSim(pmid1, pmid2, method="SP", frame="node", major=FALSE, env=NULL)
}
\arguments{
    \item{pmid1, pmid2}{  pmids of two articles whose similarity is needed to be calculated.  }
    \item{method}{  similarity measurment method, see Details for available methods.  }
    \item{frame}{  available options are node and heading, decide whether using node-based or heading-based methods.  }
    \item{env}{  the dataset to use.  }
    \item{major}{  whether the calculation should only be based on major terms  }
}
\value{
    Semantic similarity between two MEDLINE documents, whose value is between 0 and 1.
}
\examples{
    docSim("1111113","1111111")
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
\note{
    Network connection is required for using this function.  
}
\keyword{MeSH}
