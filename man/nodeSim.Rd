\name{nodeSim}
\alias{nodeSim}
\title{Similarity between nodes}
\description{
    Calculate similarity between two MeSH nodes.
}
\usage{
    nodeSim(node1, node2, method="SP", frame="node", env=NULL)
}
\arguments{
    \item{node1, node2}{  two nodes or two lists of nodes  }
    \item{method}{  similarity measurment method, see Details for available methods.  }
    \item{frame}{  available options are node and heading, decide whether using node-based or heading-based methods.  }
    \item{env}{  the dataset to use.  }
}
\value{
    Semantic similarity between two MeSH nodes, whose value is between 0 and 1.
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
\examples{
    nodeSim("B03.440.400.425.340.590", "B03.440.400.425.117.800.200")
}
