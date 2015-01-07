\name{mnodeSim}
\alias{mnodeSim}
\title{Similarity between node lists}
\description{
    Calculate similarity matrix between two MeSH node lists.
}
\usage{
    mnodeSim(nodeList1, nodeList2, method="SP", frame="node", env=NULL)
}
\arguments{
    \item{nodeList1, nodeList2}{  two nodes or two lists of nodes  }
    \item{method}{  similarity measurment method, see Details for available methods.  }
    \item{frame}{  available options are node and heading, decide whether using node-based or heading-based methods.  }
    \item{env}{  the dataset to use.  }
}
\value{
    Semantic similarity matrix between two MeSH node lists.
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
    nodeList1<-c("B03.440.400.425.340.590", "B03.440.400.425.117.800.200")
    nodeList2<-c("B03.440.400.425.340.590", "B03.440.400.425.117.800.200", "B03.440.400.425.127.100")
    mnodeSim(nodeList1,nodeList2)
}
