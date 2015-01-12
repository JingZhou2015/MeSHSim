\name{nodeInfo}
\alias{nodeInfo}
\title{Details of nodes}
\description{
    This function returns a tree contains the given node.
}
\usage{
    nodeInfo(node, brief, env=NULL)
}
\arguments{
    \item{node}{  a node name.  }
    \item{brief}{brief model for nodeInfo}
    \item{env}{  the dataset to use.  }
}
\value{
    Hierarchy information of node a
}
\examples{
	nodeInfo("B03.440.400.425.127")
    nodeInfo("B03.440.400", brief=TRUE)
}
