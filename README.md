MeSHSim
=======

MeSHSim is an R package for semantic similarities calculation among MeSH headings and MEDLINE documents. As shown in table 1 Five path-based measures and four Information Content (IC)- based measures are implemented in MeSHSim. It also supports querying the hierarchy information of a MeSH heading and information of a given document including title, abstraction and MeSH headings. The package can be easily integrated into pipelines for other biomedical text analysis tasks. 

## Installation Instructions

### R Package Dependencies
MeSHSim needs three R packages: bitops, XML and RCurl, where bitops, used by RCurl, supports bitwise operations of integer vectors, XML supports reading XML documents and RCurl is to fetch document information from PubMed. They are all freely available at CRAN (Comprehensive R Archive Network).

#### bitops
Package bitops is available at CRAN. To install it call:
<code>
	install.packages("bitops")
</code>

#### XML
Package XML is available at CRAN. To install it call:
<code>
	install.packages("XML")
</code>

####RCurl
For Linux users, make sure you Linux has library libcurl installed. Check out:

<code>
	locate libcurl
	locate curl-config
</code>

If there are not found, most Linux user will be able to fix this by running:
<code>
    sudo apt-get install libcurl4-openssl-dev
</code>

You will now be able to install R package RCurl. In R console:
<code>
    install.packages("RCurl")
</code>


### MeSHSim Installation
If everything went well you will be able to install the MeSHSim package.

### From github (latest version, temporary)
I R console, type:
```R
library(devtools)
install_github("radcheb/MeSHSim")
```

#### From Bioconductor
In R console, type:
<code>
	source("http://bioconductor.org/biocLite.R")
	bioLite("MeSHSim")
</code>

#### Install offline
In OS console, checkout the source codes from R,
<code>
	git clone git@github.com:JingZhou2015/MeSHSim.git
</code>

Build MeSHSim package in OS console,
<code>
	R CMD build MeSHSim
</code>

Install MeSHSim package in OS console,
<code>
	R CMD install MeSHSim_{#version}.tar.gz
</code>
where {#version} is the version number of MeSHSim package

Load MeSHSim pacakge in R console
<code>
	library('MeSHSim')
</code>
