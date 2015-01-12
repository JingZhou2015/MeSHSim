eutils_base<-"http://www.ncbi.nlm.nih.gov/entrez/eutils/"
eutils_fetch<-function(a){
    url<-paste(eutils_base, "esearch.fcgi?db=pubmed&usehistory=y&term=", a, sep="")
    doc<-getURL(url)
    doc<-xmlRoot(xmlTreeParse(doc))
    if(xmlValue(doc[["Count"]])!=1){
        print("Inapropriate PMID")
        return(NULL)
    }
    url<-paste(eutils_base, "efetch.fcgi?retmode=xml&db=pubmed&query_key=", xmlValue(doc[["QueryKey"]]), "&WebEnv=", xmlValue(doc[["WebEnv"]]), sep="")
    doc<-getURL(url)
    doc<-xmlRoot(xmlTreeParse(doc))
    return(doc)
}
eutils_parse<-function(a, verbose){
    xml_mesh_list<-a[["PubmedArticle"]][["MedlineCitation"]][["MeshHeadingList"]]$children
    if(verbose){
        print(paste("Title", xmlValue(a[["PubmedArticle"]][["MedlineCitation"]][["Article"]][["ArticleTitle"]]), sep=": "));
        print(paste("Abstract", xmlValue(a[["PubmedArticle"]][["MedlineCitation"]][["Article"]][["Abstract"]][["AbstractText"]]), sep=": "));
        print("MeSH Headings:");
    }
    headings<-c()
    major_headings<-c()
    get_headings<-function(b){
        xmlValue(b[["DescriptorName"]])
    }
    get_major_headings_bool<-function(b){
		i<-1;
		is_major = FALSE;
		while(i <= length(b)) {
			if(xmlAttrs(b[[i]])[["MajorTopicYN"]] == "Y") {
				is_major = TRUE;
				break;
			}
			i<-i+1;
		}
        return(is_major)
    }
    headings<-as.vector(lapply(xml_mesh_list, get_headings), mode="character")
    major_headings<-headings[as.vector(lapply(xml_mesh_list, get_major_headings_bool), mode="logical")]
	headings<-headings[!headings=="Female"]
	headings<-headings[!headings=="Male"]
	major_headings<-major_headings[!major_headings=="Female"]
	major_headings<-major_headings[!major_headings=="Male"]
    return(list(all=headings, major=major_headings))
}
get_headings_by_PMID<-function(a, verbose=FALSE, ismajor=FALSE){
    doc<-eutils_fetch(a)
    tmp<-eutils_parse(doc, verbose)
    if(ismajor) return(tmp$major)
    return(tmp$all)
}
