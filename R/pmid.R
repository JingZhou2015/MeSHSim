docSim<-function(pmid1, pmid2, method="SP", frame="node", major=FALSE, env=NULL){
    env=check_and_set_data(env);
    ha<-get_headings_by_PMID(pmid1, FALSE, major)
    hb<-get_headings_by_PMID(pmid2, FALSE, major)
    if(length(ha)==0||length(hb)==0)
    	return(0.0);
    return(`_docSim`(ha, hb, method, frame, env))
}
docInfo<-function(pmid, verbose=FALSE, major=FALSE){
    return(get_headings_by_PMID(pmid,verbose,major))
}
