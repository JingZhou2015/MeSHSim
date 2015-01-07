hierarchies<-c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "V", "Z")
new_term<-function(des="", ns=list(), ic=0){
    list(descriptor=des, nodes=ns)
}
new_node<-function(des="", n=c(), h="", m=c(0,0,0), wm=c(0,0,1/length(n)), ic=0){
    list(descriptor=des, node=n, hierarchy=h, max=m, wmax=wm)
}
new_tree<-function(){
    list(terms=c())
}
get_term_var_name<-function(a){
    paste("_mesh_term_descriptor_", a, sep="")
}
get_node_var_name<-function(a){
    paste("_mesh_node_", a, sep="")
}
get_node_var_name_v<-function(a, len=length(a)){
    if(len <= 1){
        return(get_node_var_name(a[[1]]))
    }
    b<-paste(a[2:len], collapse=".")
    get_node_var_name(b)
}
get_hierarchy_var_prefix<-function(a){
    paste(paste("_mesh_hierarchy_", a, sep=""), "_", sep="")
}
get_tree_path<-function(a){
    tmp<-strsplit(a, ".", fixed=TRUE)[[1]]
    c(strsplit(tmp[[1]], "", fixed=TRUE)[[1]][1], tmp)
}
update_node_max_length<-function(a, env){
    node_name<-get_node_var_name_v(a)
    if(!exists(node_name, envir=env)){
        print(c("error", node_name, "not found"))
        return(NULL)
    }
    node<-get(node_name, envir=env)
    node$max[[3]]<-max(node$max[[3]], node$max[[1]]+node$max[[2]])
    node$wmax[[3]]<-max(node$wmax[[3]], node$wmax[[1]]+node$wmax[[2]]+1/length(a))
    assign(node_name, node, envir=env)
    if(length(a)<=1){
        return()
    }
    parent_name<-get_node_var_name_v(a, length(a)-1)
    if(!exists(parent_name, envir=env)){
        print(c("error", parent_name, "not found"))
        return(NULL)
    }
    parent<-get(parent_name, envir=env)
    if(node$max[[3]] > parent$max[[3]]){
        parent$max[[3]]<-node$max[[3]]
    }
    if(node$max[[1]]+1 > parent$max[[1]]){
        parent$max[[2]]<-parent$max[[1]]
        parent$max[[1]]<-node$max[[1]]+1
    }else if(node$max[[1]]+1 > parent$max[[2]]){
        parent$max[[2]]<-node$max[[1]]+1
    }
    if(node$wmax[[3]] > parent$wmax[[3]]){
        parent$wmax[[3]]<-node$wmax[[3]]
    }
    if(node$wmax[[1]]+1/length(a) > parent$wmax[[1]]){
        parent$wmax[[2]]<-parent$wmax[[1]]
        parent$wmax[[1]]<-node$wmax[[1]]+1/length(a)
    }else if(node$wmax[[1]]+1/length(a) > parent$wmax[[2]]){
        parent$wmax[[2]]<-node$wmax[[1]]+1/length(a)
    }
    assign(parent_name, parent, envir=env)
    return()
}
update_hierarchy_max_length<-function(h, max, wmax, env){
    #Update max length as used in method SP
    name<-paste(get_hierarchy_var_prefix(h), "MAX_LENGTH", sep="")
    if(!exists(name, envir=env))
        assign(name, max[[3]], envir=env)
    else
        if(max[[3]]>get(name, envir=env))
            assign(name, max[[3]], envir=env)
    #Update max wlength as used in method WL
    wname<-paste(get_hierarchy_var_prefix(h), "MAX_WLENGTH", sep="")
    if(!exists(wname, envir=env))
        assign(wname, wmax[[3]], envir=env)
    else
        if(wmax[[3]]>get(wname, envir=env))
            assign(wname, wmax[[3]], envir=env)

    #Update max depth of MeSH as used in method LC
    dname<-"_mesh_hierarchy_ANY_MAX_DEPTH"
    if(!exists(dname, envir=env))
        assign(dname, max[[1]], envir=env)
    else
        if(max[[1]]>get(dname, envir=env))
            assign(dname, max[[1]], envir=env)
}
get_hierarchy_max_length<-function(a, env){
    name<-paste(get_hierarchy_var_prefix(a), "MAX_LENGTH", sep="")
    get(name, env)
}
get_hierarchy_max_wlength<-function(a, env){
    name<-paste(get_hierarchy_var_prefix(a), "MAX_WLENGTH", sep="")
    get(name, env)
}
get_hierarchy_max_ic<-function(a, env){
    name<-paste(get_hierarchy_var_prefix(a), "MAX_IC", sep="")
    get(name, env)
}
get_hierarchy_max_depth<-function(env){
    get("_mesh_hierarchy_ANY_MAX_DEPTH", envir=env)
}
#read mesh tree and gather info about it (MAX, WMAX, etc...), store into env
read_and_parse_mtrees<-function(path, env){
    tmp<-scan(path, "", sep="\n")
    tmp<-strsplit(tmp, ";")
    tree_node_list<-list()
    tree_node_list_len<-c()
    parse_record<-function(i){
        tname<-get_term_var_name(i[[1]])
        if(!exists(tname, envir=env))
            assign(tname, new_term(des=i[[1]]), envir=env)
        tmp2<-get(tname, envir=env)
        tree_path<-get_tree_path(i[[2]])
        tmp2$nodes<-c(tmp2$nodes, list(tree_path))
        assign(tname, tmp2, envir=env)
        tmph<-tree_path[[1]]
        assign(get_node_var_name(i[[2]]), new_node(n=tree_path, des=i[[1]], h=tmph), envir=env)
    }
    lapply(hierarchies,
        function(a){
            assign(get_node_var_name(a), new_node(n=a, des=a, h=a), envir=env)
            assign(get_term_var_name(a), new_term(des=a, ns=list(a)), envir=env)
        })
    parse_node<-function(a){
        get_tree_path(a[[2]])
    }
    lapply(tmp, parse_record)
    tree_node_list<-c(sapply(tmp, parse_node), hierarchies)
    tree_node_list_len<-sapply(tree_node_list, length)
    tree_node_list<-tree_node_list[sort(tree_node_list_len, decreasing=TRUE, index.return=TRUE)$ix]
    tree_node_list_len<-NULL
    for(i in tree_node_list){
        update_node_max_length(i, env)
        if(length(i) == 1){
            node<-get(get_node_var_name_v(i), envir=env)
            update_hierarchy_max_length(i[[1]], node$max, node$wmax, env)
        }
    }
}

#read information content, ignore nodes not exist in env
read_information_content<-function(path, env){
    term_tmp<-scan(path[[1]], "", sep="\n")
    term_tmp<-strsplit(term_tmp, ";", fixed=TRUE)
    record_term_ic<-function(i){
        if(!exists(get_term_var_name(i[[1]]), envir=env)){
            print(c("error term not found?", get_term_var_name(i[[1]])))
            return(NULL)
        }
        tmp<-get(get_term_var_name(i[[1]]), envir=env)
        tmp$ic<-as.numeric(i[[2]])
        assign(get_term_var_name(i[[1]]), tmp, envir=env)
    }
    lapply(term_tmp, record_term_ic)
    node_tmp<-scan(path[[2]], "", sep="\n")
    node_tmp<-strsplit(node_tmp, ";", fixed=TRUE)
    record_node_ic<-function(i){
        if(!exists(get_node_var_name(i[[1]]), envir=env)){
            print(c("error node not found?", get_node_var_name(i[[1]])))
            return(NULL)
        }
        tmp<-get(get_node_var_name(i[[1]]),envir=env)
        tmp$ic<-as.numeric(i[[2]])
        #Update max information content as used in method Resnik
        h<-get_tree_path(i[[1]])[[1]]
        icname<-paste(get_hierarchy_var_prefix(h), "MAX_IC", sep="")
        if(!exists(icname, envir=env))
            assign(icname, tmp$ic, envir=env)
        else
            if(tmp$ic>get(icname, envir=env))
                assign(icname, tmp$ic, envir=env)

        assign(get_node_var_name(i[[1]]), tmp, envir=env)
    }
    lapply(node_tmp, record_node_ic)
    return(NULL)
}
`_update_db`<-function(mtree="mtrees2011.bin", terminfo="termInfo_2011.txt", nodeinfo="nodeInfo_2011.txt"){
    y<-new.env()
    read_and_parse_mtrees(mtree, y)
    read_information_content(path=c(terminfo, nodeinfo), env=y)
    return(y)
}
update_db<-function(fn="mesh.rda", mtree="mtrees2011.bin", terminfo="termInfo_2011.txt", nodeinfo="nodeInfo_2011.txt"){
    MeshSimData<-`_update_db`(mtree, terminfo, nodeinfo)
    save(MeshSimData, file=fn, compress="xz")
}
get_node_information_content<-function(a, f, env){
    if(f == 'heading'){
        tmp<-get(get_node_var_name_v(a), envir=env)
        tmp<-get(get_term_var_name(tmp$descriptor), envir=env)
        return(tmp$ic)
    }
    tmp<-get(get_node_var_name_v(a), envir=env)
    return(tmp$ic)
}
print_node<-function(a, env){
    for(i in 2:length(a)){
        tmp<-get(get_node_var_name_v(a[1:i]), envir=env)
        print(paste(tmp$descriptor, " [", paste(a[2:i], collapse="."), "] "))
    }
    return(NULL)
}
print_heading<-function(ap, env){
    a<-get_term_var_name(ap)
    a<-get(a, envir=env)
    a<-a$nodes
    tmp<-function(a){
        print_node(a, env)
    }
    lapply(a, tmp)
    return(NULL)
}
check_and_set_data<-function(a=NULL){
    if(is.null(a)){
        if(exists("MeshSimData")){
            return(get("MeshSimData"))
        }else{
            data("mesh")
            if(exists("MeshSimData"))
                return(get("MeshSimData"))
        }
    }else
        return(a)
}
get_node_name<-function(a, len){
    if(len<=1)return(a[[1]])
    paste(a[2:len], collapse=".")
}
`_build_node_tree`<-function(a, env){
    tmp<-strsplit(a, ".", fixed=TRUE)[[1]]
    nodes<-ls(env, pattern=get_node_var_name(tmp[[1]]))
    tmpe<-new.env()
    for(i in 1:length(nodes)){
        ttmp <- get(nodes[[i]], envir=env)$node;
        p<-paste("_build_node_tree_",get_node_name(ttmp, length(ttmp)-1))
        qname<-get_node_name(ttmp, length(ttmp))
        q<-paste("_build_node_tree_",qname)
        u<-list(child=list());
        if(exists(p, envir=tmpe))
            u<-get(p, envir=tmpe);
        u$child[[qname]] <- qname;
        assign(p, envir=tmpe, value=u);
        if(!exists(q, envir=tmpe)){
            u<-list(child=list(),term=get(nodes[[i]], envir=env)$descriptor);
            assign(q, envir=tmpe, value=u);
        }
    }
    np <- get(paste("_build_node_tree_",a), tmpe)$child
    np$term = get(paste("_build_node_tree_",a), tmpe)$term
    oname <- a
    for(i in (length(tmp)-1):1){
        noname <- paste(tmp[1:i], collapse=".")
        name <- paste("_build_node_tree_", noname)
        np2 <- get(name, tmpe)$child
        np2$term <- get(name, tmpe)$term
        np2[[oname]] <- np
        np <- np2
        oname <-noname
    }
    np2 <- list()
    np2[[tmp[[1]]]]<-np
    return(np2)
}
