`_nodeSim`<-function(a, b, method="SP", frame="node", env=NULL){
    if(strsplit(a[[1]], "")[[1]][[1]] != strsplit(b[[1]], "")[[1]][[1]]){
        return(0)
    }
    func<-paste("mesh_node_similarity_", method, sep="")
    func<-get(func)
    func(a, b, frame, env)
}

`_mnodeSim`<-function(a, b, method="SP", frame="node", env=NULL){
    run<-function(na){
        run2<-function(nb){
            `_nodeSim`(na, nb, method, frame, env)
        }
        return(sapply(b, run2))
    }
    res<-as.vector(sapply(a, run))
    dim(res)<-c(length(a), length(b))
    return(res)
}
nodeSim<-function(node1, node2, method="SP", frame="node", env=NULL){
    env=check_and_set_data(env);
    `_nodeSim`(get_tree_path(node1), get_tree_path(node2), method, frame, env)
}
mnodeSim<-function(nodeList1, nodeList2, method="SP", frame="node", env=NULL){
    run<-function(na){
        run2<-function(nb){
            nodeSim(na, nb, method, frame, env)
        }
        return(sapply(nodeList2, run2))
    }
    res<-as.vector(sapply(nodeList1, run))
    dim(res)<-c(length(nodeList1), length(nodeList2))
    return(res)
}
headingSim<-function(heading1, heading2, method="SP", frame="node", env=NULL){
    env=check_and_set_data(env);
    an<-get(get_term_var_name(heading1), envir=env)$nodes
    bn<-get(get_term_var_name(heading2), envir=env)$nodes
    res<-`_mnodeSim`(an, bn, method, frame, env)
    if(frame == "heading")
        return(max_method1(res))
    else
        return(max_method2(res))
}
mheadingSim<-function(headingList1, headingList2, method="SP", frame="node", env=NULL){
    run<-function(na){
        run2<-function(nb){
            headingSim(na, nb, method, frame, env)
        }
        return(sapply(headingList2, run2))
    }
    res<-as.vector(sapply(headingList1, run))
    dim(res)<-c(length(headingList1), length(headingList2))
    return(res)
}
nodeInfo<-function(node, env=NULL){
    env=check_and_set_data(env);
    `_build_node_tree`(node, env);
}
termInfo<-function(term, env=NULL){
    env=check_and_set_data(env);
    run<-function(na){
        `_build_node_tree`(paste(na[2:length(na)],collapse="."), env);
    }
    an<-get(get_term_var_name(term), envir=env)$nodes
    lapply(an, run);
}
`_docSim`<-function(a, b, method="SP", frame="node", env=NULL){
    res<-mheadingSim(a, b, method, frame, env)
    return(max_method2(res))
}
headingSetSim<-function(headingSet1, headingSet2, method="SP", frame="node", env=NULL){
    a<-headingSet1[headingSet1!=""]
    b<-headingSet2[headingSet2!=""]
    res<-mheadingSim(a, b, method, frame, env)
    return(max_method2(res))
}

