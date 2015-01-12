lambda<-3
lca<-function(a, b){
    i<-1
    while((i <= min(length(a), length(b))) && (a[i] == b[i])){
        i<-i+1
    }
    i-1
}
path_length<-function(a, b){
    length(a)+length(b)-2*lca(a,b)
}
path_wlength<-function(a, b){
    tmp<-lca(a,b)
    sum(1/seq(tmp, length(a)))+sum(1/seq(tmp, length(b)))-1/tmp
}
mesh_node_similarity_SP<-function(a, b, f, env){
    #return 0 if hierarchy not match
    ha<-strsplit(a[[1]], "", fixed=TRUE)[[1]][1]
    max<-get_hierarchy_max_length(ha, env)
    (max-path_length(a, b))/max
}
mesh_node_similarity_WL<-function(a, b, f, env){
    ha<-strsplit(a[[1]], "", fixed=TRUE)[[1]][1]
    max<-get_hierarchy_max_wlength(ha, env)
    (max-path_wlength(a, b))/max
}
mesh_node_similarity_WP<-function(a, b, f, env){
    (2*lca(a, b))/(length(a)+length(b))
}
mesh_node_similarity_LC<-function(a, b, f, env){
    maxd<-get_hierarchy_max_depth(env)
    1-log(path_length(a, b)+1)/log(maxd*2+1)
}
mesh_node_similarity_Li<-function(a, b, f, env){
    H<-lca(a, b)
    L<-path_length(a, b)
    exp(-0.2*L)*(exp(0.6*H)-exp(-0.6*H))/(exp(0.6*H)+exp(-0.6*H))
}
mesh_node_similarity_Lord<-function(a, b, f, env){
    lca_path<-a[1:lca(a,b)]
    ic<-get_node_information_content(a, f, env)
    1-exp(-ic)
}
mesh_node_similarity_Resnik<-function(a, b, f, env){
    hic<-get_node_information_content(a[1:lca(a,b)], f, env)
    maxc<-get_hierarchy_max_ic(strsplit(a[[1]], "")[[1]][[1]], env)
    hic/maxc
}
mesh_node_similarity_Lin<-function(a, b, f, env){
    hic<-get_node_information_content(a[1:lca(a,b)], f, env)
    ic1<-get_node_information_content(a, f, env)
    ic2<-get_node_information_content(b, f, env)
    hic*2/(ic1+ic2)
}
mesh_node_similarity_JC<-function(a, b, f, env){
    hic<-get_node_information_content(a[1:lca(a,b)], f, env)
    ic1<-get_node_information_content(a, f, env)
    ic2<-get_node_information_content(b, f, env)
    dist<-ic1+ic2-2*hic
    exp(-dist/lambda)
    
}
max_method1<-function(a){
    m1<-function(i){
        max(as.vector(a[i,]))
    }
    max(sapply(1:dim(a)[1], m1))
}
max_method2<-function(a){
    m1<-function(i){
        max(as.vector(a[i,]))
    }
    lx<-dim(a)[1]
    x<-sum(sapply(1:lx, m1))
    a<-t(a)
    ly<-dim(a)[1]
    y<-sum(sapply(1:ly, m1))
    (x+y)/(lx+ly)
}

