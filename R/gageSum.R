gageSum <- function(rawRes, ref = NULL, test4up = TRUE, 
    same.dir = TRUE, compare = "paired", use.fold = TRUE, weights = NULL, 
    full.table = FALSE, use.stouffer = TRUE, ...) {
    if (test4up) 
        p.results = rawRes$p.results
    else p.results = rawRes$ps.results
    results = rawRes$results
    cns = unique(colnames(results))
    ord = rank(cns)
    nc <- ncol(p.results)
    if (is.null(weights)) 
        sg.glob <- apply(-log(p.results), 1, sum)
    else sg.glob <- apply(-log(p.results), 1, function(x) sum(x * 
        weights))
    if (use.stouffer) {
        q.results = qnorm(p.results)
        if (is.null(weights)) 
            q.glob <- apply(q.results, 1, sum)
        else q.glob <- apply(q.results, 1, function(x) sum(x * 
            weights))
        if (compare == "unpaired" & use.fold & length(ref)>0) 
            mod = (nc * length(ref))^(-0.5)
        else if (!is.null(weights)) 
            mod = sum(weights)^(-0.5)
        else mod = nc^(-0.5)
        p.val <- pnorm(q.glob * mod)
        if(length(cns)!=nc & length(cns)>1){
          p.results = as.matrix(t(aggregate(t(q.results), list(colnames(q.results)), 
            function(x) pnorm(mean(x)))[, -1]))[, ord]
          colnames(p.results) = cns
        }
      }
    else {
        if (compare == "unpaired" & use.fold & length(ref)>0) {
            p.val <- pgamma(sg.glob/length(ref), shape = nc/length(ref), 
                rate = 1, lower.tail = FALSE)
        }
        else if (!is.null(weights)) {
            p.val <- pgamma(sg.glob, shape = sum(weights), rate = 1, 
                lower.tail = FALSE)
        }
        else {
            p.val <- pgamma(sg.glob, shape = nc, rate = 1, lower.tail = FALSE)
        }
        if(length(cns)!=nc & length(cns)>1){
          p.results = as.matrix(t(aggregate(t(p.results), list(colnames(p.results)), 
            function(x) exp(mean(log(x))))[, -1]))[, ord]
          colnames(p.results) = cns
        }
      }
        if(length(cns)!=nc & length(cns)>1){
    results = as.matrix(t(aggregate(t(results), list(colnames(results)), 
        mean, na.rm = TRUE)[, -1]))[, ord]
    colnames(results) = cns
}
    mstat = apply(results, 1, mean, na.rm = TRUE)
    Fdr = NULL
    
    if (require(multtest)) {
        q.val = p.val
        res <- mt.rawp2adjp(p.val[!is.na(p.val)], "BH")
        q.val[!is.na(p.val)] = res$adjp[order(res$index), "BH"]
        Fdr = cbind(Fdr, q.val)
    }
    else print("package multtest not available")
    
    p.glob = cbind(exp(-sg.glob/ifelse(is.null(weights), nc, 
        sum(weights))), mstat, p.val, Fdr, rawRes$setsizes, p.results)
    colnames(p.glob) = c("p.geomean", "stat.mean", "p.val", colnames(Fdr), 
        "set.size", colnames(p.results))
    results = cbind(mstat, results)
    colnames(results) = c("stat.mean", cns)
    ind <- order(p.glob[, "p.val"])
    return(list(p.glob = p.glob[ind, ], results = results[ind, 
        ]))
}
 
