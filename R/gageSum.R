gageSum <-
function(p.results, ref = NULL, mstat = NULL, 
    setsizes = NULL, same.dir = TRUE, compare = "paired", use.fold = TRUE, 
    weights = NULL, full.table = FALSE, ...) {
    nc <- ncol(p.results)
    if (is.null(weights)) 
        sg.glob <- apply(-log(p.results), 1, sum)
    else sg.glob <- apply(-log(p.results), 1, function(x) sum(x * 
        weights))
    if (compare == "unpaired" & use.fold) {
        p.erlang <- pgamma(sg.glob/length(ref), shape = nc/length(ref), 
            rate = 1, lower.tail = FALSE)
    }
    else if (!is.null(weights)) {
        p.erlang <- pgamma(sg.glob, shape = sum(weights), rate = 1, 
            lower.tail = FALSE)
    }
    else {
        p.erlang <- pgamma(sg.glob, shape = nc, rate = 1, lower.tail = FALSE)
    }
    
    
    Fdr = NULL
    
    if (require(fdrtool)) {
        q.fdrtool = p.erlang
        q.fdrtool[!is.na(p.erlang)] = fdrtool(p.erlang[!is.na(p.erlang)], 
            statistic = "pvalue", plot = FALSE, verbose = options()$verbose)$qval
        Fdr = cbind(Fdr, q.fdrtool)
    }
    else print("package fdrtool not available")
    
    if (require(multtest)) {
        q.BH = p.erlang
        res <- mt.rawp2adjp(p.erlang[!is.na(p.erlang)], "BH")
        q.BH[!is.na(p.erlang)] = res$adjp[order(res$index), "BH"]
        Fdr = cbind(Fdr, q.BH)
    }
    else print("package multtest not available")
    
    p.glob = cbind(exp(-sg.glob/ifelse(is.null(weights), nc, 
        sum(weights))), mstat, p.erlang, Fdr, setsizes)
    colnames(p.glob) = c("P.geomean", "stat.mean", "P.erlang", 
        colnames(Fdr), "set.size")
    return(p.glob)
}

