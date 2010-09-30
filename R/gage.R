gage <-
function(exprs, gsets, ref = NULL, samp = NULL, 
    set.size = c(10, 500), same.dir = TRUE, compare = "paired", 
    rank.test = FALSE, use.fold = TRUE, FDR.adj = TRUE, weights = NULL, 
    full.table = FALSE, saaPrep = gagePrep, saaTest = gs.tTest, 
    saaSum = gageSum, ...) {
    exprs = saaPrep(exprs, ref = ref, samp = samp, same.dir = same.dir, 
        compare = compare, rank.test = rank.test, use.fold = use.fold, 
        weights = weights)
    rawRes = saaTest(exprs, gsets = gsets, set.size = set.size, 
        same.dir = same.dir)
    p.results = rawRes$p.results
    mstat = rawRes$mstat
    setsize = rawRes$setsize
    p.glob = saaSum(p.results, ref = ref, mstat = mstat, setsize = setsize, 
        same.dir = same.dir, compare = compare, use.fold = use.fold, 
        weights = weights, full.table = full.table)
    ind <- order(p.glob[, 1])
    nc <- ncol(exprs)
    if (same.dir) {
        ps.results = rawRes$ps.results
        ps.glob = saaSum(ps.results, ref = ref, mstat = mstat, 
            setsize = setsize, same.dir = same.dir, compare = compare, 
            use.fold = use.fold, weights = weights, full.table = full.table)
        inds <- order(ps.glob[, 1])
        return(list(greater = (if (nc <= 10 | full.table) cbind(p.glob, 
            p.results) else p.glob)[ind, ], less = (if (nc <= 
            10 | full.table) cbind(ps.glob, ps.results) else ps.glob)[inds, 
            ]))
    }
    else return((if (nc <= 10 | full.table) cbind(p.glob, p.results) else p.glob)[ind, 
        ])
}

