sigGeneSet <-
function(setp, cutoff = 0.1, dualSig = (0:2)[2], 
    qpval = c("q.BH", "P.erlang")[1]) {
    if (is.list(setp)) {
        gs.name = rownames(setp$greater)
        setp$less = setp$less[gs.name, ]
        greater = setp$greater[, "P.erlang"] < setp$less[, "P.erlang"]
        sel.greater.0 = !is.na(setp$greater[, qpval]) & setp$greater[, 
            qpval] < cutoff
        sel.less.0 = !is.na(setp$less[, qpval]) & setp$less[, 
            qpval] < cutoff
        if (dualSig == 0) {
            sel.greater = sel.greater.0 & !sel.less.0
            sel.less = sel.less.0 & !sel.greater.0
        }
        else if (dualSig == 1) {
            sel.greater = sel.greater.0 & greater
            sel.less = sel.less.0 & !greater
        }
        else if (dualSig == 2) {
            sel.greater = sel.greater.0
            sel.less = sel.less.0
        }
        else {
            print("incorrect value for dualSig arguement")
        }
        ord = order(setp$less[sel.less, "P.erlang"])
        setp.sig = rbind(setp$greater[sel.greater, ], rbind(setp$less[sel.less, 
            ])[ord, ])
        rownames(setp.sig) = c(gs.name[sel.greater], gs.name[sel.less][ord])
        print(paste("there are", sum(sel.greater), "signficantly up-regulated gene sets"))
        print(paste("there are", sum(sel.less), "signficantly down-regulated gene sets"))
        return(setp.sig)
    }
    else {
        gs.name = rownames(setp)
        sel = !is.na(setp[, qpval]) & setp[, qpval] < cutoff
        setp.sig = rbind(setp[sel, ])
        rownames(setp.sig) = gs.name[sel]
        print(paste("there are", sum(sel), "signficantly two-direction perturbed gene sets"))
        return(setp.sig)
    }
}

