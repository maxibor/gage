require(gage)
print(0)
### gage 1-d #implement weights later
gage.res=gage(exprs=exprs, gsets=gsets, ref = ref, samp = samp, 
  set.size = set.size, same.dir = TRUE, compare = compare,
  rank.test = rank.test, use.fold = use.fold, saaTest = gs.test)
write.table(rbind(gage.res$greater, gage.res$less), 
            file = "gage.res.txt", sep = "\t")

print(1)

### significant.genesets
gage.res.sig<-sigGeneSet(gage.res, outname="gage.res", cutoff=cutoff)
sig.gs=unique(c(rownames(gage.res.sig$greater), rownames(gage.res.sig$less)))
nsig=length(sig.gs)
if(nsig>0) {
  write.table(rbind(gage.res.sig$greater, gage.res.sig$less), 
              file = "gage.res.sig.txt", sep = "\t")
} else print("No gene set selected in 1d-test!")

print(2)

### gage 2-d
if(test.2d & gs.type!="go"){
  gage.res.2d=gage(exprs=exprs, gsets=gsets, ref = ref, samp = samp, 
    set.size = set.size, same.dir = FALSE, compare = compare,
    rank.test = rank.test, use.fold = use.fold, saaTest = gs.test)
  write.table(gage.res.2d$greater, file = "gage.res.2d.txt", 
              sep = "\t")
  gage.res.2d.sig<-sigGeneSet(gage.res.2d, outname="gage.res", cutoff=cutoff)
  sig.gs.2d=rownames(gage.res.2d.sig$greater)
  nsig.2d=length(sig.gs.2d)
  if(nsig.2d>0) {
    write.table(gage.res.2d.sig$greater, file = "gage.res.2d.sig.txt", 
                sep = "\t")
  } else print("No gene set selected in 2d-test!")
} else sig.gs.2d=NULL
  sig.gs.all=unique(c(sig.gs,sig.gs.2d))
  nsig.all=length(sig.gs.all)

print(3)


### output
  if(nsig.all>0){

### geneData
    if(gs.type!="user") outnames =sapply(strsplit(sig.gs.all, " "), "[", 1)
    else outnames=sig.gs.all
    outnames = gsub(" |:|/", "_", outnames)
    for (i in (1:nsig.all)[1:3]) {
      geneData(genes = gsets[[sig.gs.all[i]]], exprs = exprs, ref = ref, 
               samp = samp, outname = outnames[i], txt = T, heatmap = T, 
               Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
    }

print(4)


### pathview
    if(do.pathview & gs.type=="kegg"){
      require(pathview)
      if(compare=="paired") exprs.d=exprs[,samp]-exprs[,ref]
      else exprs.d=exprs[,samp]-rowMeans(exprs[,ref])
      if(data.type=="gene")
        pv.out.list <- sapply(outnames[1:3], function(pid) pathview(gene.data = exprs.d, pathway.id = pid, species = species))
      else pv.out.list <- sapply(outnames[1:3], function(pid) pathview(cpd.data = exprs.d, pathway.id = pid, species = species))
    }
  } else print("No gene set selected by GAGE, you may relax the cutoff q-value!")

print(5)
  save.image("workenv.RData")
