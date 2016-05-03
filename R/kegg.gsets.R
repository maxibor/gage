kegg.gsets<-function(species="hsa", id.type="kegg"){
  use.egid=FALSE
  if(species!="ko") {
    species.info=kegg.species.code(species, na.rm=T, code.only=F)
    species=species.info[1]
    use.egid=species.info[4]=="1"
  } else if(tolower(id.type)!="kegg"){
    id.type="kegg"
  }
  path.list <- keggLink("pathway",species)
  genes=gsub(paste(species,":",sep=""), "", names(path.list))
  if(tolower(id.type)=="entrez" & !use.egid){
    gid.map=keggConv("ncbi-geneid",species)
    kegg.ids=gsub(paste(species, ":", sep=""), "", names(gid.map))
    eg.ids=gsub("ncbi-geneid:", "", gid.map)
    id.idx=match(genes, kegg.ids)
    genes[!is.na(id.idx)]=eg.ids[id.idx[!is.na(id.idx)]]
  }
  paths=gsub(paste("path:",species,sep=""), "", path.list)
  kg.sets=split(genes, paths)

  data(khier, package="gage")
  kh.idx=match(names(kg.sets), substr(khier[,3], 1,5))
  kg.sets=kg.sets[!is.na(kh.idx)]
  kh.idx=kh.idx[!is.na(kh.idx)]
  ks.names=khier[kh.idx,3]

  aterms=unique(khier[,1])
  signals=khier[which(khier[,1] %in% aterms[2:5]),3]
  metabs=khier[which(khier[,1] == aterms[1]),3]
  diseases=khier[which(khier[,1] == aterms[6]),3]

  sig.idx=which(ks.names %in% signals)
  met.idx=which(ks.names %in% metabs)
  dise.idx=which(ks.names %in% diseases)
  sigmet.idx=c(sig.idx,met.idx)
  names(kg.sets)=paste(species, ks.names, sep="")
  res=list(kg.sets=kg.sets, sigmet.idx=sigmet.idx, sig.idx
    =sig.idx, met.idx=met.idx, dise.idx=dise.idx)
  return(res)
}
