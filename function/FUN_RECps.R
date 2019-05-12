RECps <- function(x, ps) {
  TaxTab2 <- as.data.frame(x)
  
  list.g <- as.character(TaxTab2$Genus)
  list.f <- as.character(TaxTab2$Family)
  list.o <- as.character(TaxTab2$Order)
  list.c <- as.character(TaxTab2$Class)
  list.p <- as.character(TaxTab2$Phylum)
  list.k <- as.character(TaxTab2$Kingdom)
  list.REC <- character(length(list.g))
  
  for(i in 1:dim(TaxTab2)[1]){
    G = which(TaxTab2$Genus[i] == "" | TaxTab2$Genus[i] == "None")
    F = which(TaxTab2$Family[i] == "" | TaxTab2$Family[i] == "None")
    O = which(TaxTab2$Order[i] == "" | TaxTab2$Order[i] == "None")
    C = which(TaxTab2$Class[i] == "" | TaxTab2$Class[i] == "None")
    P = which(TaxTab2$Phylum[i] == "" | TaxTab2$Phylum[i] == "None")
    K = which(TaxTab2$Kingdom[i] == "" | TaxTab2$Kingdom[i] == "None")
    if(length(G) == 0){
      list.REC[i] <- list.g[i]
    } else if(length(F) == 0){
      list.REC[i] <- list.f[i]
    } else if(length(O) == 0){
      list.REC[i] <- list.o[i]
    } else if(length(C) == 0){
      list.REC[i] <- list.c[i]
    } else if(length(P) == 0){
      list.REC[i] <- list.p[i]
    } else if(length(K) == 0){
      list.REC[i] <- list.k[i]
    } else
      list.REC[i] <- "meow"
  }
  
  TaxTab2$REC <- list.REC
  TaxTab2$REC <- factor(TaxTab2$REC)
  merge_phyloseq(ps, TaxTab2 %>% as.matrix() %>% tax_table())
}