source('geneset.R')

sets = loadBroadSets("data/msigdb_v3.0.xml")

setsEntrez = getGmt("data/msigdb.v3.0.entrez.gmt", geneIdType = EntrezIdentifier())

library("hom.Hs.inp.db")
library("org.Mm.eg.db")
library("org.Hs.eg.db")

## Map to mouse ids
setsMapped = mapSetsToSpecies(sets, setsEntrez = setsEntrez, inProt = org.Hs.egENSEMBLPROT, outProt = org.Mm.egENSEMBLPROT2EG, hom = hom.Hs.inpMUSMU)

saveSetCache(setsMapped, "data/msigdb.v3.0.mouse.c2c3c5.noCGP.RData")