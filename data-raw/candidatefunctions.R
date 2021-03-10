



round(zoneDDR$exploitB[40:47,,1],2)

round(zoneDDR$Nt[,nyrs,1,1:5])



cbind(round(zoneDD$matureB[,1:4],3),round(zoneDDR$matureB[,1:4,5],3))











pLML <- projC$projLML[1] # needs development to allow variation in projLML
for (pop in 1:numpop) {
  selL50 <- popdefs["SelP1",pop]
  selL95 <- popdefs["SelP2",pop]
  projSel[,pop] <- logistic((pLML + selL50),selL95,midpts)
  projSelWt[,pop] <- projSel[,pop] * zoneC[[pop]]$WtL
}



















