

filename <- filenametopath(resdir,"regiontable.html")



capt <- "Table Regional Population Biological Properties."
filename <- filenametopath(resdir,"regiontable.html")
htmltable(results,filename,capt)



cat('<caption> </caption> \n',file=filename,append=TRUE)
cat('<tr> \n',file=filename,append=TRUE)
for (cl in 1:numcol) {
  cat('<th>',columns[cl],'</th> \n',file=filename,append=TRUE)
}
for (rw in 1:numrow) {
  cat('<tr> \n',file=filename,append=TRUE)
  for (cl in 1:numcol) { cat('<th>',res[rw,cl],'</th> \n',file=filename,append=TRUE) }
  cat('</tr> \n',file=filename,append=TRUE)
}
cat('</table>',file=filename,append=TRUE)


tmp <- NULL
for (cl in 1:numcol) tmp <- paste(tmp,paste('<th>',columns[cl],'</th>',collapse=""),collapse="")
tmp
