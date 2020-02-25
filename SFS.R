
setwd("C:\\Users\\Anoob\\Desktop\\SFS_EcoGen")

list.files()

SFS <- scan("XGL_outFold.sfs")

sumSFS <- sum(SFS)

pctPoly <- 100*(1-(SFS[1]/sumSFS))

plotSFS <- SFS[-c(1,length(SFS))]

barplot(plotSFS)

div <- read.table("XGL_folded_allsites.thetas.idx.pestPG")

colnames(div) <- c("window", "chrname","wincenter","tW","tP","tF","tH","tL",
                   "tajD","fulif","fuliD","fayH","zengsE","numSites")

div$tWpersite <- div$tW/div$numSites
div$tPpersite <- div$tP/div$numSites

par(mfrow=c(2,2))
hist(div$tWpersite, col="gray",xlab="Theta-W",main="")
hist(div$tPpersite, col="gray",xlab="Theta-P",main="")
hist(div$tajD, col="gray",xlab="Tajima's D",main="")
barplot(plotSFS, xlab="SFS",)
dev.off() # closes the pdf and writes it to file

summary(div)
