
# run ordisurf on the NMDS plots
all.equal(rownames(myMDS$points)[-285],allmeta0922$id)
surf <- vegan::ordisurf(x = myMDS, y = c(allmeta0922$bio2, NA))
test = summary(surf)
