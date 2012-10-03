library(bio3d)
library(cluster)

plot.pca2 <-
function(x, pch=16, col=par("col"), cex=0.8, mar=c(4, 4, 1, 1), ...) {

  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))

  par(mfrow=c(2, 2), cex=cex, mar=mar)
  par(pty="s")
  plot(x$z[,1],x$z[,2], type="p", pch=pch, xlab="PC1", ylab="PC2", col=col, ...)
  abline(h=0,col="gray",lty=2);abline(v=0,col="gray",lty=2)
  plot(x$z[,3],x$z[,2], type="p", pch=pch, xlab="PC3", ylab="PC2", col=col,...)
  abline(h=0,col="gray",lty=2);abline(v=0,col="gray",lty=2)
  plot(x$z[,1],x$z[,3], type="p", pch=pch, xlab="PC1", ylab="PC3", col=col,... )
  abline(h=0,col="gray",lty=2);abline(v=0,col="gray",lty=2)
  plot.pca.scree(x$L,ylim=c(0,30))
}



##-- Read a trajectory file
trj <- read.ncdf("prod_CA.nc")
#trj <- trj[1:1000,]  ### FOR TESTING!
##-- Read a PDB file with the CA from the minimized structure
pdb <- read.pdb("prod_CA.pdb.1")
ca_inds <- atom.select(pdb, "calpha")

# anything "2" is the crystal structure
pdb2 <- read.pdb("prot.pdb")
ca_inds2 <- atom.select(pdb2, "calpha")$atom
ca_seq <- aa321(pdb2$atom[ca_inds2,'resid'])
ca_resno <- pdb2$atom[ca_inds2,'resno']
ca_names <- paste(ca_seq,ca_resno) # x-axis labels
ca_chain <- pdb2$atom[ca_inds2,'chain']
ca_len <- length(ca_inds2)
ca_breaks <- list()
for (i in 2:ca_len) { # get breaks in the chain and plot vertical line
	if ((as.integer(ca_resno[i]) - as.integer(ca_resno[i-1])) != 1 || ca_chain[i] != ca_chain[i-1]) {
		ca_breaks <- c(ca_breaks,i-0.5)
	}
}




##-- Secondary structure assignment with DSSP/usr/local/lib/vmd
sse <- dssp(pdb2, exepath = "/opt/local/bin/", resno=FALSE)
#sse <- stride(pdb2, exepath = "/usr/local/lib/vmd")
# +/- 0.5 to extend SS rectangles
if ( ! is.null(sse$helix$start)) {
	sse$helix$start <- sse$helix$start - 0.5
	sse$helix$end <- sse$helix$end + 0.5
}
if ( ! is.null(sse$sheet$start)) {
	sse$sheet$start <- sse$sheet$start - 0.5
	sse$sheet$end <- sse$sheet$end + 0.5
}

# Fit trajectory onto pdb with CA
xyz <- fit.xyz(pdb$xyz, trj, fixed.inds=ca_inds$xyz, mobile.inds=ca_inds$xyz)

##-- RMSD
ave <- apply(xyz,2,mean)
r_ave <- rmsd(a=ave, b=xyz)
r_first <- rmsd(a=xyz[1,], b=xyz)
r_crystal <- rmsd(a=c(as.numeric(t(pdb2$atom[ca_inds2,c('x','y','z')]))), b=xyz, fit=TRUE)
png("rmsd.png", pointsize=20, width=1200, height=1200)
par(mfrow=c(3,1),cex=1.2)
hist(r_first,xlim=c(0,2.5))
#hist(r_first)
time <- seq(1,dim(trj)[1]*2,2)
plot.bio3d(time,r_first, typ="l", ylab="RMSD from first frame",xlab="time (ps)",ylim=c(0,2.5))
plot.bio3d(time,r_ave, typ="l", ylab="RMSD fr ave/crystal struct",xlab="time (ps)",ylim=c(0,2.5))
par(new=TRUE)
plot.bio3d(time,r_crystal, typ="l", ylab="",xlab="",ylim=c(0,2.5),col='blue')
dev.off()

##-- RMSF

fluc=rmsf(xyz)
#fluc <- rmsf(trj)
#write.pdb(pdb, b=fluc, file="fluct.pdb")

png("rmsf.png", pointsize=20, width=1200, height=1200)

par(mar=c(5,4,4,4) + 0.1, cex=1.0)
ca_gapped <- seq(1,ca_len,5)

plot.bio3d(fluc, sse=sse, typ="l", ylim=c(0,3), ylab="RMSF (A)", xlab="Residues",helix.col = "red", sheet.col = "yellow", sse.border = FALSE, axes=FALSE)
axis(2)
axis(1,at=ca_gapped,labels=ca_names[ca_gapped],las=2)
axis(3,at=ca_gapped,labels=ca_names[ca_gapped],las=2)
abline(v=ca_breaks)

par(new=TRUE)
plot(pdb2$atom[pdb2$calpha,"b"],type='l',col='blue',xlab='',ylab='',xaxt='n',yaxt='n')
axis(4, labels=FALSE)
at = axTicks(4)
mtext(side = 4, text = at, at = at, col = "blue", line = 1) 
mtext('B-Factor',side=4,line=3,col='blue')
box()

dev.off()

##-- DCCM

cij <- dccm(xyz)
png("dccm2.png", pointsize=24, width=1200, height=1200)
library(lattice)
x <- y <- c(1:ca_len)

ca_gapped <- seq(1,ca_len,5)
filled.contour(x, y, cij,color.palette=bwr.colors,xlab="Residue Number", ylab="Residue Number",main="DCCM: dynamic cross-correlation map",zlim=c(-1,1),plot.axes={axis(1,at=ca_gapped,labels=ca_names[ca_gapped],las=2); axis(2,at=ca_gapped,labels=ca_names[ca_gapped],las=2); abline(v=ca_breaks,h=ca_breaks); box()})

dev.off()

png("dccm.png", pointsize=24, width=1200, height=1200)
library(lattice)
x <- y <- c(1:ca_len)

ca_gapped <- seq(1,ca_len,5)
levels=c(-1, -0.7, -0.4, 0.4, 0.7, 1)
colors=c('darkblue','skyblue','white','pink','red2')
filled.contour(x, y, cij,col=colors,levels=levels,xlab="Residue Number", ylab="Residue Number",main="DCCM: dynamic cross-correlation map",zlim=c(-1,1),plot.axes={axis(1,at=ca_gapped,labels=ca_names[ca_gapped],las=2); axis(2,at=ca_gapped,labels=ca_names[ca_gapped],las=2); abline(v=ca_breaks,h=ca_breaks); box()})

dev.off()

##-- Do PCA
trj.pca<- pca.xyz(xyz)
png("pca.png", pointsize=8, width=1200, height=1200)
plot.pca2(trj.pca,col=rainbow(length(trj.pca$z[,1])), cex=3.0, xlim=c(-25,25), ylim=c(-25,25))

dev.off()

## Plot atom wise loadings
# it is included within pca.loadings
png("pca_rmsf.png", pointsize=24, width=1200, height=1200)
par(mfrow=c(3,1),mar=c(3,4,1,1))
plot.bio3d(trj.pca$au[, 1], typ="l", sse=sse, ylim=c(0,0.30),ylab="PC1 (A)",helix.col = "red", sheet.col = "yellow", sse.border = TRUE)
plot.bio3d(trj.pca$au[, 2], typ="l", sse=sse, ylim=c(0,0.30),ylab="PC2 (A)",helix.col = "red", sheet.col = "yellow", sse.border = TRUE)
plot.bio3d(trj.pca$au[, 3], typ="l", sse=sse, ylim=c(0,0.30),ylab="PC3 (A)",helix.col = "red", sheet.col = "yellow", sse.border = TRUE)
dev.off()

# Conformer Clustering in PC space. Clustering of structures from the trajectory in the PC1-3 planes
#hc <- hclust(dist(trj.pca$z[, 1:3]))
#grps <- cutree(hc, k=3)
#png("pc_clustering.png",res=300, pointsize=2, width=1200, height=1200)
#plot.pca2(trj.pca, col = c("red", "black", "blue")[grps], cex=3.0, xlim=c(-10,10), ylim=c(-10,10))
###plot(hbclust(dist(trj.pca$z[, 1:3])))
#dev.off()

# Write PC trajectory
#a <- mktrj.pca(trj.pca, pc=1, file="pc1.pdb", resno = ca_resno, resid = ca_seq)
#b <- mktrj.pca(trj.pca, pc=2, file="pc2.pdb", resno = ca_resno, resid = ca_seq)
#c <- mktrj.pca(trj.pca, pc=3, file="pc3.pdb", resno = ca_resno, resid = ca_seq)


