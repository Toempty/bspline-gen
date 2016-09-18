A <- commandArgs( trailingOnly=T );
dat <- read.table( A[[1]], header=F );
spline.dat <- split( dat[,2:3], dat[,1]);
png( A[[2]] );
plot( range(dat[,2]), range(dat[,3]), xlab="X", ylab="Y" );
for( i in 1:length(spline.dat) ) lines( spline.dat[[i]], col=rainbow(length(spline.dat))[i], lwd=2 )
dev.off()
