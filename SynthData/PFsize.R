## Calculate size of the place fields in synthetic data
## Figure 3g
## Requires generation of the synthetic datasets (sim_thetaseq_kalman.R) and preprocessing them (preproc_syntetic_data.R as called in Decode_theta_sim.R)
dir.create('./SimFigs/01PlaceFieldSize/', showW=F)


library(gplots)
Ncells <- 200
codes <- c('MAP', 'sampling', 'PPC_var', 'DDC')
PFsize <- array(NA, dim=c(4, 3, 200))

for(i_code in 1:4){
	mapsname <- paste('./3maps_Tmax1800_Ncells200_', codes[i_code], '_spPth32irreg_0.RData', sep='')
	load(mapsname)
	for (i_cell in 1:200){
		emap <- maps$ratemaps_early[i_cell,,]
		mmap <- maps$ratemaps_mid[i_cell,,]
		lmap <- maps$ratemaps_late[i_cell,,]
		PFsize[i_code, 1, i_cell] <- sum(emap > 0.1*max(emap)) / 1600 * 4
		PFsize[i_code, 2, i_cell] <- sum(mmap > 0.1*max(mmap)) / 1600 * 4
		PFsize[i_code, 3, i_cell] <- sum(lmap > 0.1*max(lmap)) / 1600 * 4
	}
}

## late / early
mm <- apply(PFsize[,3,] / PFsize[,1,], 1, mean)
ss <- apply(PFsize[,3,] / PFsize[,1,], 1, sd)

## late / mid
# mm <- apply(PFsize[,3,] / PFsize[,2,], 1, mean)
# ss <- apply(PFsize[,3,] / PFsize[,2,], 1, sd)

# ## mid / early
# mm <- apply(PFsize[,2,] / PFsize[,1,], 1, mean)


pdf('./01PlaceFieldSize/PFsize_irreg_MAP.pdf', 6, 3, useD=F)
par(mfcol=c(1,2))
par(mar=c(6,6,4,4))
plot(PFsize[1,1,], PFsize[1,3,], pch=21, bg=grey(0.75), axes=F, xlim=c(0, 4), ylim=c(0, 4), xlab='early (m2)', ylab='late (m2)', main='MAP')
# plot(PFsize[2,1,], PFsize[2,3,], pch=21, bg=grey(0.75), axes=F, xlim=c(0, 4), ylim=c(0, 4), xlab='early (m2)', ylab='late (m2)')
abline(0, 1, lwd=2)
points(median(PFsize[2,1,]), median(PFsize[2,3,]), col=2, pch=3, cex=3, lwd=3)
axis(1)
axis(2, las=2)

par(mar=c(4,5,1,5))

barplot2(mm[c(3,4,2,1)], 0.8, plot.ci = TRUE, ci.l = mm[c(3,4,2,1)]+ss[c(3,4,2,1)], ci.u = mm[c(3,4,2,1)]-ss[c(3,4,2,1)], axes=F)
axis(1, 1:4-0.5, codes[c(3,4,2,1)], las=2, tick=F)
axis(2, las=2)
dev.off()

t.test(PFsize[1,3,]-PFsize[1,1,])
t.test(PFsize[2,3,]-PFsize[2,1,])
t.test(PFsize[3,3,]-PFsize[3,1,])
t.test(PFsize[4,3,]-PFsize[4,1,])


############################################################
## violin plots

library(ggplot2)
LpE <- NA
code <- NA
for (i_code in 1:4){
	LpEi <- PFsize[i_code,3,] / PFsize[i_code,1,]
	LpE <- c(LpE, LpEi[!is.nan(LpEi)])
	code <- c(code, rep(codes[i_code], length(LpEi)))
}
LpE <- LpE[-1]
code <- code[-1]

field_size <- data.frame(LpE)
field_size$codes <- as.factor(code)

FieldSize_plot_file = './01PlaceFieldSize/PFsize_rel_viol.pdf'
pdf(file= FieldSize_plot_file, 3, 3, useD=F)

p <-ggplot(field_size, aes(x=codes, y=LpE, col=codes), main='') + geom_violin()		
p + stat_summary(fun.data='mean_sdl', na.rm=T, fun.args = list(mult = 1), geom='pointrange', color=1) + theme_minimal() + ylim(0, 4.5)
dev.off()

