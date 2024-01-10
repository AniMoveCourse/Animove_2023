
## for this script to work, run first the script "Lesson5A_2023_Segmentation.R", until line 272 has to be run (i.e. NSD, FTP, BPMM, EmbC, BCPA, MRW)

## this is just one way to extract the information of the segmentation algorithms, there are probably more and smarter ones
## first I create a graph of the segmentation
## than I transfer this information on to the track to be able to plot it

############## plot segments in box #######################################
pdf("segmentationPlots.pdf", width=14,height=(4*8))
par(mfrow=c(8,1))
#######
# NSD #
#######
# finding break points in the NSD using the lavielle function
lvNSD <- lavielle(drop_units(SieritNSD), Lmin=9, Kmax=12, type = "var")
brksNSDl <-  findpath(lvNSD, 12, plotit = F)
brksNSD <- unlist(brksNSDl)
plot(Sierit1day$timestamp, SieritNSD, type="l", xlab="", ylab="NSD (Km²)", main="NSD")
points(Sierit1day$timestamp, SieritNSD,pch=20,cex=.4)
abline(v = Sierit1day$timestamp[brksNSD], col = "red")

#######
# FTP #
#######
brksL <- findpath(Sierit.lav[[4]], Sierit.optseg[4],plotit=F)
plot(attributes(Sierit.fpt[[1]])$date,Sierit.fpt[[1]]$r4, type="l",xlab="", ylab="Days", main="FTP")
NAsBeggining <- rle(is.na(Sierit.fpt[[1]]$r4))
brksFTP <- unlist(brksL)+NAsBeggining$lengths[1]
abline(v = attributes(Sierit.fpt[[1]])$date[brksFTP], col = "red")

########
# BPMM #
########
plot.segments(SieritT.partition, xlab="", ylab="Step length (m)", main="BPMM")

########
# EMbC #
########
plot(bc@pth$dTm, rep(1,length(bc@pth$dTm)), ylab='',ylim=c(0,1),yaxt="n", type="n",main="EMbC")
rect(bc@pth$dTm[-length(bc@pth$dTm)], 0, bc@pth$dTm[-1], 1, col=bc@C[bc@A], border=NA)
legend("topleft",fill=bc@C[sort(unique(bc@A))],legend=sort(unique(bc@A)))

########
# BCPA #
########
plot(Sierit.ws, type="smooth", threshold = 10, legend = FALSE, main="smooth BCPA")
plot(Sierit.ws, type="flat",clusterwidth = 6, legend = FALSE, main="flat BCPA")

# ## alternative way for the flat bcpa
# cps <- ChangePointSummary(Sierit.ws, clusterwidth = 6, tau = F)
# str(cps[[1]])
# plot(Sierit.ws$t.POSIX,Sierit.ws$x, pch=20)
# abline(v = cps[[1]]$middle.POSIX, col = "red")


################
# HMM - 2state #
################
states2 <- viterbi(doubleStateModel)
plot(SieritT.reg$Time, rep(1,length(SieritT.reg$Time)), ylab='',ylim=c(0,1),yaxt="n", type="n", main="HMM - 2 state")
rect(SieritT.reg$Time[-length(SieritT.reg$Time)], 0, SieritT.reg$Time[-1], 1, col=c("gold","blue")[states2], border=NA)
legend("topleft",fill=c("gold","blue"), legend=unique(states2))

################
# HMM - 3state #
################
states3 <- viterbi(tripleStateModel)
plot(SieritT.reg$Time, rep(1,length(SieritT.reg$Time)), ylab='',ylim=c(0,1),yaxt="n",xlab="", type="n", main="HMM - 2 state")
rect(SieritT.reg$Time[-length(SieritT.reg$Time)], 0, SieritT.reg$Time[-1], 1, col=c("gold","blue","green3")[as.factor(states3)], border=NA)
legend("topleft",fill=c("gold","blue","green3"), legend=sort(unique(states3)))

dev.off()





################ plot segments on tracks ############################################

library(colorRamps)

pdf("segmentationPlots_Tracks.pdf", width=(4*2),height=(4*4))
par(mfrow=c(4,2))

#######
# NSD #
#######
# finding break points in the NSD using the lavielle function
lvNSD <- lavielle(drop_units(SieritNSD), Lmin=9, Kmax=12, type = "var")
brksNSDl <-  findpath(lvNSD, 12, plotit = F)
brksNSD <- unlist(brksNSDl)
tscutNSD <- mt_time(Sierit1day)[brksNSD]
# assigning unique ID to each segment
Sierit1day$nsdcat <- NA
for(x in seq(1,length(tscutNSD),2)){
  Sierit1day$nsdcat[mt_time(Sierit1day)>=tscutNSD[x] & mt_time(Sierit1day)<=tscutNSD[x+1]] <- paste0("Seg_",x)
}
# coloring segments by category
Sierit1day$segments <- mt_segments(Sierit1day)
ggplot(Sierit1day)+
  theme_void() +
  geom_sf(aes(geometry=segments,color=nsdcat))


#######
# FPT #
#######
# getting the break points
brksFTPL <- findpath(Sierit.lav[[4]], Sierit.optseg[4],plotit=F)
# finding those pts of the beggining of the track not included in the FTP
NAsBeggining <- rle(is.na(Sierit.fpt[[1]]$r4))
brksFTP <- unlist(brksFTPL)+NAsBeggining$lengths[1]
tscutFTP <- attributes(Sierit.fpt[[1]])$date[brksFTP]
# assigning "excluded" to those pts that were excluded from the FTP
Sierit.prj$ftpcat <- NA
Sierit.prj$ftpcat[mt_time(Sierit.prj)<=tscutFTP[1]] <- "excluded"
# assigning unique ID to each segment
for(x in seq(1,length(tscutFTP),2)){
  Sierit.prj$ftpcat[mt_time(Sierit.prj)<=tscutFTP[x+1] & mt_time(Sierit.prj)>tscutFTP[x]] <- "Seg1" <- paste0("Seg_",x)
}
# coloring segments by category
Sierit.prj$segments <- mt_segments(Sierit.prj)
ggplot(Sierit.prj)+
  theme_void() +
  geom_sf(aes(geometry=segments,color=ftpcat))


########
# BPMM #
########
# getting the break points
splits <- SieritT.partition$Phases.bpm
# assigning unique ID to each segment
Sierit.prj$bpmmcat <- NA
for(x in 1:nrow(splits)){
  Sierit.prj$bpmmcat[mt_time(Sierit.prj)<=splits$date.end[x] & mt_time(Sierit.prj)>splits$date.begin[x]] <- paste0("Seg_",x)
}
# coloring segments by category
# Sierit.prj$segments <- mt_segments(Sierit.prj)
ggplot(Sierit.prj)+
  theme_void() +
  geom_sf(aes(geometry=segments,color=bpmmcat))

########
# EMbC #
########
# bursting move object by segment ID
Sierit$embccat <- as.factor(bc@A)
# coloring segments by category
Sierit$segments <- mt_segments(Sierit)
ggplot(Sierit)+
  theme_void() +
  geom_sf(aes(geometry=segments,color=embccat))


###############
# flat BCPA # ## I did not find an easy why of extracting the break points for the smooth bcpa, but forsure the is one, just keep digging into the functions of the package
##############
# getting the break points
cps <- ChangePointSummary(Sierit.ws, clusterwidth = 6, tau = F)
tsbcpa <- mt_time(Sierit)
brkts <- cps[[1]]$middle.POSIX
# assigning unique ID to each segment
Sierit$bcpsSegmts <- NA
for(x in 1:length(brkts)){
  if(x==1){
    Sierit$bcpsSegmts[mt_time(Sierit)<=brkts[x]] <- paste0("Seg_",x)
  }
  if(x%in%c(2:(length(brkts)))){
    Sierit$bcpsSegmts[mt_time(Sierit)>brkts[x-1] & mt_time(Sierit)<=brkts[x]] <- paste0("Seg_",x)
  }
  if(x==length(brkts)){
    Sierit$bcpsSegmts[mt_time(Sierit)>=brkts[x]] <- paste0("Seg_",x+1)
  }
}
# coloring segments by category
# Sierit$segments <- mt_segments(Sierit)
ggplot(Sierit)+
  theme_void() +
  geom_sf(aes(geometry=segments,color=bcpsSegmts))

################
# HMM - 2state #
################
# creating a move2 object with the data used for this method
SieritTLL.reg$track <- "Sierit"
SieritTLL.regM <- mt_as_move2(SieritTLL.reg,
                              coords = c("X","Y"),
                              time_column = "Time",
                              track_id_column="track",
                              crs = "EPSG:4326")
# getting the break points
states2 <- viterbi(doubleStateModel)
SieritTLL.regM$state2 <- states2

# coloring segments by category
SieritTLL.regM$segments <- mt_segments(SieritTLL.regM)
ggplot(SieritTLL.regM)+
  theme_void() +
  geom_sf(aes(geometry=segments,color=state2))


################
# HMM - 3state #
################
# getting the break points
states3 <- viterbi(tripleStateModel)
SieritTLL.regM$state3 <- states3
# coloring segments by category
SieritTLL.regM$segments <- mt_segments(SieritTLL.regM)
ggplot(SieritTLL.regM)+
  theme_void() +
  geom_sf(aes(geometry=segments,color=states3))
dev.off()

