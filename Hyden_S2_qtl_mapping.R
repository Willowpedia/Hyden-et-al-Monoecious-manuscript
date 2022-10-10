##########################################################################################################
##                                                 ________________________________________________     ##
##   Craig H. Carlson                             /    _                                           \    ##
##   Research Geneticist                         |   /   \                           ,-.            |   ##
##   Cereal Crops Research                       |   \ _ /                           oo )           |   ##
##   United States Dept of Agriculture           |                                   \=  _          |   ##
##   Agricultural Research Service               |                             ``, /`_ _'.\         |   ##
##   Edward T. Schafer Ag Research Center        |           ,----.             \\//|/.\| \\        |   ##
##   1616 Albrect Blvd N                         |          /[O]   \             \/  \ /  ||        |   ##
##   Fargo, ND 58102-2765                        |        _:[]][o][]_   ~           |\"/| ||        |   ##
##                                               |       | |=======| |    ~         \ _ / ,#        |   ##
##   Office: 701-239-1344                        |       |_|  [0]  |_|     ~~       || ||           |   ##
##   Mobile: 701-205-5715                        |       |||  [0]  |||     ~~~~     || ||  ~        |   ##
##   E-mail: craig.h.carlson@usda.gov            |       |||__[_]__|||   ~~~ ~      [| []  ~~       |   ##
##                                               |       |~|\\___//|~|       ~~     || ||  ~~~      |   ##
##   github.com/craighcarlson                    |       /==\ /=\ /==\ ~~~     _ __/_]_[]__ ~       |   ##
##   researchgate.net/profile/Craig-Carlson-2    |  _  __[__]_[ ]_[__]___ _   ~~  ~~    ~           |   ##
##                                               | ~ ~~      ~~~   ~~   ~ ~~~   ~   ~~~   ~~~~~  ~~ |   ##
##   "Do right and feed everyone"                 \________________________________________________/    ##
##                                                                                                      ##
##   Project: scanPlotQTL.R                                                                             ##
##   Date: 2021-09-27                                                                                   ##
##   R version: 4.1.0                                                                                   ##
##                                                                                                      ##  
##########################################################################################################

setwd("C:/Users/Craig.H.Carlson/Desktop/Cornell/2021_Carlson/Manuscripts/Hyden_TBD")

library(qtl)
library(RColorBrewer)

# Read in R/qtl cross to pull map
mapthis <- read.cross(
  format="csv",
  file="S2_8.csv",
  crosstype="f2",
  na.strings=c("-"),
  genotypes=c("AA", "AB", "BB"),
  alleles=c("A", "B"),
  estimate.map=FALSE,
  convertXdata=FALSE
  )

# Format: qtl.fit.list[, c(trait, chr, pos, lod, ci.lo, ci.hi)]
qtl.fit.list <- read.delim("S2_8.qtl.fit.list.txt")

# Prepare some things:
map <- qtl::pull.map(mapthis)
n.chrom <- length(map)
c.list <- seq(1, n.chrom*2, by=2)
map <- lapply(map, function(y) y - y[1])
max.length <- max(unlist(lapply(map, max)))
xlim <- c(1, round(c(n.chrom + NROW(qtl.fit.list)*0.75), 0)) + c(-0.5, 0.5)
ylim <- c(max.length + 30, -20)
chrList <- gtools::mixedsort(as.character(unique(qtl.fit.list$chr)))

# Optional: Color QTL intervals for specific traits
qtlColors <- brewer.pal(length(unique(qtl.fit.list$Trait)), "Dark2")
qtl.fit.list$qtlIntColor <- NA
qtl.fit.list$qtlIntColor[qtl.fit.list$Trait == "LFA"] <- qtlColors[1]
qtl.fit.list$qtlIntColor[qtl.fit.list$Trait == "LFL"] <- qtlColors[2]
qtl.fit.list$qtlIntColor[qtl.fit.list$Trait == "LFP"] <- qtlColors[3]
qtl.fit.list$qtlIntColor[qtl.fit.list$Trait == "SPAD"] <- qtlColors[4]
qtl.fit.list$qtlIntColor[qtl.fit.list$Trait == "RATIO"] <- qtlColors[5]
qtl.fit.list$qtlIntColor[qtl.fit.list$Trait == "SEX"] <- qtlColors[6]

# Make an empty plot
par(mar=c(1.2, 4, 2, 1))
plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="Location (cM)", yaxt='n', xaxt='n', yaxs="i", bty='n')

# Add y-axis elements
x.seqs <- seq(0, round(max(unlist(lapply(pull.map(mapthis), max)))), 50)
axis(side=2, labels=FALSE, at=x.seqs, tck=c(-0.015), las=1, lwd=1, lend=2)
for (i in x.seqs) { mtext(side=2, at=i, text=i, line=0.75, las=1, cex=0.85, adj=1, padj=0.425) }

# Loop through chromosomes
for (i in 1:length(chrList)) {
  # Subset by chromosome
  c.pheno <- subset(qtl.fit.list, chr == chrList[i])
  c.pheno <- c.pheno[order(c.pheno$pos, decreasing=FALSE), ]
  # Figure out where to place linkage groups based on number of previous LGs and QTL
  nrow.to.chr <- NROW(qtl.fit.list[qtl.fit.list$chr %in% chrList[1:c(i-1)], ])
  chr.loc <- ifelse(i > 1, c(nrow.to.chr + c(i * 1.25)), 2)
  # Prints plot x-axis location
  print(paste0("LG", chrList[i], " @ ", chr.loc, "cM"))
  # Axis elements
  axis(side=3, at=chr.loc, labels=chrList[i], tick=FALSE, line=c(-1), font=2)
  segments(x0=(chr.loc - 0.4), y0=map[[chrList[i]]], x1=(chr.loc + 0.4), y1=map[[chrList[i]]], lwd=0.5, col=adjustcolor("black", 0.85), lend=1)
  rect(xleft=(chr.loc - 0.4), ybottom=min(map[[chrList[i]]]), xright=(chr.loc + 0.4), ytop=max(map[[chrList[i]]]), lwd=1, lend=2, border="black")
  # Loop through phenotypes
  for (j in 1:nrow(c.pheno)) {
    # Set positions
    x.left <- (chr.loc + (j*0.75))
    x.right <- (x.left + 0.5)
    y.bottom <- (c.pheno$ci.hi[j] + (j * 0.75))
    y.top <- c.pheno$ci.lo[j]
    y.peak <- c.pheno$pos[j]
    # Make sure QTL intervals don't exceed chr len
    if (y.bottom > max.length) { y.bottom <- max.length }
    # QTL LOD support interval
    lod.color <- c.pheno$qtlIntColor[j] 
    rect(xleft=x.left, ybottom=y.bottom, xright=x.right, ytop=y.top, col=lod.color, lend=2)
    # Peak QTL
    par(lend=1)
    segments(x0=x.left, y0=y.peak, x1=x.right, y1=y.peak, lwd=2, lend=1)
    par(lend=2)
    # Name of phenotype
    text(x=(x.left + 0.3), y=(y.bottom + 5), labels=c.pheno$Trait[j], offset=0, cex=0.6, srt=90, pos=2)
    # LOD value
    text(x=(x.left + 0.25), y=(y.top - 3), labels=round(c.pheno$lod[j], 1), offset=0, cex=0.5, srt=90, pos=4)
  }
}

# Optional:
legend(
  legend=c("Leaf area (LFA)", 
           "Leaf length (LFL)", 
           "Leaf perimeter (LFP)", 
           "Chlorophyll content (SPAD)", 
           "Catkin sex ratio (RATIO)", 
           "Individual sex (SEX)"), 
  fill=qtlColors, x.intersp=0.65, y.intersp=0.9, 
  y=360, x=37, bty="n", border="black", cex=0.6, ncol=2, pt.lwd=1
  )

##########################################################################################################
