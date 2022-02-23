#-------------------------------
# R funtions for OTU analysis.# MASS，epicalc and clusterSim were deleted after consulting with Shi (03072020).
#-------------------------------
p <- c("ade4","calibrate","reshape2","randomForest","candisc","cluster","fpc","plyr","lattice","gplots","ggplot2","squash","RColorBrewer","pheatmap","pROC","ROCR","vegan","fmsb")
usePackage <- function(p) {
if (!is.element(p, installed.packages()[,1]))
install.packages(p, dep = TRUE, repos = "http://cran.us.r-project.org")
suppressWarnings(suppressMessages(invisible(require(p, character.only = TRUE))))}
invisible(lapply(p, usePackage))

#clustering for "types"
#-------------------------------
## KL/JS divergence measure for relative-abundance(density/frequency) data
#-------------------------------
JSD<-function(object, eps=10^-4, overlap=TRUE,...)
{
    if(!is.numeric(object))
        stop("object must be a numeric matrix\n")
    
    z <- matrix(NA, nrow=ncol(object), ncol=ncol(object))
    colnames(z) <- rownames(z) <- colnames(object)
    
    w <- object < eps
    if (any(w)) object[w] <- eps
## If you takes as input a matrix of density values 
## with one row per observation and one column per 
## distribution, add following statement below.
 # object <- sweep(object, 2, colSums(object) , "/")
    
    for(k in seq_len(ncol(object)-1)){
      for(l in 2:ncol(object)){
        ok <- (object[, k] > eps) & (object[, l] > eps)
        if (!overlap | any(ok)) {
          m=0.5*(object[,k]+object[,l])
          z[k,l] <- sqrt(0.5*sum(object[,k] *(log(object[,k]) - log(m)))+0.5*sum(object[,l] *(log(object[,l]) - log(m))))
          z[l,k] <- sqrt(0.5*sum(object[,l] *(log(object[,l]) - log(m)))+0.5*sum(object[,k] *(log(object[,k]) - log(m))))
        }
      }
    }
    diag(z)<-0
    z
}

#-------------------------------
# Log2 transformation of data matrix
#-------------------------------
log.mat<-function(mat,base=2){
mat[mat==0]<-0.000001
log.mat<-log(mat,base)
return(log.mat)
}
#-------------------------------
PAM.best<-function(matrix,dist.mat){
#-------------------------------
# matrix: 
#         row.names	Sample_id
#         col.names	Varibles
# For example, data should be organized like this:
# Sample_id	V1	V2	V3	etc...
# sample_0001	15	6	25
# sample_0002	7	9	32
# etc...
#-------------------------------
 if(!is.numeric(matrix))
        stop("matrix must be a numeric matrix\n")
 if(!is.numeric(dist.mat) && class(dist.mat)=="dist")
        stop("dist.mat must be numeric distance matrix\n")
#-------------------------------
# nc - number_of_clusters
#-------------------------------
min_nc=2
if(nrow(matrix)>20){
max_nc=20} else {
max_nc=nrow(matrix)-1}
res <- array(0,c(max_nc-min_nc+1, 2))
res[,1] <- min_nc:max_nc
siavgs <- array(0,c(max_nc-min_nc+1, 2))
siavgs[,1] <- min_nc:max_nc
clusters <- NULL
for (nc in min_nc:max_nc)
{
cl <- pam(dist.mat, nc, diss=TRUE)
res[nc-min_nc+1,2] <- CH <- index.G1(matrix,cl$cluster,d=dist.mat,centrotypes="medoids")
siavgs[nc-1,2]<-cl$silinfo$avg.width
clusters <- rbind(clusters, cl$cluster)
}
CH_nc<-(min_nc:max_nc)[which.max(res[,2])]
Si_nc<-(min_nc:max_nc)[which.max(siavgs[,2])]
print(paste("max CH for",CH_nc,"clusters=",max(res[,2])))
print(paste("max Si for",Si_nc,"clusters=",max(siavgs[,2])))

CH_cluster<-clusters[which.max(res[,2]),]
Si_cluster<-clusters[which.max(siavgs[,2]),]
objectList      <- list()
    objectList$min_nc <- min_nc
    objectList$max_nc <- max_nc
    objectList$CH     <- CH_cluster
    objectList$Si     <- Si_cluster
    objectList$CH_nc  <- CH_nc
    objectList$Si_nc  <- Si_nc
    objectList$res    <- res
    objectList$siavgs <- siavgs
    return(objectList)
}

#--------------------------------------------------
rangeScaling <- function(v) {
themin <- min(v)
themax <- max(v) 
v <- (v- themin) / (themax - themin)
return(v)
}
#--------------------------------------------------
paretoScaling <- function(v) {
themean <- mean(v)
thesd <- sd(v) 
v <- (v- themean) / sqrt(thesd)
return(v)
}
#--------------------------------------------------
PlotCorrHeatMap<-function(cor.method, colors, data){
    main <- xlab <- ylab <- NULL;
    #print (dim(data));
    #write.csv(data,file="test.csv"); # debug
    
    # Decide size of margins and length of labels
    labelLen<-max(nchar(colnames(data)))
    margins<-c(12, 12)

    if(ncol(data) > 500){
        filter.val <- apply(data.matrix(data), 2, IQR, na.rm=T);
        rk <- rank(-filter.val, ties.method='random');
        data <- as.data.frame(data[,rk < 500]);
        
        print("Data is reduced to 500 vars ..");
    }

    #colnames(data)<-substr(colnames(data), 1, labelLen);
     corr.mat<-cor(data, method=cor.method);

    # set up parameter for heatmap
    suppressMessages(require(RColorBrewer));
    suppressMessages(require(gplots));
    if(colors=="jet"){
        colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(256) 
    }else if(colors=="gbr"){
        colors <- colorRampPalette(c("green", "black", "red"), space="rgb")(256);
    }else if(colors == "heat"){
        colors <- heat.colors(256);
    }else if(colors == "topo"){
        colors <- topo.colors(256);
    }else{
        colors <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256));
    }
    
    heatmap<-heatmap.2(corr.mat,
             Rowv=TRUE,
             Colv=TRUE,
            #dendrogram= c("none"),
             distfun = dist,
             hclustfun = hclust,
             xlab = xlab,
             ylab = ylab,
             key=TRUE,
             keysize=0.8, # size of the key in the chart
             trace="none",
             density.info=c("none"),
             margins=margins,
            #col=brewer.pal(10,"PiYG")
             main = main,
             col=colors
             )
    objectList      <- list()
    objectList$heatmap  <- heatmap
    objectList$corr.mat <- corr.mat
    return(objectList)
}
#--------------------------------------------------
CleanData <-function(bdata, removeNA=T, removeNeg=T){
    if(sum(bdata==Inf)>0){
        inx <- bdata == Inf;
        bdata[inx] <- NA;
        bdata[inx] <- max(bdata, na.rm=T)*2
    }
    if(sum(bdata==-Inf)>0){
        inx <- bdata == -Inf;
        bdata[inx] <- NA;
        bdata[inx] <- min(bdata, na.rm=T)/2
    }
    if(removeNA){
        if(sum(is.na(bdata))>0){
            bdata[is.na(bdata)] <- min(bdata, na.rm=T)/2
        }
    }
    if(removeNeg){
        if(sum(bdata<=0) > 0){
            inx <- bdata <= 0;
            bdata[inx] <- NA;
            bdata[inx] <- min(bdata, na.rm=T)/2
        }
    }
    bdata;
}
#--------------------------------------------------
#--------------------------------------------------
PlotRF.VIP<-function(rf, metadata){
	vip.score <- rev(sort(rf$importance[,"MeanDecreaseAccuracy"]));
	names(vip.score) = substr(names(vip.score), 1, 12);
    if(nrow(rf$importance)>15){
          cmpd.num <- 15;
          cap <- "Top 15 Important Features";
    }else{
          cmpd.num <- nrow(rf$importance);
          cap <- "Important Features";
    }
 	op<-par(mar=c(5,8,4,2));
	dotchart(sort(vip.score[1:cmpd.num]), xlab="MeanDecreaseAccuracy", bg="limegreen", main=cap);
  #mtext(GetVariableLabel(), 2, 7);
	par(op); # return to default setting
}
#--------------------------------------------------
PlotRF.Classify<-function(rf, metadata){
	op<-par(mar=c(4,4,3,2));
    cols <- rainbow(length(levels(factor(metadata)))+1);
    plot(rf, main="Random Forest classification", col=cols);
    legend("topright", legend = c("Overall", levels(factor(metadata))),lty=2, lwd=1, col=cols);
	par(op); # return to default setting
}
#--------------------------------------------------
PlotRF.Outlier<-function(rf, metadata){
    cols <- rainbow(nlevels(as.factor(metadata)))[unclass=metadata];

    # Use a subset of returned colors to exclude samples whose metadata value
    # is NA, and also subset the metadata to be used for the legend.
    # Added by David Arndt, May 4, 2012.
    cols <- cols[!is.na(cols)];
   
    uniq.cols <- unique(cols);
    legend.nm <- unique(as.character(metadata));

    distres <- outlier(rf);

    layout(matrix(c(1,2), 1, 2, byrow = TRUE), width=c(4,1));

 	op<-par(mar=c(5,5,4,0));
    plot(distres, type="h", col=cols, xlab="Samples", xaxt="n", ylab="Outlying Measures", bty="n");

    # add sample names to top 5
    rankres <- rank(-abs(distres), ties.method="random");

    inx.x <- which(rankres < 6);
    inx.y <- distres[inx.x];
    nms <- names(distres)[inx.x];
    text(inx.x, inx.y, nms, cex=0.8,pos=ifelse(inx.y >= 0, 3, 1), xpd=T)
 	op<-par(mar=c(5,0,4,1));
    plot.new();
    plot.window(c(0,1), c(0,1));

    legend("center", legend =legend.nm, pch=15, col=uniq.cols);
	par(op); # return to default setting
}
#--------------------------------------------------
PlotRF.mislabeling<-function(mislabeling, metadata){
    cols <- rainbow(nlevels(as.factor(metadata)))[unclass=metadata];
    cols <- cols[!is.na(cols)];
    uniq.cols <- unique(cols);
    legend.nm <- unique(as.character(metadata));
    layout(matrix(c(1,2), 1, 2, byrow = TRUE), width=c(4,1));

	op<-par(mar=c(5,5,4,0));
    plot(mislabeling, type="h", col=cols, xlab="Samples", xaxt="n", ylab="Outlying Measures", bty="n");

    # add sample names to top 5
    rankres <- rank(-abs(mislabeling), ties.method="random");

    inx.x <- which(rankres < 6);
    inx.y <- mislabeling[inx.x];
    nms <- names(mislabeling)[inx.x];
    text(inx.x, inx.y, nms,cex=0.8, pos=ifelse(inx.y >= 0, 3, 1), xpd=T)
 	op<-par(mar=c(5,0,4,1));
    plot.new();
    plot.window(c(0,1), c(0,1));

    legend("center", legend =legend.nm, pch=15, col=uniq.cols);
	par(op); # return to default setting
}
#--------------------------------------------------

BetweenGroup.test <-function(data, group, p.adj.method="bonferroni",paired=FALSE){
# p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")

n_group<-nlevels(group)
        if(!is.numeric(n_group) | n_group==1)
        stop("group must be a numeric and up to two levels\n")
if(n_group==2){
               output1<-matrix(NA,ncol=9,nrow=ncol(data))
               rownames(output1)<-colnames(data)
               colnames(output1)<-c(paste("mean_",levels(group)[1],sep=""),paste("mean_",levels(group)[2],sep=""),
                                   paste("sd_",levels(group)[1],sep=""),paste("sd_",levels(group)[2],sep=""),"Var.test","T-test","Wilcoxon.test",
                                   paste("T-test_",p.adj.method,sep=""),paste("Wilcoxon.test_",p.adj.method,sep=""))
               for(i in 1:ncol(data))
               {
               output1[i,1]<-mean(data[which(group==levels(group)[1]),i])
               output1[i,2]<-mean(data[which(group==levels(group)[2]),i])
               output1[i,3]<-sd(data[which(group==levels(group)[1]),i])
               output1[i,4]<-sd(data[which(group==levels(group)[2]),i])
               output1[i,5]<-var.test(data[,i]~group)$p.value
               if(output1[i,5]<0.01)
               output1[i,6]<-t.test(data[,i]~group,paired=paired)$p.value
               else
               output1[i,6]<-t.test(data[,i]~group, var.equal=T,paired=paired)$p.value
               output1[i,7]<-wilcox.test(data[,i]~group, paired=paired, conf.int=TRUE, exact=FALSE, correct=FALSE)$p.value
               output1[i,8]<-NA
               output1[i,9]<-NA
               }
               output1[,8]<-p.adjust(output1[,6], method = p.adj.method, n = ncol(data))
               output1[,9]<-p.adjust(output1[,7], method = p.adj.method, n = ncol(data))
               
               return(data.frame(output1))
}else{
      output2<-matrix(NA,ncol=n_group+5,nrow=ncol(data))
      rownames(output2)<-colnames(data)
      colnames.output2<-array(NA)
      for(j in 1:ncol(output2)){
      if(j<=n_group){
      colnames.output2[j]<-c(paste("mean_",levels(group)[j],sep=""))
      }else{
      colnames.output2[(n_group+1):(n_group+5)]<-c("Var.test","Oneway-test","Kruskal.test",
                                                    paste("Oneway-test_",p.adj.method,sep=""),paste("Kruskal.test_",p.adj.method,sep=""))
                                                    }
      }
      colnames(output2)<-colnames.output2
      for(i in 1:ncol(data))
      {
      for(j in 1:n_group)
      {
      output2[i,j]<-mean(data[which(group==levels(group)[j]),i])
      }
      output2[i,(n_group+1)]<-bartlett.test(data[,i]~group)$p.value
      if(output2[i,(n_group+1)]<0.01)
      output2[i,(n_group+2)]<-oneway.test(data[,i]~group)$p.value
      else
      output2[i,(n_group+2)]<-oneway.test(data[,i]~group, var.equal=T)$p.value
      output2[i,(n_group+3)]<-kruskal.test(data[,i]~group)$p.value
      output2[i,(n_group+4)]<-NA
      output2[i,(n_group+5)]<-NA
      }
      output2[ ,(n_group+4)]<-p.adjust(output2[,(n_group+2)], method = p.adj.method, n = ncol(data))
      output2[ ,(n_group+5)]<-p.adjust(output2[,(n_group+3)], method = p.adj.method, n = ncol(data))
      return(data.frame(output2))
      }
      
      
}

########################################
#### Scatterplot3D
#### adapted for better visualization
#######################################

Plot3D <- function(x, y = NULL, z = NULL, color = par("col"), pch = NULL,
     main = NULL, sub = NULL, xlim = NULL, ylim = NULL, zlim = NULL,
     xlab = NULL, ylab = NULL, zlab = NULL, scale.y = 1, angle = 40,
     axis = TRUE, tick.marks = TRUE, label.tick.marks = TRUE,
     x.ticklabs = NULL, y.ticklabs = NULL, z.ticklabs = NULL,
     y.margin.add = 0, grid = TRUE, box = TRUE, lab = par("lab"),
     lab.z = mean(lab[1:2]), type = "p", highlight.3d = FALSE,
     mar = c(5, 3, 4, 3) + 0.1, col.axis = par("col.axis"),
     col.grid = "grey", col.lab = par("col.lab"), cex.symbols = par("cex"),
     cex.axis = 0.8 * par("cex.axis"), cex.lab = par("cex.lab"),
     font.axis = par("font.axis"), font.lab = par("font.lab"),
     lty.axis = par("lty"), lty.grid = 2, lty.hide = 1,
     lty.hplot = par("lty"), log = "", ...)
     # log not yet implemented
{
    ## Uwe Ligges <ligges@statistik.tu-dortmund.de>,
    ## http://www.statistik.tu-dortmund.de/~ligges
    ##
    ## For MANY ideas and improvements thanks to Martin Maechler!!!
    ## Parts of the help files are stolen from the standard plotting functions in R.

    mem.par <- par(mar = mar)
    x.scal <- y.scal <- z.scal <- 1
    xlabel <- if (!missing(x)) deparse(substitute(x))
    ylabel <- if (!missing(y)) deparse(substitute(y))
    zlabel <- if (!missing(z)) deparse(substitute(z))
    ## verification, init, ...
    if(highlight.3d && !missing(color))
        warning("color is ignored when highlight.3d = TRUE")


    ## color as part of `x' (data.frame or list):
    if(!is.null(d <- dim(x)) && (length(d) == 2) && (d[2] >= 4))
        color <- x[,4]
    else if(is.list(x) && !is.null(x$color))
        color <- x$color

    ## convert 'anything' -> vector
    xyz <- xyz.coords(x=x, y=y, z=z, xlab=xlabel, ylab=ylabel, zlab=zlabel,
                      log=log)
    if(is.null(xlab)) { xlab <- xyz$xlab; if(is.null(xlab)) xlab <- "" }
    if(is.null(ylab)) { ylab <- xyz$ylab; if(is.null(ylab)) ylab <- "" }
    if(is.null(zlab)) { zlab <- xyz$zlab; if(is.null(zlab)) zlab <- "" }

    if(length(color) == 1)
        color <- rep(color, length(xyz$x))
    else if(length(color) != length(xyz$x))
        stop("length(color) ", "must be equal length(x) or 1")

    angle <- (angle %% 360) / 90
    yz.f <- scale.y * abs(if(angle < 1) angle else if(angle > 3) angle - 4 else 2 - angle)
    yx.f <- scale.y * (if(angle < 2) 1 - angle else angle - 3)
    if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
        temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
        temp <- xlab;  xlab <- ylab;   ylab <- temp
        temp <- xlim;  xlim <- ylim;   ylim <- temp
    }
    angle.1 <- (1 < angle && angle < 2) || angle > 3
    angle.2 <- 1 <= angle && angle <= 3
    dat <- cbind(as.data.frame(xyz[c("x","y","z")]), col = color)

    ## xlim, ylim, zlim -- select the points inside the limits
    if(!is.null(xlim)) {
        xlim <- range(xlim)
        dat <- dat[ xlim[1] <= dat$x & dat$x <= xlim[2] , , drop = FALSE]
    }
    if(!is.null(ylim)) {
        ylim <- range(ylim)
        dat <- dat[ ylim[1] <= dat$y & dat$y <= ylim[2] , , drop = FALSE]
    }
    if(!is.null(zlim)) {
        zlim <- range(zlim)
        dat <- dat[ zlim[1] <= dat$z & dat$z <= zlim[2] , , drop = FALSE]
    }
    n <- nrow(dat)
    if(n < 1) stop("no data left within (x|y|z)lim")

    y.range <- range(dat$y[is.finite(dat$y)])

### 3D-highlighting / colors / sort by y
    if(type == "p" || type == "h") {
        y.ord <- rev(order(dat$y))
        dat <- dat[y.ord, ]
        if(length(pch) > 1)
            if(length(pch) != length(y.ord))
                stop("length(pch) ", "must be equal length(x) or 1")
            else pch <- pch[y.ord]
        daty <- dat$y
        daty[!is.finite(daty)] <- mean(daty[is.finite(daty)])
        if(highlight.3d && !(all(diff(daty) == 0)))
            dat$col <- rgb(seq(0, 1, length = n) * (y.range[2] - daty) / diff(y.range), g=0, b=0)
    }

### optim. axis scaling
    p.lab <- par("lab")
    ## Y
    y.range <- range(dat$y[is.finite(dat$y)], ylim)
    y.prty <- pretty(y.range, n = lab[2],
        min.n = max(1, min(.5 * lab[2], p.lab[2])))
    y.scal <- round(diff(y.prty[1:2]), digits = 12)
    y.add <- min(y.prty)
    dat$y <- (dat$y - y.add) / y.scal
    y.max <- (max(y.prty) - y.add) / y.scal
    if(!is.null(ylim)) y.max <- max(y.max, ceiling((ylim[2] - y.add) / y.scal))
#    if(angle > 2) dat$y <- y.max - dat$y  ## turn y-values around
    ## X
    x.range <- range(dat$x[is.finite(dat$x)], xlim)
    x.prty <- pretty(x.range, n = lab[1],
        min.n = max(1, min(.5 * lab[1], p.lab[1])))
    x.scal <- round(diff(x.prty[1:2]), digits = 12)
    dat$x <- dat$x / x.scal
    x.range <- range(x.prty) / x.scal
    x.max <- ceiling(x.range[2])
    x.min <-   floor(x.range[1])
    if(!is.null(xlim)) {
        x.max <- max(x.max, ceiling(xlim[2] / x.scal))
        x.min <- min(x.min,   floor(xlim[1] / x.scal))
    }
    x.range <- range(x.min, x.max)
    ## Z
    z.range <- range(dat$z[is.finite(dat$z)], zlim)
    z.prty <- pretty(z.range, n = lab.z,
        min.n = max(1, min(.5 * lab.z, p.lab[2])))
    z.scal <- round(diff(z.prty[1:2]), digits = 12)
    dat$z <- dat$z / z.scal
    z.range <- range(z.prty) / z.scal
    z.max <- ceiling(z.range[2])
    z.min <-   floor(z.range[1])
    if(!is.null(zlim)) {
        z.max <- max(z.max, ceiling(zlim[2] / z.scal))
        z.min <- min(z.min,   floor(zlim[1] / z.scal))
    }
    z.range <- range(z.min, z.max)

### init graphics
    plot.new()
    if(angle.2) {x1 <- x.min + yx.f * y.max; x2 <- x.max}
    else        {x1 <- x.min; x2 <- x.max + yx.f * y.max}
    plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
    temp <- strwidth(format(rev(y.prty))[1], cex = cex.axis/par("cex"))
    if(angle.2) x1 <- x1 - temp - y.margin.add
    else        x2 <- x2 + temp + y.margin.add
    plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max))
    if(angle > 2) par("usr" = par("usr")[c(2, 1, 3:4)])
    usr <- par("usr") # we have to remind it for use in closures
    title(main, sub, ...)

### draw axis, tick marks, labels, grid, ...
    xx <- if(angle.2) c(x.min, x.max) else c(x.max, x.min)
    if(grid) {
	## grids
	###################
	# XY wall
        i <- x.min:x.max;
        segments(i, z.min, i + (yx.f * y.max), yz.f * y.max + z.min,
                 col = col.grid, lty = lty.grid);

        i <- 0:y.max;
        segments(x.min + (i * yx.f), i * yz.f + z.min,
                 x.max + (i * yx.f), i * yz.f + z.min,
                 col = col.grid, lty = lty.grid);

	######################
	# XZ wall
	# verticle lines
        temp <- yx.f * y.max;
        temp1 <- yz.f * y.max;
        i <- (x.min + temp):(x.max + temp);
        segments(i, z.min + temp1, i, z.max + temp1,
                 col = col.grid, lty = lty.grid);

	# horizontal lines
        i <- (z.min + temp1):(z.max + temp1);
        segments(x.min + temp, i, x.max + temp, i,
                 col = col.grid, lty = lty.grid)


	##################
	# YZ wall
	# horizontal lines
        i <- xx[2]:x.min;
	mm <- z.min:z.max;
        segments(i, mm, i + temp, mm + temp1,
                 col = col.grid, lty = lty.grid);
	# verticle lines
        i <- 0:y.max;
        segments(x.min + (i * yx.f), i * yz.f + z.min,
                 xx[2] + (i * yx.f), i * yz.f + z.max,
                 col = col.grid, lty = lty.grid)


	# make the axis into solid line
        segments(x.min, z.min, x.min + (yx.f * y.max), yz.f * y.max + z.min,
                 col = col.grid, lty = lty.hide);
        segments(x.max, z.min, x.max + (yx.f * y.max), yz.f * y.max + z.min,
                 col = col.axis, lty = lty.hide);
        segments(x.min + (y.max * yx.f), y.max * yz.f + z.min,
                 x.max + (y.max* yx.f), y.max * yz.f + z.min,
                col = col.grid, lty = lty.hide);
        segments(x.min + temp, z.min + temp1, x.min + temp, z.max + temp1,
                col = col.grid, lty = lty.hide);
        segments(x.max + temp, z.min + temp1, x.max + temp, z.max + temp1,
                col = col.axis, lty = lty.hide);
        segments(x.min + temp, z.max + temp1, x.max + temp, z.max + temp1,
                col = col.axis, lty = lty.hide);
        segments(xx[2], z.max, xx[2] + temp, z.max + temp1,
                col = col.axis, lty = lty.hide);
    }
    if(axis) {
        if(tick.marks) { ## tick marks
            xtl <- (z.max - z.min) * (tcl <- -par("tcl")) / 50
            ztl <- (x.max - x.min) * tcl / 50
            mysegs <- function(x0,y0, x1,y1)
                segments(x0,y0, x1,y1, col=col.axis, lty=lty.axis)
            ## Y
            i.y <- 0:y.max
            mysegs(yx.f * i.y - ztl + xx[1], yz.f * i.y + z.min,
                   yx.f * i.y + ztl + xx[1], yz.f * i.y + z.min)
            ## X
            i.x <- x.min:x.max
            mysegs(i.x, -xtl + z.min, i.x, xtl + z.min)
            ## Z
            i.z <- z.min:z.max
            mysegs(-ztl + xx[2], i.z, ztl + xx[2], i.z)

            if(label.tick.marks) { ## label tick marks
                las <- par("las")
                mytext <- function(labels, side, at, ...)
                    mtext(text = labels, side = side, at = at, line = -.5,
                          col=col.lab, cex=cex.axis, font=font.lab, ...)
                ## X
                if(is.null(x.ticklabs))
                    x.ticklabs <- format(i.x * x.scal)
                mytext(x.ticklabs, side = 1, at = i.x)
                ## Z
                if(is.null(z.ticklabs))
                    z.ticklabs <- format(i.z * z.scal)
                mytext(z.ticklabs, side = if(angle.1) 4 else 2, at = i.z,
                       adj = if(0 < las && las < 3) 1 else NA)
                ## Y
                temp <- if(angle > 2) rev(i.y) else i.y ## turn y-labels around
                if(is.null(y.ticklabs))
                    y.ticklabs <- format(y.prty)
                else if (angle > 2)
                    y.ticklabs <- rev(y.ticklabs)
                text(i.y * yx.f + xx[1],
                     i.y * yz.f + z.min, y.ticklabs,
                     pos=if(angle.1) 2 else 4, offset=1,
                     col=col.lab, cex=cex.axis/par("cex"), font=font.lab)
            }
        }

        ## axis and labels

        mytext2 <- function(lab, side, line, at)
            mtext(lab, side = side, line = line, at = at, col = col.lab,
                  cex = cex.lab, font = font.axis, las = 0)
        ## X
        lines(c(x.min, x.max), c(z.min, z.min), col = col.axis, lty = lty.axis)
        mytext2(xlab, 1, line = 1.5, at = mean(x.range))
        ## Y
        lines(xx[1] + c(0, y.max * yx.f), c(z.min, y.max * yz.f + z.min),
              col = col.axis, lty = lty.axis)
        mytext2(ylab, if(angle.1) 2 else 4, line= 0.5, at = z.min + y.max * yz.f)

        ## Z
        lines(xx[c(2,2)], c(z.min, z.max), col = col.axis, lty = lty.axis)
        mytext2(zlab, if(angle.1) 4 else 2, line= 1.5, at = mean(z.range))

    }

### plot points
    x <- dat$x + (dat$y * yx.f)
    z <- dat$z + (dat$y * yz.f)
    col <- as.character(dat$col)
    if(type == "h") {
        z2 <- dat$y * yz.f + z.min
        segments(x, z, x, z2, col = col, cex = cex.symbols, lty = lty.hplot, ...)
        points(x, z, type = "p", col = col, pch = pch, cex = cex.symbols, ...)
    }
    else points(x, z, type = type, col = col, pch = pch, cex = cex.symbols, ...)

### box-lines in front of points (overlay)
    if(axis && box) {
        lines(c(x.min, x.max), c(z.max, z.max),
              col = col.axis, lty = lty.axis)
        lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + z.max,
              col = col.axis, lty = lty.axis)
        lines(xx[c(1,1)], c(z.min, z.max), col = col.axis, lty = lty.axis)
    }


    # par(mem.par) # we MUST NOT set the margins back
    ### Return Function Object
    ob <- ls() ## remove all unused objects from the result's enviroment:
    rm(list = ob[!ob %in% c("angle", "mar", "usr", "x.scal", "y.scal", "z.scal", "yx.f",
        "yz.f", "y.add", "z.min", "z.max", "x.min", "x.max", "y.max",
        "x.prty", "y.prty", "z.prty")])
    rm(ob)
    invisible(list(
        xyz.convert = function(x, y=NULL, z=NULL) {
            xyz <- xyz.coords(x, y, z)
            if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
                temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
            }
            y <- (xyz$y - y.add) / y.scal
            return(list(x = xyz$x / x.scal + yx.f * y,
                y = xyz$z / z.scal + yz.f * y))
        },
        points3d = function(x, y = NULL, z = NULL, type = "p", ...) {
            xyz <- xyz.coords(x, y, z)
            if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
                temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
            }
            y2 <- (xyz$y - y.add) / y.scal
            x <- xyz$x / x.scal + yx.f * y2
            y <- xyz$z / z.scal + yz.f * y2
            mem.par <- par(mar = mar, usr = usr)
            on.exit(par(mem.par))
            if(type == "h") {
                y2 <- z.min + yz.f * y2
                segments(x, y, x, y2, ...)
                points(x, y, type = "p", ...)
            }
            else points(x, y, type = type, ...)
        },
        plane3d = function(Intercept, x.coef = NULL, y.coef = NULL,
            lty = "dashed", lty.box = NULL, ...){
            if(!is.atomic(Intercept) && !is.null(coef(Intercept))) Intercept <- coef(Intercept)
            if(is.null(lty.box)) lty.box <- lty
            if(is.null(x.coef) && length(Intercept) == 3){
                x.coef <- Intercept[if(angle > 2) 3 else 2]
                y.coef <- Intercept[if(angle > 2) 2 else 3]
                Intercept <- Intercept[1]
            }
            mem.par <- par(mar = mar, usr = usr)
            on.exit(par(mem.par))
            x <- x.min:x.max
            ltya <- c(lty.box, rep(lty, length(x)-2), lty.box)
            x.coef <- x.coef * x.scal
            z1 <- (Intercept + x * x.coef + y.add * y.coef) / z.scal
            z2 <- (Intercept + x * x.coef +
                (y.max * y.scal + y.add) * y.coef) / z.scal
            segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, lty = ltya, ...)
            y <- 0:y.max
            ltya <- c(lty.box, rep(lty, length(y)-2), lty.box)
            y.coef <- (y * y.scal + y.add) * y.coef
            z1 <- (Intercept + x.min * x.coef + y.coef) / z.scal
            z2 <- (Intercept + x.max * x.coef + y.coef) / z.scal
            segments(x.min + y * yx.f, z1 + y * yz.f,
                x.max + y * yx.f, z2 + y * yz.f, lty = ltya, ...)
        },

        wall3d = function(Intercept, x.coef = NULL, y.coef = NULL,
            lty = "dashed", lty.box = NULL, ...){
            if(!is.atomic(Intercept) && !is.null(coef(Intercept))) Intercept <- coef(Intercept)
            if(is.null(lty.box)) lty.box <- lty
            if(is.null(x.coef) && length(Intercept) == 3){
                x.coef <- Intercept[if(angle > 2) 3 else 2]
                y.coef <- Intercept[if(angle > 2) 2 else 3]
                Intercept <- Intercept[1]
            }
            mem.par <- par(mar = mar, usr = usr)
            on.exit(par(mem.par))
            x <- x.min:x.max
            ltya <- c(lty.box, rep(lty, length(x)-2), lty.box)
            x.coef <- x.coef * x.scal
            z1 <- (Intercept + x * x.coef + y.add * y.coef) / z.scal
            z2 <- (Intercept + x * x.coef +
                (y.max * y.scal + y.add) * y.coef) / z.scal
            segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, lty = ltya, ...)
            y <- 0:y.max
            ltya <- c(lty.box, rep(lty, length(y)-2), lty.box)
            y.coef <- (y * y.scal + y.add) * y.coef
            z1 <- (Intercept + x.min * x.coef + y.coef) / z.scal
            z2 <- (Intercept + x.max * x.coef + y.coef) / z.scal
            segments(x.min + y * yx.f, z1 + y * yz.f,
                x.max + y * yx.f, z2 + y * yz.f, lty = ltya, ...)
        },
        box3d = function(...){
            mem.par <- par(mar = mar, usr = usr)
            on.exit(par(mem.par))
            lines(c(x.min, x.max), c(z.max, z.max), ...)
            lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + z.max, ...)
            lines(c(0, y.max * yx.f) + x.min, c(0, y.max * yz.f) + z.max, ...)
            lines(c(x.max, x.max), c(z.min, z.max), ...)
            lines(c(x.min, x.min), c(z.min, z.max), ...)
            lines(c(x.min, x.max), c(z.min, z.min), ...)
        }
    ))
}

