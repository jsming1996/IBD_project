
cat("
 Discription: This R script implements Breiman's random forest algorithm (based
              on Breiman and Cutler's original Fortran code) for classification 
              and regression.
               
       Input: (1) Training data
              (2) Metadata of training data
              (3) Test data (The data you want to predict)
              (4) Metadata of test data
              (5) The factor category (i.e. status) both in the two metadata you 
                  provided
              
              Note: The non-overlaping variables between training and test data will
              be screen out for modeling. 
              
      Output: (1) The variables whose importance (mean) were over 0 in the RF 
                  classification.
              (2) The variables whose importances were all over 0.001 in the RF 
                  classification. Heatmap will visualize the correlation between 
                  the selected variables and the factor category you provided.
              (3) The predicted probability of next status of host in test data. 
              
              \n")
cat("
       Usage: setwd(\"$PWD\"); soucre(\"ECCK1.tax.randomforest.pred-prob.R\") \n\n")

# Last update: 20150305
#--------------------------------------------------
# Run randomForest selecting taxa
# Require source("util.R")
#--------------------------------------------------
setwd("/home/jsming/IBD/r_code")
source("util.R")
source("randomforests_util.R")

#--------------------------------------------------
file.opts<-list(
                prefix="gene5",
                filename="gene5.csv", #"2abundace_snp_discovery/abundance_snp.csv" #"1abundance_discovery/abundance.csv",
                metadata.filename="gene_meta.txt", ## "1abundance_discovery/abundance_meta.txt"
                group.type="Status",
                               
                pred.prefix="check", #NA,
                pred.filename="gene5_1.csv", #NA, 
                pred.metadata.filename="gene_meta1.txt", #NA,
                pred.group.type="Status", #NA, #"Status"
                
                fdr_cutoff=0.1,
                rf.cv=TRUE
               )
#--------------------------------------------------
rf.opts<-list(outdir=NA,ntree=5000,verbose=FALSE,errortype="oob",n.oob=10, nfolds=5)

#--------------------------------------------------
filename<-file.opts$filename
metadata.filename<-file.opts$metadata.filename
pred.filename<-file.opts$pred.filename
pred.metadata.filename<-file.opts$pred.metadata.filename

#-------------------------------

#-------------------------------
# data format
#-------------------------------
# matrix: 
#         row.names	Sample_id
#         col.names	Varibles
# For example, data should be organized like this:
# Sample_id	group	V1	V2	etc...
# sample_0001	A	6	25
# sample_0002	B	9	32
# etc...
#-------------------------------

if(!is.na(file.opts$pred.prefix)){
dir.create(outpath<-paste("./",file.opts$prefix,"-",file.opts$pred.prefix,".",file.opts$group.type,"/",sep=""))
con <- file(paste(outpath,file.opts$prefix,"-",file.opts$pred.prefix,".",file.opts$group.type,".",rf.opts$n.oob,".log",sep=""))}else{
dir.create(outpath<-paste("./",file.opts$prefix,".",file.opts$group.type,"/",sep=""))
con <- file(paste(outpath,file.opts$prefix,".",file.opts$group.type,".",rf.opts$n.oob,".log",sep=""))}
#sink(con, append=TRUE)
#sink(con, append=TRUE, type="message")
#-------------------------------
# Training data input
#-------------------------------
    g<-read.csv(filename,header=T,row.names=1)
    g<-g[order(rownames(g)),]
    print(paste("The number of variables : ", ncol(g) ,sep=""))
    #-------------------------------filtering taxa with zero variance
    g<-g[,which(apply(g,2,var)!=0)]
    print(paste("The number of variables (removed variables with zero variance) : ", ncol(g) ,sep=""))
    #-------------------------------
    gmat<-data.matrix(g)
#-------------------------------
# Metadata of training data input
#-------------------------------
    metadata<-read.table(metadata.filename,header=T,sep="\t",row.names=1)
    metadata<-metadata[order(rownames(metadata)),]
    group<-metadata[,file.opts$group.type]; names(group)<-rownames(metadata)

if(!is.na(file.opts$pred.prefix)){
#-------------------------------
# Test data input
#-------------------------------
    pred.g<-read.csv(pred.filename,header=T,row.names=1)
    pred.g<-pred.g[order(rownames(pred.g)),]
    print(paste("The number of variables of prediction data: ", ncol(pred.g) ,sep=""))
    #-------------------------------filtering taxa with X% zero
    NonZero.p<-1
    pred.g<-pred.g[,which(colSums(pred.g==0)<NonZero.p*nrow(pred.g))]
    print(paste("The number of variables (removed variables containing over ", NonZero.p," zero) of prediction data: ", ncol(pred.g) ,sep=""))
    #-------------------------------
    pred.gmat<-data.matrix(pred.g)
#-------------------------------
# Metadata of test data input
#-------------------------------
    pred.metadata<-read.table(pred.metadata.filename,header=T,sep="\t",row.names=1)
    pred.metadata<-pred.metadata[order(rownames(pred.metadata)),]
    pred.group<-pred.metadata[,file.opts$pred.group.type]; names(pred.group)<-rownames(pred.metadata)
#-------------------------------
# To get shared taxa between training and prediction data set
#-------------------------------
    gmatO<-gmat[,colnames(gmat) %in% colnames(pred.gmat)]
    pred.gmatO<-pred.gmat[,colnames(pred.gmat) %in% colnames(gmat)]
    print(paste("The number of variables for prediction in test data: ", ncol(gmatO) ,sep=""))
    #--------------------------------------------------
    x<-gmatO
    y<-group
    }else{
    x<-gmat
    y<-group
    }
    
#--------------------------------------------------
# Ten-folds cross validation in traing data
#--------------------------------------------------
    result<-rf.cross.validation(x,y,nfolds=rf.opts$nfolds, verbose=rf.opts$verbose,ntree=rf.opts$ntree)
    print(paste("The ", rf.opts$nfolds, "-folds cross validation for randomforest classification",sep=""))
    print(result$confusion.matrix)
    kappa.result.cv<-Kappa.test(result$confusion.matrix)
    acc.cv<-sum(diag(result$confusion.matrix))/sum(result$confusion.matrix)
    cat(paste("Cohen's kappa statistics for agreement in CV: ", kappa.result.cv$Result$estimate," -- ",kappa.result.cv$Judgement,"\n",sep=""))
    cat(paste("Accuracy in CV: ",acc.cv ,"\n",sep=""))
    
    #---------------------------   
    sink(paste(outpath,file.opts$prefix,".",file.opts$group.type,".CV",rf.opts$nfolds,".Prob.xls",sep="")); cat("\t"); write.table(result$probabilities,sep="\t",quote=FALSE); sink()
    #----------------------------probability_boxplot
    pdf(paste(outpath,file.opts$prefix,".",file.opts$group.type,".CV",rf.opts$nfolds,".probability_boxplot.pdf",sep=""), height=8, width=4)
    par(mar=c(4,4,3,3))
    boxplot(result$probabilities[,1]~y, ylab=paste("Probability of ",colnames(result$probabilities)[1],sep=""),ylim=c(0,1),outline=FALSE)
    abline(h=0.5,lty=2)
    stripchart(result$probabilities[,1]~y,method="jitter",jitter=.05,vertical=T,add=T,pch=20) 
    dev.off()
    
    if(nlevels(group)==2){
    sens.cv<-diag(result$confusion.matrix)[1]/sum(result$confusion.matrix[,1])
    spec.cv<-diag(result$confusion.matrix)[2]/sum(result$confusion.matrix[,2])
    cat(paste("Sensitivity in CV: ", sens.cv ,"\n",sep=""))
    cat(paste("Specificity in CV: ", spec.cv ,"\n",sep=""))
    #--------------------------------------------------
    #  ROC plot using "ROCR" package
    #--------------------------------------------------
    pred<-prediction(result$probabilities[,2],y)
    perf<-performance(pred,"tpr","fpr")
    auc.tmp <- performance(pred,"auc"); auc <- as.numeric(auc.tmp@y.values)
    cat(paste("AUC in CV: ", auc ,"\n",sep="")) 
    #----------------------------
    pdf(paste(outpath,file.opts$prefix,".CV",rf.opts$nfolds,".Prob.ROCR_plot.pdf",sep=""), height=4, width=4)
    par(mar=c(4,4,3,3))
    plot(perf,main="Training-CV",col=2,lwd=2)
    abline(a=0,b=1,lwd=2,lty=2,col="gray")
    text(0.8,0.2,paste("AUC: ", formatC(auc,digits=2,format="f"),sep=""))
    dev.off()
    #----------------------------
    pdf(paste(outpath,file.opts$prefix,".CV",rf.opts$nfolds,".ROCR.accuracy_VS_cutoff.plot.pdf",sep=""), height=4, width=4)
    par(mar=c(4,4,3,3))
    acc.perf <- performance(pred, "acc")
    plot(acc.perf, avg= "vertical")
    dev.off()
    #--------------------------------------------------
    #  ROC plot using "pROC" package
    #--------------------------------------------------
    pdf(paste(outpath,file.opts$prefix,".CV",rf.opts$nfolds,".Prob.pROC.ci.pdf",sep=""),width=5,height=5)
    rocobj <- plot.roc(y, result$probabilities[,2],main="", percent=TRUE,ci=TRUE) # print the AUC (will contain the CI)
    ciobj <- ci.se(rocobj,specificities=seq(0, 100, 5)) # over a select set of specificities
    plot(ciobj, type="shape", col="#1c61b6AA", print.auc=TRUE) # plot as a blue shape
    text(50,20,paste("AUC = ",formatC(rocobj$auc,digits=2,format="f"),sep=""),pos=4)
    ci.lower<-formatC(rocobj$ci[1],digits=2,format="f")
    ci.upper<-formatC(rocobj$ci[3],digits=2,format="f")
    text(50,10,paste("95% CI: ",ci.lower,"-",ci.upper,sep=""),pos=4)
    dev.off()
    }
#--------------------------------------------------
# Out-of-bag classification in traing data and prediction in test data
#--------------------------------------------------
#--------------------------------------------------
# (1) Out-of-bag classification in traing data 
#--------------------------------------------------
    oob.result <- rf.out.of.bag(x, y, verbose=rf.opts$verbose, ntree=rf.opts$ntree)
    print(paste("The out-of-bag (3-fold CV) result for randomforest classification",sep=""))
    #----------------------------
    print(oob.result$confusion.matrix)
    kappa.result.oob<-Kappa.test(oob.result$confusion.matrix)
    cat(paste("Cohen's kappa statistics for agreement in out-of-bag classification: ", kappa.result.oob$Result$estimate," -- ",kappa.result.oob$Judgement,"\n",sep=""))
    #---------------------------   
    sink(paste(outpath,file.opts$prefix,".",file.opts$group.type,".oob.Prob.xls",sep="")); cat("\t"); write.table(oob.result$probabilities,sep="\t",quote=FALSE); sink()
#--------------------------------------------------
# (2) Prediction in test data
#--------------------------------------------------
    if(!is.na(file.opts$pred.prefix)){
    #----------------------------
    predicted.prob<-predict(oob.result$rf.model,pred.gmatO,type="prob")      
    sink(paste(outpath,file.opts$pred.prefix,".",file.opts$pred.group.type,".oob.Pred_Prob.xls",sep="")); cat("\t"); write.table(predicted.prob,sep="\t",quote=FALSE); sink()
    #----------------------------probability_boxplot
    pdf(paste(outpath,file.opts$pred.prefix,".",file.opts$pred.group.type,".oob.prob_boxplot.pdf",sep=""), height=8, width=4)
    par(mar=c(4,4,3,3))
    print(predicted.prob[,1])
    print(pred.group)
    boxplot(predicted.prob[,1]~pred.group,xlab=file.opts$pred.prefix, ylab=paste("Probability of ",colnames(predicted.prob)[1],sep=""),ylim=c(0,1),outline=FALSE)
    abline(h=0.5,lty=2)
    stripchart(predicted.prob[,1]~pred.group,method="jitter",jitter=.05,vertical=T,add=T,pch=20) 
    dev.off()
    
    if(nlevels(group)==2){
    #----------------------------probability_xyplot
    pdf(paste(outpath,file.opts$pred.prefix,".",file.opts$pred.group.type,".oob.prob_xyplot.pdf",sep=""), height=6, width=6)
    par(mar=c(4,4,3,3))
    plot(predicted.prob[,1],predicted.prob[,2],xlab=paste("Probability of ",colnames(predicted.prob)[1],sep=""),
                                               ylab=paste("Probability of ",colnames(predicted.prob)[2],sep=""),
                                               main=file.opts$pred.prefix, xlim=c(0,1),ylim=c(0,1),col="white",pch=21,bg="black")
    abline(v=0.5,h=0.5,lty=2)
    dev.off()
    #----------------------------
    if(nlevels(pred.group)==2){
    #--------------------------------------------------
    #  ROC plot using "ROCR" package
    #--------------------------------------------------
    pred<-prediction(predicted.prob[,2],pred.group)
    perf<-performance(pred,"tpr","fpr")
    auc.tmp <- performance(pred,"auc"); auc <- as.numeric(auc.tmp@y.values)
    #----------------------------
    pdf(paste(outpath,file.opts$pred.prefix,".",file.opts$pred.group.type,".oob.ROCR_plot.pdf",sep=""), height=4, width=4)
    par(mar=c(4,4,3,3))
    plot(perf,main="Predicted",col=2,lwd=2)
    abline(a=0,b=1,lwd=2,lty=2,col="gray")
    text(0.8,0.2,paste("AUC: ", formatC(auc,digits=2,format="f"),sep=""))
    dev.off()
    #----------------------------
    pdf(paste(outpath,file.opts$prefix,".",file.opts$pred.group.type,".oob.ROCR.accuracy_VS_cutoff.plot.pdf",sep=""), height=4, width=4)
    par(mar=c(4,4,3,3))
    acc.perf <- performance(pred, "acc")
    plot(acc.perf, avg= "vertical")
    dev.off()
    }
    }
    #--------------------------------------------------
    #  probability_boxplot of both train and test data
    #--------------------------------------------------
    if(file.opts$pred.group.type==file.opts$pred.group.type){
    comb.group<-as.factor(c(as.character(group),paste("relative",as.character(pred.group),sep="")))
    comb.group<-factor(comb.group,levels=c(levels(comb.group)[2],levels(comb.group)[3],levels(comb.group)[1]))
    comb.prob<-rbind(oob.result$probabilities,predicted.prob)
    pdf(paste(outpath,file.opts$prefix,"-",file.opts$pred.prefix,".",file.opts$group.type,".probability_boxplot.pdf",sep=""), height=8, width=4)
    par(mar=c(4,4,3,3))
    boxplot(comb.prob[,1]~comb.group, ylab=paste("Probability of ",colnames(comb.prob)[1],sep=""),ylim=c(0,1),outline=FALSE)
    abline(h=0.5,lty=2)
    stripchart(comb.prob[,1]~comb.group,method="jitter",jitter=.05,vertical=T,add=T,pch=20) 
    dev.off()
    
    sink(paste(outpath,file.opts$prefix,"-",file.opts$pred.prefix,".",file.opts$group.type,".Prob.xls",sep="")); cat("\t"); write.table(comb.prob,sep="\t",quote=FALSE); sink()
    }
    }
    
#--------------------------------------------------
# Barplot of importances of important taxa
#--------------------------------------------------
    imps<-oob.result$importances
    imps<-imps[order(imps)]
    imps.cutoff<-imps[which(imps>0)]
    sink(paste(outpath,file.opts$prefix,".",file.opts$group.type,".imps_All_sorted.xls",sep=""));cat("\t");write.table(imps,quote=FALSE,sep="\t");sink()
    sink(paste(outpath,file.opts$prefix,".",file.opts$group.type,".imps_cutoff_sorted.xls",sep=""));cat("\t");write.table(imps.cutoff,quote=FALSE,sep="\t");sink()
    #--------------------------------------------------
    pdf(paste(outpath,file.opts$prefix,".impVars_cutoff.ordered.",length(imps.cutoff),".pdf",sep=""), height=length(imps.cutoff)/4, width=6)
    print(barchart(imps.cutoff,xlab="Mean decrease in accuracy",box.ratio=4))
    dev.off()
    
#--------------------------------------------------
if(file.opts$rf.cv){
#--------------------------------------------------
#--------------------------------------------------
# Estimate the minErr of RF model and the top discriminatory taxa
#--------------------------------------------------
    results <- replicate(10, rfcv(x, y, step=0.9, cv.fold=10), simplify=FALSE)
    err.cv <- sapply(results, "[[", "error.cv")
    #--------------------------------------------------matplot
    pdf(paste(outpath,file.opts$prefix,".CV.MeanErrPlot.pdf",sep=""), height=6, width=6)
    matplot(results[[1]]$n.var, cbind(rowMeans(err.cv), err.cv), type="p",log="x",
            col=c(2, rep("grey60", ncol(err.cv))), 
            pch=c(19, rep(20, ncol(err.cv))), 
            xlab="Number of variables", ylab="CV Error")
    lines(results[[1]]$n.var,rowMeans(err.cv),col=2)
    breaks<-axTicks(side=1)
    dev.off()
    #--------------------------------------------------ggplot
    ErrSumm<-data.frame(NumVars=results[[1]]$n.var,MeanErr=apply(err.cv,1,mean),MedianErr=apply(err.cv,1,median),SdErr=apply(err.cv,1,sd),SeErr=apply(err.cv,1,function(x) sd(x)/sqrt(length(x))))
    #-------------------------------------------------MeanErr
    p<-ggplot(ErrSumm,aes(x=NumVars,y=MeanErr))+
       geom_point(size=2)+geom_line(alpha=0.8)+
       xlab("Number of taxa")+ylab("Ten-fold CV error")+
       geom_errorbar(aes(ymin=MeanErr-SeErr,ymax=MeanErr+SeErr),width=0.05)+
       #coord_flip()+
       scale_x_continuous(trans = "log",breaks=breaks)+
       geom_hline(aes(yintercept=min(MeanErr)),linetype="longdash",alpha=0.2)
    ggsave(filename=paste(outpath,file.opts$prefix,".CV.MeanErrPlot.ggplot.pdf",sep=""),plot=p,width=6,height=4)
    #-------------------------------------------------MedianErr
    p<-ggplot(ErrSumm,aes(x=NumVars,y=MedianErr))+
       geom_point(size=2)+geom_line(alpha=0.8)+
       xlab("Number of taxa")+ylab("Ten-fold CV error")+
       geom_errorbar(aes(ymin=MedianErr-SeErr,ymax=MedianErr+SeErr),width=0.05)+
       scale_x_continuous(trans = "log",breaks=breaks)+
       geom_hline(aes(yintercept=min(MeanErr)),linetype="longdash",alpha=0.2)
    ggsave(filename=paste(outpath,file.opts$prefix,".CV.MedianErrPlot.ggplot.pdf",sep=""),plot=p,width=6,height=4)
    
    #-------------------------------------------------
    minErr_n_imps<-names(which.min(rowMeans(err.cv)))
    imps.minErr<-imps[1:minErr_n_imps]
    #-------------------------------------------------
    sink(paste(outpath,file.opts$prefix,".",file.opts$group.type,".imps_minErr_sorted.xls",sep=""));cat("\t");write.table(imps.minErr,quote=FALSE,sep="\t");sink()
    #-------------------------------------------------
    # Barplot of importances of all important taxa that minimized the error rate of model
    #-------------------------------------------------
    pdf(paste(outpath,file.opts$prefix,".impVars_minErr.ordered.",length(imps.minErr),".pdf",sep=""), height=length(imps.minErr)/4, width=6)
    print(barchart(imps.minErr,xlab="Mean decrease in accuracy",box.ratio=4))
    dev.off()
}

#--------------------------------------------------
#   BetweenGroup test
#--------------------------------------------------
    test.results<-BetweenGroup.test(x,y,p.adj.method = "fdr")
    if(nlevels(group)==2){
    sink(paste(outpath,file.opts$prefix,".",file.opts$group.type,".imps_wilcox.",nrow(test.results),".xls",sep=""))}else{
    sink(paste(outpath,file.opts$prefix,".",file.opts$group.type,".imps_kruscal.",nrow(test.results),".xls",sep=""))}
    cat("\t");write.table(test.results,quote=FALSE,sep="\t");sink()
    
#--------------------------------------------------
if(nlevels(group)==2){
#--------------------------------------------------
    fdr_cutoff<-file.opts$fdr_cutoff
    meanTab<-apply(x,2,function(x) tapply(x,y,mean))
    log.meanTab<-data.matrix(log.mat(meanTab[1,]/meanTab[2,])); AllEnr<-as.factor(log.meanTab>0)
    if(sum(AllEnr==TRUE)>0 & sum(AllEnr==FALSE)>0){
    AllEnr<-revalue(AllEnr,c("FALSE"="Decreased","TRUE"="Increased")) }else if(sum(AllEnr==TRUE)>0 & sum(AllEnr==FALSE)==0){
    AllEnr<-rep("Increased",length(AllEnr))
    }else{
    AllEnr<-rep("Decreased",length(AllEnr))
    }
    names(AllEnr)<-colnames(meanTab)
    All.enr<-data.frame(Log.ratio=log.meanTab,Enriched=AllEnr)
    #--------------------------------------------------
    IfSig<-as.factor(test.results$Wilcoxon.test_fdr<fdr_cutoff)
    if(sum(IfSig==TRUE)>0 & sum(IfSig==FALSE)>0){
    IfSig<-revalue(IfSig,c("FALSE"="NotSig","TRUE"="Sig")) }else if(sum(IfSig==TRUE)>0 & sum(IfSig==FALSE)==0){
    IfSig<-rep("Sig",length(IfSig))
    }else{
    IfSig<-rep("NotSig",length(IfSig))
    }
    names(IfSig)<-rownames(test.results)
    SigEnr<-with(data.frame(IfSig,AllEnr), interaction(IfSig,AllEnr))
    SigEnr<-gsub("NotSig.*","Neutral",SigEnr,perl=TRUE)
    SigEnr<-as.factor(gsub("Sig.","",SigEnr))
    All.Wilcox.enr<-data.frame(Log.ratio=log.meanTab,Enriched=SigEnr)
    #--------------------------------------------------
     
    #--------------------------------------------------
    wilcoxEnr<-All.Wilcox.enr[rownames(test.results[which(test.results$Wilcoxon.test_fdr<fdr_cutoff),]),]
    gmat.wilcox<-gmat[,rownames(test.results[which(test.results$Wilcoxon.test_fdr<fdr_cutoff),])]
    sink(paste(outpath,file.opts$prefix,".",file.opts$group.type,".imps_wilcoxSig.",nrow(wilcoxEnr),".xls",sep=""));cat("\t");write.table(wilcoxEnr,quote=FALSE,sep="\t");sink()
    sink(paste(outpath,file.opts$prefix,".",file.opts$group.type,".wilcoxSig.",ncol(gmat.wilcox),".xls",sep=""));cat("\t");write.table(gmat.wilcox,quote=FALSE,sep="\t");sink()
    #--------------------------------------------------BetweenGroup Plot
    Log10Meds<-rbind(apply(log.mat(x,base=10),2,function(x) tapply(x,group,median)),All=apply(log.mat(x,base=10),2,median),Enr=as.factor(test.results$Wilcoxon.test_fdr<fdr_cutoff))
    pdf(paste(outpath,file.opts$prefix,".",file.opts$pred.group.type,".",levels(group)[1],"_VS_",levels(group)[2],".Wilcox.",nrow(wilcoxEnr),".LogMedAbundPlot.pdf",sep=""), height=5, width=5)
        par(mar=c(4,4,3,3))
    plot(Log10Meds[1,],Log10Meds[2,],pch=21,cex=(-1/Log10Meds[3,]*2), col="white",bg=c("grey","red")[unclass=as.factor(Log10Meds["Enr",])],bty='n',
         xlab=paste( levels(group)[1]," median relative abundance log10"),
         ylab=paste( levels(group)[2]," median relative abundance log10"))
    abline(0,1, lty=3,col="grey")
    dev.off()

#-------------------------------
# AUC of individual features in train data
#-------------------------------
    AUCs<-apply(x,2,function(x) auc(y,x));AUCs[AUCs<0.5]<-1-AUCs[AUCs<0.5]
    pdf(paste(outpath,file.opts$prefix,".AUC_distribution.",ncol(x),".pdf",sep=""), height=4, width=4)
    hist(AUCs,xlim=c(0.5,1),xlab="AUC",main="AUC distribution of all features",breaks=10,col="grey60")
    dev.off()
    #-------------------------------
    AUCs=AUCs[order(names(AUCs))]
    ImpScore=imps[order(names(imps))]
    MeanAbund<-apply(x,2,mean);MeanAbund<-MeanAbund[order(names(MeanAbund))]
    MedianAbund<-apply(x,2,median);MedianAbund<-MedianAbund[order(names(MedianAbund))]
    TaxaSumm<-data.frame(AUCs,ImpScore,MeanAbund,MedianAbund,All.Wilcox.enr[order(rownames(All.Wilcox.enr)),])
    TaxaSumm$TaxaNames<-rownames(TaxaSumm)
    #-------------------------------
    p<-ggplot(TaxaSumm,aes(x=AUCs,y=ImpScore))+xlab("AUC")+ylab("Importance score")+xlim(0.5,1)+
       geom_point(aes(colour=Enriched,size=MeanAbund))+
       geom_vline(xintercept=auc,linetype = "longdash")+
       geom_text(data=subset(TaxaSumm,Enriched!="Neutral"),aes(x=AUCs,y=ImpScore,label=TaxaNames),size=2,angle=0,hjust=0, vjust=0,alpha = 0.25)
    ggsave(filename=paste(outpath,file.opts$prefix,".AUC_VS_ImpScore.",ncol(x),".ColByWilcox",nrow(wilcoxEnr),".pdf",sep=""),plot=p,width=6,height=4)
    
#-------------------------------
# Occrruence rate of individual features in train data
#-------------------------------
    AllOccRate=t(apply(x,2,function(x) tapply(x,group,function(x) sum(x!=0)/length(x)))); colnames(AllOccRate)<-paste("OccRate",colnames(AllOccRate),sep="_");AllOccRate<-AllOccRate[order(rownames(AllOccRate)),]
    AllMedianAbund=t(apply(x,2,function(x) tapply(x,group,median))); colnames(AllMedianAbund)<-paste("MedianAbund",colnames(AllMedianAbund),sep="_");AllMedianAbund<-AllMedianAbund[order(rownames(AllMedianAbund)),]
    AllMedianAbundLog10<-log.mat(AllMedianAbund,base=10); colnames(AllMedianAbundLog10)<-paste("Log10",colnames(AllMedianAbundLog10),sep="")
    TaxaSumm<-data.frame(TaxaSumm,AllOccRate,AllMedianAbund,AllMedianAbundLog10)
    TaxaColor<-c(rgb(158,11,11,maxColorValue=255),rgb(33,99,48,maxColorValue=255),"grey60")
    #-------------------------------
    if("Increased" %in% levels(All.Wilcox.enr$Enriched) & summary(All.Wilcox.enr$Enriched)["Increased"]>1){
    TaxaSumm.Inc<-subset(TaxaSumm,Enriched=="Increased")[,c(7,8,12)]; colnames(TaxaSumm.Inc)<-c("TaxaNames","OccRate","Log10MedianAbund")
    p<-ggplot(data=TaxaSumm.Inc,aes(x=OccRate,y=Log10MedianAbund))+xlim(0,1)+
       geom_point(colour=TaxaColor[1],alpha=0.6)+
       geom_text(data=TaxaSumm.Inc,aes(x=OccRate,y=Log10MedianAbund,label=TaxaNames),size=2,angle=0,hjust=1, vjust=0,alpha = 0.25)+
       xlab(paste("Occurrence Rate in ",levels(group)[1]," ( n=",summary(group)[1]," )",sep=""))+
       ylab(paste("Median abundance in ",levels(group)[1]," (log10 scale)",sep=""))     
    ggsave(filename=paste(outpath,file.opts$prefix,".OccRate_VS_MedAbund.wilcox_",nrow(TaxaSumm.Inc),".",levels(group)[1],"-",summary(group)[1],".pdf",sep=""),plot=p,width=6,height=6)
    }
    if("Decreased" %in% levels(All.Wilcox.enr$Enriched) & summary(All.Wilcox.enr$Enriched)["Decreased"]>1){
    TaxaSumm.Dec<-subset(TaxaSumm,Enriched=="Decreased")[,c(7,8,12)]; colnames(TaxaSumm.Dec)<-c("TaxaNames","OccRate","Log10MedianAbund")
    p<-ggplot(data=TaxaSumm.Dec,aes(x=OccRate,y=Log10MedianAbund))+xlim(0,1)+
       geom_point(colour=TaxaColor[2],alpha=0.6)+
       geom_text(data=TaxaSumm.Dec,aes(x=OccRate,y=Log10MedianAbund,label=TaxaNames),size=2,angle=0,hjust=1, vjust=0,alpha = 0.25)+
       xlab(paste("Occurrence Rate in ",levels(group)[2]," ( n=",summary(group)[2]," )",sep=""))+
       ylab(paste("Median abundance in ",levels(group)[2]," (log10 scale)",sep=""))     
    ggsave(filename=paste(outpath,file.opts$prefix,".OccRate_VS_MedAbund.wilcox_",nrow(TaxaSumm.Dec),".",levels(group)[2],"-",summary(group)[2],".pdf",sep=""),plot=p,width=6,height=6)
    }
    #-------------------------------
    sink(paste(outpath,file.opts$prefix,".",file.opts$group.type,".All_AUCs-Imps-OccRates.",ncol(x),".xls",sep=""))
    cat("\t");write.table(TaxaSumm,quote=FALSE,sep="\t");sink()
#-------------------------------
# Barplot of importances of all important taxa that minimized the error rate of model
#-------------------------------
    if(file.opts$rf.cv){
    minErr_n_imps<-names(which.min(rowMeans(err.cv)))
    TaxaSumm.ordered<-TaxaSumm[order(TaxaSumm$ImpScore,decreasing=T),]
    TaxaSumm.minErr<-TaxaSumm.ordered[1:minErr_n_imps,]
    sink(paste(outpath,file.opts$prefix,".",file.opts$group.type,".RFminErr_AUCs-Imps.",nrow(TaxaSumm.minErr),".xls",sep=""));cat("\t");write.table(TaxaSumm.minErr,quote=FALSE,sep="\t");sink()
    #-------------------------------
    byimps <- with(TaxaSumm.minErr, reorder(TaxaNames, ImpScore))
    TaxaSumm.minErr$TaxaNames<-factor(TaxaSumm.minErr$TaxaNames,levels=rev(byimps),ordered=TRUE)
    plot<- ggplot(TaxaSumm.minErr, aes(x=TaxaNames, y=ImpScore)) + ylab('Mean decrease in accuracy')+
           geom_bar(stat="identity",fill = I("grey50")) +
           coord_flip()
    ggsave(filename=paste(outpath,file.opts$prefix,".RFminErr_impVars.",nrow(TaxaSumm.minErr),".ggplot.pdf",sep=""),plot=plot,limitsize=FALSE,width=6) #height=nrow(imps.minErr)/2,
    }
#-------------------------------
# Derived index based important features
#-------------------------------
    Index<-function(mat,A=Increased,B=Decreased){
           #-------------------------------Data check
           if(!all(is.na(A))){
           Mis_A<-A[which(!A %in% colnames(mat))]
           if(length(Mis_A)>0){
           cat(length(Mis_A),' specified "Incresed" varible(s) NOT found in data: ',Mis_A,'\n',sep=' ')
           A<-A[which(A!=Mis_A)]}
           }
           if(!all(is.na(B))){
           Mis_B<-B[which(!B %in% colnames(mat))]
           if(length(Mis_B)>0){
           cat(length(Mis_B),' specified "Decresed" varible(s) NOT found in data: ',Mis_B,'\n',sep=' ')
           B<-B[which(B!=Mis_B)]}
           }
           #-------------------------------Index calculation
           if(!all(is.na(A)) & !all(is.na(B))){
           Index<-apply(mat,1,function(x){(sum(x[A])/length(A)-sum(x[B])/length(B))*10})}else 
           if(all(is.na(A)) & !all(is.na(B))){Index<-apply(mat,1,function(x){sum(x[B])/length(B)*10})
           }else{Index<-apply(mat,1,function(x){sum(x[A])/length(A)*10})
           }
           
           return(Index)
    }
    #--------------------------------------------------
    #  Microbial Index derived based on Wilcoxon taxa
    #--------------------------------------------------
    if(nrow(wilcoxEnr)>0){
    if(sum((wilcoxEnr$Enriched)=="Increased")){Increased=rownames(wilcoxEnr)[which(wilcoxEnr$Enriched=="Increased")]}else{Increased=NA}
    if(sum((wilcoxEnr$Enriched)=="Decreased")){Decreased=rownames(wilcoxEnr)[which(wilcoxEnr$Enriched=="Decreased")]}else{Decreased=NA}
    #-------------------------------Microbial Index Calculation in Train data set
    IndexInTrain<-Index(gmat,Increased,Decreased)
    #--------------------------------------------------
    #  ROC plot of Microbial Index using "pROC" package
    #--------------------------------------------------
    pdf(paste(outpath,file.opts$prefix,".",file.opts$group.type,".AllWilcox.",nrow(wilcoxEnr),".ROC.ci.index.pdf",sep=""),width=5,height=5)
    rocobj <- plot.roc(y, IndexInTrain,main="Microbial Index", percent=TRUE,ci=TRUE) # print the AUC (will contain the CI)
    ciobj <- ci.se(rocobj,specificities=seq(0, 100, 5)) # over a select set of specificities
    plot(ciobj, type="shape", col="#1c61b6AA", print.auc=TRUE) # plot as a blue shape
    text(50,20,paste("AUC = ",round(rocobj$auc,2),sep=""),pos=4)
    ci.lower<-formatC(rocobj$ci[1],digits=2,format="f")
    ci.upper<-formatC(rocobj$ci[3],digits=2,format="f")
    text(50,10,paste("95% CI: ",ci.lower,"-",ci.upper,sep=""),pos=4)
    dev.off()
    #--------------------------------------------------
    #  Box plot of Microbial Index
    #--------------------------------------------------
    if(!is.na(file.opts$pred.prefix)){
    IndexInTest<-Index(pred.gmat,Increased,Decreased)
    if(file.opts$pred.group.type==file.opts$pred.group.type){
        comb.group<-as.factor(c(as.character(group),paste("relative",as.character(pred.group),sep="")))
        comb.group<-factor(comb.group,levels=c(levels(comb.group)[2],levels(comb.group)[3],levels(comb.group)[1]))
        comb.Index<-c(IndexInTrain,IndexInTest)
        pdf(paste(outpath,file.opts$prefix,"-",file.opts$pred.prefix,".",file.opts$group.type,".AllWilcox.",nrow(wilcoxEnr),".Index_boxplot.pdf",sep=""), height=8, width=4)
        par(mar=c(4,4,3,3))
        boxplot(comb.Index~comb.group, ylab=paste("Microbial Index of ",colnames(comb.prob)[1],sep=""),outline=FALSE)
        abline(h=0.5,lty=2)
        stripchart(comb.Index~comb.group,method="jitter",jitter=.05,vertical=T,add=T,pch=20) 
        dev.off()
        sink(paste(outpath,file.opts$prefix,"-",file.opts$pred.prefix,".",file.opts$group.type,".AllWilcox.",nrow(wilcoxEnr),".Index.xls",sep="")); cat("\t"); write.table(comb.Index,sep="\t",quote=FALSE); sink()
        }else{
        pdf(paste(outpath,file.opts$prefix,".",file.opts$group.type,".AllWilcox.",nrow(wilcoxEnr),".Index_boxplot.pdf",sep=""), height=8, width=4)
        par(mar=c(4,4,3,3))
        boxplot(IndexInTrain~group, ylab=paste("Microbial Index of ",colnames(predicted.prob)[1],sep=""),outline=FALSE)
        abline(h=0.5,lty=2)
        stripchart(IndexInTrain~group,method="jitter",jitter=.05,vertical=T,add=T,pch=20) 
        dev.off()
        sink(paste(outpath,file.opts$prefix,"-",file.opts$pred.prefix,".",file.opts$group.type,".AllWilcox.",nrow(wilcoxEnr),".Index.xls",sep="")); cat("\t"); write.table(IndexInTrain,sep="\t",quote=FALSE); sink()
        }
        
    }
    

#-------------------------------
# Barplot of importances of all important taxa that picked by wilcox test using fdr cutoff
#-------------------------------
    TaxaSumm.wilcox<-TaxaSumm[which(rownames(TaxaSumm) %in% rownames(wilcoxEnr)),]
    sink(paste(outpath,file.opts$prefix,".",file.opts$group.type,".wilcox_AUCs-Imps.",nrow(TaxaSumm.wilcox),".xls",sep=""));cat("\t");write.table(TaxaSumm.wilcox,quote=FALSE,sep="\t");sink()
    #-------------------------------
    byimps <- with(TaxaSumm.wilcox, reorder(TaxaNames, ImpScore))
    TaxaSumm.wilcox$TaxaNames<-factor(TaxaSumm.wilcox$TaxaNames,levels=rev(byimps),ordered=TRUE)
    plot<- ggplot(TaxaSumm.wilcox, aes(x=TaxaNames, y=ImpScore)) + ylab('Mean decrease in accuracy')+
           geom_bar(stat="identity",fill = I("grey50")) +
           coord_flip()
    ggsave(filename=paste(outpath,file.opts$prefix,".wilcox_impVars.",nrow(TaxaSumm.wilcox),".ggplot.pdf",sep=""),plot=plot,limitsize=FALSE,width=6) #height=nrow(imps.minErr)/2,
    
#-------------------------------
# Heatmap of important taxa in train data
#-------------------------------
    log.gmat<-log.mat(gmat[,rownames(wilcoxEnr)],base=10)
    log.gmat<-log.gmat[order(group),order(colnames(log.gmat))]
    #-------------------------------
    annotation = data.frame(Group = as.factor(group[order(group)]))
    colnames(annotation) = "Group"
    rownames(annotation) = rownames(log.gmat)
    Group = c(rgb(107,73,242, max = 255), rgb(208,242,73, max = 255))
    names(Group) = levels(group)
    ann_colors = list(Group = Group)
    h<-pheatmap(t(log.gmat),annotation =annotation,annotation_colors = ann_colors, cluster_row=TRUE,cluster_col=F,fontsize_row=4,fontsize_col=1,cellwidth=1,cellheight=4,filename=paste(outpath,file.opts$prefix,".impVars.wilcox.",nrow(wilcoxEnr),".heatmap.log10.pdf",sep=""))

#-------------------------------
# Boxplot of abundance of important taxa in train data
#-------------------------------
    MyColor<-c(rgb(158,11,11,maxColorValue=255),rgb(33,99,48,maxColorValue=255))
    log.mat_Group<-data.frame(Group=group,log.gmat[order(rownames(log.gmat)),order(colnames(log.gmat))])
    data_melt<-melt(log.mat_Group)
    p <- ggplot(data=data_melt,aes(x=variable,y=value,fill=Group))+xlab("Taxa")+ylab("Relative abundance (log10 scale)")+
         geom_boxplot(outlier.shape = NA,aes(fill=Group))+
         scale_fill_manual(name = file.opts$group.type, values = MyColor)+
         coord_flip()+#if want to filp coordinate
         geom_jitter(position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,aes(fill=Group),alpha=0.4)#jitter???????ڻ???
    ggsave(filename=paste(outpath,file.opts$prefix,".impVars.wilcox.",nrow(wilcoxEnr),".boxplot.pdf",sep=""),plot=p,limitsize = FALSE,width=6)
    
#-------------------------------
# Scatterplot of AUC of important taxa in train data
#-------------------------------
    wilcox.TaxaSumm<-TaxaSumm[colnames(log.gmat),]
    EnrColor<-c(rgb(33,99,48,maxColorValue=255),rgb(158,11,11,maxColorValue=255))
    p <- ggplot(data=wilcox.TaxaSumm,aes(x=AUCs,y=TaxaNames))+geom_point(aes(colour=Enriched))+xlim(0.5,1)+
         geom_vline(xintercept=auc,linetype = "longdash")+
         scale_color_manual(name = "Enrichment", values = EnrColor)
    ggsave(filename=paste(outpath,file.opts$prefix,".impVars.wilcox.",nrow(wilcoxEnr),".AUC.scatterplot.pdf",sep=""),plot=p,limitsize = FALSE,width=7)

    if(!is.na(file.opts$pred.prefix)){

#-------------------------------
# Heatmap of important taxa in test data
#-------------------------------
    log.pred.gmat<-log.mat(pred.gmat[,rownames(wilcoxEnr)],base=10)
    log.pred.gmat<-log.pred.gmat[,order(colnames(log.pred.gmat))]
    log.pred.gmat<-log.pred.gmat[order(pred.group),h$tree_row$order]
    #-------------------------------
    #-------------------------------
    annotation = data.frame(Var1 = predicted.prob[,1])
    colnames(annotation) = "Predicted"
    rownames(annotation) =  rownames(log.pred.gmat)
    ann_colors = list(Predicted = c(rgb(208,242,73, max = 255),rgb(107,73,242, max = 255)))
    pheatmap(t(log.pred.gmat),annotation =annotation,annotation_colors = ann_colors, cluster_row=FALSE,cluster_col=TRUE,fontsize_row=4,fontsize_col=1,cellwidth=1,cellheight=4,filename=paste(outpath,file.opts$pred.prefix,".impVars.wilcox.",nrow(wilcoxEnr),".heatmap.log10.pdf",sep=""))
    
#-------------------------------
# Boxplot of abundance of important taxa in test data
#-------------------------------
    log.pred.mat_Group<-data.frame(Group=pred.group,log.pred.gmat[order(rownames(log.pred.gmat)),order(colnames(log.pred.gmat))])
    data_melt<-melt(log.pred.mat_Group)
    p <- ggplot(data=data_melt,aes(x=variable,y=value,fill=Group))+xlab("Taxa")+ylab("Relative abundance (log10 scale)")+
         geom_boxplot(outlier.shape = NA,aes(fill=Group))+
         scale_fill_manual(name = file.opts$group.type, values = MyColor)+
         coord_flip()+#if want to filp coordinate
         geom_jitter(position=position_jitterdodge(jitter.width= 0.2,dodge.width = 0.8),size=1,aes(fill=Group),alpha=0.4)#jitter???????ڻ???
    ggsave(filename=paste(outpath,file.opts$pred.prefix,".impVars.wilcox.",nrow(wilcoxEnr),".boxplot.pdf",sep=""),plot=p,limitsize = FALSE,width=6)
    }
}    
} 


##-------------------------------
#sink(type="message")
#sink() 
##-------------------------------

