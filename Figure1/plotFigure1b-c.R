library(SingleCellExperiment)
library(gridExtra)
library(grid)
library(foreach)

source("../Functions/QQsim_v2_SingleCellExperiment.R")
source("../Functions/DROPOUT_FUN.r")

sm <- function(x, y, x.log = FALSE,n.bins = 25){
    if(x.log){ 
        brks <- unique(quantile(x, probs = seq(0,1,len=25))) 
    } else {
        brks <- 2^unique(quantile(log2(x), probs = seq(0,1,len=n.bins))) 
    }
    mids <- (brks[-1] + brks[-length(brks)] )/ 2
    x.in <- cut(x, breaks = brks, include.lowest = TRUE)
    m <- tapply(y, x.in, mean)
    fit = lm(y~x)
    l <- predict(fit, newdata = data.frame(x  = mids))
    dat <- data.frame(x=mids, y = m, n = as.numeric(table(x.in)), pred = l)
}


dropoutfun<-function(data){
    xx<-rowMeans(data)
    yy<-apply(data,1,function(x){length(which(x==0))/length(x)})
    qq<-cbind(xx,yy)
    return(qq)
}
SIM_FUN<-function(DATA,MU,SIZE,BETA)
{
    
    nCells<-dim(DATA)[2]
    nGenes<-dim(DATA)[1]
    
    GeneMean_mat<-matrix(MU,ncol=nCells ,nrow=nGenes,byrow=F)
    
    one_bcv2<-SIZE
    
    Gamma_Means_mat <- matrix(rgamma(nGenes * nCells, shape =one_bcv2, scale = GeneMean_mat * (1/one_bcv2)),nrow = nGenes, ncol = nCells)
    
    true.counts <- matrix(rpois(nGenes * nCells, lambda =  Gamma_Means_mat ),nrow = nGenes, ncol = nCells)
    rownames(true.counts)<-rownames(DATA)
    colnames(true.counts)<-colnames(DATA)
    
    downsample.counts <-bayNorm::DownSampling(true.counts,BETA)
    
    rownames(downsample.counts)<-rownames(DATA)
    colnames(downsample.counts)<-colnames(DATA)
    
    return(list(true.counts=true.counts,downsample.counts=downsample.counts))
}


if(!file.exists("Simulations_realdata.RData")) {

	load("../RealData/Klein_study/bay_Klein.RData")

	#run "Binomial_bayNorm" simulation protocol
	BaySim_Kelin<-SIM_FUN(DATA=Real_Klein,MU=bay_Klein$PRIORS$MME_prior$MME_MU,SIZE=bay_Klein$PRIORS$MME_SIZE_adjust,BETA=bay_Klein$BETA)

	#run "Binomial" simulation protocol
	Sim_List_Input<-trytt2(counts=Real_Klein,inputBeta = NULL,inputMeanBeta = 0.06,inputtrim = 0.01)

	#run Splatter
	library(splatter)
	# argument to as.matrix was originally 'Real_data', replacing with 'Real_Klein'
	splatter_klein_params <- splatEstimate(as.matrix(Real_Klein))
	splatter_klein_sim<-splatSimulate(splatter_klein_params)

	#Prepare SCE lists, then we will put it into a function for producing comparison plots
	SCElist_Klein2<-list(
	    Real=SingleCellExperiment(assays=list(counts=as.matrix(Real_Klein))),
	    Binomial_bayNorm=SingleCellExperiment(assays=list(counts=BaySim_Kelin$downsample.counts)),
	    Binomial=SingleCellExperiment(assays=list(counts=assays(Sim_List_Input$SCE)$counts)),
	    Splatter=splatter_klein_sim)

	save.image("Simulations_realdata.RData")
} else {
	load("Simulations_realdata.RData")
}

point.size<-1
point.alpha<-0.4
linewidth<-1
linewidth.exp<-1.5

textsize<-10
legendpointsize=4
legend_key_size=0.8

theme_set(theme_grey()+theme(
  axis.title = element_text(size=12),
  plot.title=element_text(hjust=0.5, size=14)))

sces<-SCElist_Klein2
    
for (name in names(sces)) {
sce <- sces[[name]]
rowData(sce)$Dataset <- name
colData(sce)$Dataset <- name
sce <- scater::calculateQCMetrics(sce)
cpm(sce) <- scater::calculateCPM(sce, use_size_factors = FALSE)
sce <- addFeatureStats(sce, "counts")
sce <- addFeatureStats(sce, "cpm")
sce <- addFeatureStats(sce, "cpm", log = TRUE)
n.features <- colData(sce)$total_features_by_counts
colData(sce)$PctZero <- 100 * (1 - n.features / nrow(sce))
sces[[name]] <- sce
}

features <- rowData(sces[[1]])
cells <- colData(sces[[1]])

if (length(sces) > 1) {
for (name in names(sces)[-1]) {
    sce <- sces[[name]]
    features <- rbindMatched(features, rowData(sce))
    cells <- rbindMatched(cells, colData(sce))
}
}
features$Dataset <- factor(features$Dataset, levels = names(sces))
cells$Dataset <- factor(cells$Dataset, levels = names(sces))
features <- data.frame(features)
cells <- data.frame(cells)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colours <- cbbPalette[seq(1,length(sces))]

mean.var <- #ggplot()+
ggplot(features,aes_string(x = "MeanLogCPM", y = "VarLogCPM",colour = "Dataset", fill = "Dataset")) +
geom_point(size = point.size, alpha = point.alpha,shape=46) +
scale_colour_manual(values = colours) +
scale_fill_manual(values = colours) +
guides(colour = guide_legend(override.aes = list(size=legendpointsize,shape=16,alpha=1)))+
xlab("Mean expression") +
ylab("Variance of gene expression") +
ggtitle("Mean-variance relationship")+
theme(legend.position=c(0.75,0.8))


theoline_ref_name<-names(sces)[1]
theoline<-data.frame(xx =features$mean_counts[which(features$Dataset== theoline_ref_name)], yy=exp(-features$mean_counts[which(features$Dataset== theoline_ref_name)]))

mean.zeros <- ggplot() + geom_point(data=features,
   aes_string(x = "MeanCounts",
	      y = features$pct_dropout_by_counts/100,colour = "Dataset",
	      fill = "Dataset"),
   size = point.size,
   alpha = point.alpha,shape=46,show.legend=F) +
geom_line(aes(x=theoline$xx,y=theoline$yy, linetype = 'exp(-mean expression)'),size = linewidth.exp,alpha = point.alpha)+
scale_x_log10(labels = scales::comma) +
scale_colour_manual(values = colours) +
scale_linetype_manual(values='dashed','Dashed line')+
labs(x = "Mean expression",y="Dropout rate",fill='Data',colour='Dataset',title="Mean-dropout rate relationship",subtitle='')+
theme(legend.position=c(0.75,0.9))
grob<-cbind(ggplotGrob(mean.var), ggplotGrob(mean.zeros))
ggsave(grob, file="figure_1b-c.pdf", width=10, height=5)
