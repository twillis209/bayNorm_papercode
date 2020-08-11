source("../Functions/QQsim_v2_SingleCellExperiment.R")
source("../Functions/DROPOUT_FUN.r")
library(SingleCellExperiment)
library(ggplot2)

#######Bacher study (H1_P24 cells)####
load("../RealData/Bacher_study/H1_many_normalizations.RData")

theme_set(theme_grey()+theme(
  axis.title = element_text(size=12),
  plot.title=element_text(hjust=0.5, size=14)))


#Run Binomial simulation protocol
qq<-round(H1_p24/20)
system.time(H1p24_bay_sim<-trytt2(counts=qq,inputMeanBeta=1,inputtrim=0.01,inputBeta=Beta_H1[colnames(qq)]))

H1p24_real<-SingleCellExperiment(assays=list(counts=as.matrix(H1_p24)))
#H1p24_bay_DAT<-SingleCellExperiment(assays=list(counts=H1p24_bay_sim$SCE@assays@.xData$data$counts))
H1p24_bay_DAT<-SingleCellExperiment(assays=list(counts=assays(H1p24_bay_sim$SCE)$counts))
# No reference to H1p24_splatter_sim elsewhere in the repo
#H1p24_splatter_DAT<-SingleCellExperiment(assays=list(counts=assays(H1p24_splatter_sim$SCE)$counts))

#prepare list of data for further analysis
listused_H1p24<-list(Real = assays(H1p24_real)$counts,Binomial =assays(H1p24_bay_DAT)$counts)



#######Bacher study (H1_P96 cells)####
# Seems we only need the H1_p96 object
#load("H9_many_normalizations.RData")
load("../RealData/Bacher_study/RAW_INITIATE.RData")

#Scaled Non-UMI data
qq2<-round(H1_p96/10)

#Run Binomial simulation protocol
H1p96_bay_sim<-trytt2(counts=qq2,inputMeanBeta=1,inputtrim=0.01,inputBeta=Beta_H1[colnames(qq2)])

H1p96_real<-SingleCellExperiment(assays=list(counts=as.matrix(H1_p96)))
H1p96_bay_DAT<-SingleCellExperiment(assays=list(counts=assays(H1p96_bay_sim$SCE)$counts))

#prepare list of data for further analysis
listused_H1p96<-list(Real = assays(H1p96_real)$counts,Binomial =assays(H1p96_bay_DAT)$counts)

########Islam study (H1_P96 cells)####
load('../RealData/Islam_study/Islam_many_normalizations.RData')

#scaled non-UMI data
qq<-round(DAT_Jaakkola/10)

#Run Binomial simulation protocol
Islam_bay_sim<-trytt2(counts=qq,inputMeanBeta=1,inputtrim=0.01,inputBeta=BETA)

Islam_real<-SingleCellExperiment(assays=list(counts=as.matrix(DAT_Jaakkola)))
Islam_bay_DAT<-SingleCellExperiment(assays=list(counts=assays(Islam_bay_sim$SCE)$counts))

listused_Islam<-list(Real = assays(Islam_real)$counts,Binomial =assays(Islam_bay_DAT)$counts)

#begin plotting Fig S9####
names(listused_H1p24)<-c('Real','Binomial')
names(listused_H1p96)<-c('Real','Binomial')

point.size<-1
point.alpha<-0.4
linewidth<-1
linewidth.exp<-1.5

textsize<-10
legendpointsize=4
legend_key_size=0.8


plot_DROPOUT<-function(listused_N1,MAIN='',CAPTION='',legendpointsize=1,legend_key_size=1,subtitle=''){

library(foreach)
DROPOUT_DAT<-foreach(i=1:length(listused_N1),.combine=rbind)%do%{

    if(names(listused_N1)[i]=='Real') {colll<-cbPaletteee[1]}else if(names(listused_N1)[i]=='Binomial'){
        colll<-cbPaletteee[2]
    } else if(names(listused_N1)[i]=='Scaled raw'){
        colll<-cbPaletteee[3]
    } else{colll<-cbPaletteee[1]}

  dropout<-apply(listused_N1[[i]],1,function(x){length(which(x==0))/length(x)})
  meann<-rowMeans(listused_N1[[i]])
  qqq<-cbind(dropout,meann,rep(names(listused_N1)[i],length(dropout)),rep(colll,length(dropout)))
  return(qqq)
}
DROPOUT_DAT<-as.data.frame(DROPOUT_DAT)
colnames(DROPOUT_DAT)<-c('Dropout rate','Mean expression','Dataset','Colour')
DROPOUT_DAT$`Dropout rate`<-as.numeric(as.character(DROPOUT_DAT$`Dropout rate`))
DROPOUT_DAT$`Mean expression`<-as.numeric(as.character(DROPOUT_DAT$`Mean expression`))
DROPOUT_DAT$Dataset<-factor(DROPOUT_DAT$Dataset,levels=unique(DROPOUT_DAT$Dataset))
DROPOUT_DAT$Colour<-factor(DROPOUT_DAT$Colour,levels=unique(DROPOUT_DAT$Colour))

sces<-listused_N1
colours <- scales::hue_pal()(length(sces))

theoline_ref_name<-names(listused_N1)[1]

theoline<-data.frame(xx =sort(DROPOUT_DAT$`Mean expression`[which(DROPOUT_DAT$Dataset==theoline_ref_name)]), yy=exp(-sort(DROPOUT_DAT$`Mean expression`[which(DROPOUT_DAT$Dataset==theoline_ref_name)])))

mean.zeros <- ggplot() +
    geom_line(data=theoline,aes(x=xx,y=yy, lty = 'exp(-mean expression)'))+
    geom_point(data=DROPOUT_DAT,aes_string(x = DROPOUT_DAT$`Mean expression`, y = DROPOUT_DAT$`Dropout rate`,colour = DROPOUT_DAT$Colour),size = point.size, alpha = point.alpha,shape=46) +
  scale_x_log10(labels = scales::comma) +
  xlab("Mean expression") +
  ylab("Dropout rates") +
    scale_color_manual(values=as.character(unique(DROPOUT_DAT$Colour)),labels=names(listused_N1))+
    scale_linetype_manual(values=1,'Black line')+
  labs(x = "Mean expression",y="Dropout rates",fill='Dataset',colour='Dataset',caption=CAPTION,title=MAIN,subtitle=subtitle)+
    guides(colour = guide_legend(override.aes = list(size=legendpointsize,shape=16,alpha=1)))
return(mean.zeros)
}

D_SCnorm<-plot_DROPOUT(listused_H1p24,MAIN='',legendpointsize=legendpointsize,legend_key_size=legend_key_size,CAPTION='',subtitle='Black line: exp(-mean expression of raw data)')

D_SCnorm_96<-plot_DROPOUT(listused_H1p96,MAIN='',legendpointsize=legendpointsize,legend_key_size=legend_key_size,CAPTION='',subtitle='')

listused_H1p24_div<-listused_H1p24
names(listused_H1p24_div)<-c('Scaled raw','Binomial')
listused_H1p24_div$`Scaled raw`<-round(listused_H1p24$Real/20)
listused_H1p24_div$Binomial<-listused_H1p24$Binomial

D_SCnorm_div<-plot_DROPOUT(listused_H1p24_div,MAIN='',legendpointsize=legendpointsize,legend_key_size=legend_key_size,CAPTION='',subtitle='')

listused_H1p96_div<-listused_H1p96
names(listused_H1p96_div)<-c('Scaled raw','Binomial')
listused_H1p96_div$`Scaled raw`<-round(listused_H1p96$Real/10)
listused_H1p96_div$Binomial<-listused_H1p96$Binomial

D_SCnorm_div_96<-plot_DROPOUT(listused_H1p96_div,MAIN='',legendpointsize=legendpointsize,legend_key_size=legend_key_size,CAPTION='',subtitle='')

# TODO: D_Islam and D_Islam_div
D_Islam<-plot_DROPOUT(listused_Islam,MAIN='',legendpointsize=legendpointsize,legend_key_size=legend_key_size,CAPTION='',subtitle='')

listused_Islam_div<-listused_Islam
names(listused_Islam_div)<-c('Scaled raw','Binomial')
# I took division by 10 from the supplement
listused_Islam_div$`Scaled raw`<-round(listused_Islam$Real/10)
listused_Islam_div$Binomial<-listused_Islam$Binomial

D_Islam_div<-plot_DROPOUT(listused_Islam_div,MAIN='',legendpointsize=legendpointsize,legend_key_size=legend_key_size,CAPTION='',subtitle='')

cbPaletteee <- c("#999999", "#E69F00", "#56B4E9")

library(gridExtra)
library(ggpubr)
library(cowplot)

#grid.arrange(D_SCnorm,D_SCnorm_div,D_SCnorm_96,D_SCnorm_div_96,D_Islam,D_Islam_div,nrow=3,ncol=2)
pp<-plot_grid(
    D_SCnorm+ theme(legend.position="none")+ggtitle("Bacher data: H1-4M hESCs"),
    D_SCnorm_div+ theme(legend.position="none")+ggtitle("Bacher data: H1-4M hESCs"),
    D_SCnorm_96+ theme(legend.position="none")+ggtitle("Bacher data: H1-1M hESCs"),
    D_SCnorm_div_96+ theme(legend.position="none")+ggtitle("Bacher data: H1-1M hESCs"),
    D_Islam+ theme(legend.position="none")+ggtitle("Islam data"),
    D_Islam_div+ theme(legend.position="none")+ggtitle("Islam data"),
    nrow=3,ncol=2)
mm<-ggplot()+geom_point(aes(x=seq(1,3),y=seq(1,3),color=cbPaletteee))+scale_color_manual(values=cbPaletteee ,labels=c('Real','Binomial','Scaled and rounded real data'))+labs(colour='Data set')+theme(legend.position='bottom')

pp<-plot_grid( pp, get_legend(mm),axis='b',align='v',ncol=1,rel_heights=c(20,1))
ggsave(pp, filename="myFigure_S9.pdf")
