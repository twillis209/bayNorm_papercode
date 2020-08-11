source("../Functions/QQsim_v2_SingleCellExperiment.R")
source("../Functions/DROPOUT_FUN.r")
library(SingleCellExperiment)

#######Bacher study (H1_P24 cells)####
load("../RealData/Bacher_study/H1_many_normalizations.RData")


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
library(cowplot)
grid.arrange(D_SCnorm,D_SCnorm_div,D_SCnorm_96,D_SCnorm_div_96,D_Islam,D_Islam_div,nrow=3,ncol=2)
pp<-plot_grid(
    D_SCnorm+ theme(legend.position="none"),
    D_SCnorm_div+ theme(legend.position="none"),
    D_SCnorm_96+ theme(legend.position="none"),
    D_SCnorm_div_96+ theme(legend.position="none"),
    D_Islam+ theme(legend.position="none"),
    D_Islam_div+ theme(legend.position="none"),
    nrow=3,ncol=2)
mm<-ggplot()+geom_point(aes(x=seq(1,3),y=seq(1,3),color=cbPaletteee))+scale_color_manual(values=cbPaletteee ,labels=c('Real','Binomial','Scaled and rounded raw data'))+labs(colour='Dataset')+theme(legend.position='bottom')

plot_grid( pp, get_legend(mm),axis='b',align='v',ncol=1,rel_heights=c(20,1))
