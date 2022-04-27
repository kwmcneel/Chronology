##the file SST012 would look like this...
#Site 	Species   ID 	Yearoto annulus  width JanSST0  FebSST0  MarSST0 ...
#1 01.BIGS      BR BRBIGS45    2006       2 154.2262 12.97074 12.58376 12.23892 ...
#2 01.BIGS      BR BRBIGS45    2007       3 110.0932 12.18314 12.64134 12.09585 ...
#3 01.BIGS      BR BRBIGS45    2008       4 118.9007 12.05271 11.78405 11.45861 ... 
#4 01.BIGS      BR BRBIGS45    2009       5 147.5249 12.50552 12.53182 11.87801 ...
#5 01.BIGS      BR BRBIGS46    2007       2 191.6801 12.18314 12.64134 12.09585 ...
#6 01.BIGS      BR BRBIGS46    2008       3 138.8029 12.05271 11.78405 11.45861 ...
#



#new variables
SST012$logannulus<-log(SST012$annulus)
SST012$logwidth<-log(SST012$width)

#scale variables so that coefficients are comparable across models
SST012$FallSST0<-scale((OctSST0+NovSST0+DecSST0)/3)
SST012$WinterSST0<-scale((JanSST0+FebSST0+MarSST0)/3)
SST012$SpringSST0<-scale((AprSST0+MaySST0+JunSST0)/3)
SST012$SummerSST0<-scale((JulSST0+AugSST0+SepSST0)/3)
SST012$FallSST1<-scale((OctSST1+NovSST1+DecSST1)/3)
SST012$WinterSST1<-scale((JanSST1+FebSST1+MarSST1)/3)
SST012$SpringSST1<-scale((AprSST1+MaySST1+JunSST1)/3)
SST012$SummerSST1<-scale((JulSST1+AugSST1+SepSST1)/3)
SST012$FallSST2<-scale((OctSST2+NovSST2+DecSST2)/3)
SST012$WinterSST2<-scale((JanSST2+FebSST2+MarSST2)/3)
SST012$SpringSST2<-scale((AprSST2+MaySST2+JunSST2)/3)
SST012$SummerSST2<-scale((JulSST2+AugSST2+SepSST2)/3)



#split data by species
BR<-SST012[SST012$Species=="BR",]
KG<-SST012[SST012$Species=="KG",]


require(lme4)
require(nlme)
require(mgcv)
require(plyr)
require(effects)

######################################
#BLACK ROCKFISH both currents
######################################

#clean up data to remove na values and only include individuals >4 years old
BR2<-na.omit(BR)
summary<-ddply(BR2, ~ID, summarise,age=max(annulus) )
BR3<-merge(BR2,summary,by="ID")
BR4<-BR3[BR3$age>4,]


#base model without any environmental effects
base<-lme(logwidth ~logannulus +Current ,
random=~logannulus|Site/ID, correlation=corAR1(form=~annulus|Site/ID), BR4, method="ML")
BR4$base.resid<-residuals(base)

#models with environmental effects
Fall0<-lme(logwidth ~logannulus +Current/FallSST0-1,
random=~logannulus|Site/ID, correlation=corAR1(form=~annulus|Site/ID), BR4, method="ML") 
Winter0<-lme(logwidth ~logannulus +Current/WinterSST0-1,
random=~logannulus|Site/ID, correlation=corAR1(form=~annulus|Site/ID), BR4, method="ML") 
Spring0<-lme(logwidth ~logannulus +Current/SpringSST0-1,
random=~logannulus|Site/ID, correlation=corAR1(form=~annulus|Site/ID), BR4, method="ML") 
Summer0<-lme(logwidth ~logannulus +Current/SummerSST0-1,
random=~logannulus|Site/ID, correlation=corAR1(form=~annulus|Site/ID), BR4, method="ML") 

Fall1<-lme(logwidth ~logannulus +Current/FallSST1-1,
random=~logannulus|Site/ID, correlation=corAR1(form=~annulus|Site/ID), BR4, method="ML") 
Winter1<-lme(logwidth ~logannulus +Current/WinterSST1-1,
random=~logannulus|Site/ID, correlation=corAR1(form=~annulus|Site/ID), BR4, method="ML") 
Spring1<-lme(logwidth ~logannulus +Current/SpringSST1-1,
random=~logannulus|Site/ID, correlation=corAR1(form=~annulus|Site/ID), BR4, method="ML") 
Summer1<-lme(logwidth ~logannulus +Current/SummerSST1-1,
random=~logannulus|Site/ID, correlation=corAR1(form=~annulus|Site/ID), BR4, method="ML") 
Fall2<-lme(logwidth ~logannulus +Current/FallSST2-1,
random=~logannulus|Site/ID, correlation=corAR1(form=~annulus|Site/ID), BR4, method="ML") 
Winter2<-lme(logwidth ~logannulus +Current/WinterSST2-1,
random=~logannulus|Site/ID, correlation=corAR1(form=~annulus|Site/ID), BR4, method="ML") 
Spring2<-lme(logwidth ~logannulus +Current/SpringSST2-1,
random=~logannulus|Site/ID, correlation=corAR1(form=~annulus|Site/ID), BR4, method="ML") 
Summer2<-lme(logwidth ~logannulus +Current/SummerSST2-1,
random=~logannulus|Site/ID, correlation=corAR1(form=~annulus|Site/ID), BR4, method="ML") 


#make a dataframe of the AIC values 
BRaic<-as.data.frame(AIC(base,Fall2,Winter2,Spring2,Summer2,
Fall1,Winter1,Spring1,Summer1,Fall0,Winter0,Spring0,Summer0))
BRaic$dAIC<-AIC(base)-BRaic$AIC #dAIC shows improvement from the base model (aic improvment of 2 or more)
BRaic$DAIC<-BRaic$AIC-AIC(Summer0) #Summer0 is the best model in this example and used to calculate DAIC
BRaic<-BRaic[2:13,] #remove the line for the base model, which won't be plotted
BRaic
order<-c(1:12) #needed to plot the results in the right temporal order
species<-rep("BR",12) #needed if I have plots with more than one species
BRaic$names<-rownames(BRaic) 
BRaic<-cbind(BRaic,order,species) #add the order and species columns to the existing dataframe


#extract coefficients from models and add to dataframe
L<-list(Fall2,Winter2,Spring2,Summer2,
Fall1,Winter1,Spring1,Summer1,Fall0,Winter0,Spring0,Summer0)
BRresults<-lapply(L,function(x) summary(x)$tTable[4:5,c(1,5)] )  #all the coefficients and pvalues
BRresults2<-as.data.frame(do.call(rbind, BRresults))  #set up as dataframe
BRresults2$order<-rep(1:12,each=2)
BRresults2$group<-rep(c("BRACC","BRCC"),times=12) #identified the results from each current system

#merge the dataframes of the AIC values with the coefficient values
BRaic2<-merge(BRresults2,BRaic, by="order")
names(BRaic2)[2]<-"coef" 
BRaic2


############
###plot coefficients
#####

full<-BRaic2
names(full)[3]<-"pvalue" 
full
#this column will only show the coefficient value if there is AIC improvment from base model
full$coef.mod<-ifelse(full$dAIC<2,full$coef*0,full$coef)
#this column will only show coefficient value if the pvalue of the coefficient is <0.05 and AIC improvement of model >2
full$coef.mod2<-ifelse(full$pvalue<0.05,full$coef.mod,full$coef.mod*0)
#this column will only show coefficient value if the p value of the coefficients is <0.01 and AIC improvement is >2
full$coef.mod3<-ifelse(full$pvalue<0.01,full$coef.mod,full$coef.mod*0)
full<-full[order(full$group) , ]
full

#color and black/white pallets
CL<-c("dodgerblue3","skyblue1","tomato4","tomato1")
BW<-c("black","darkgrey","white","lightgrey")

library(ggplot2)

plot<-ggplot(data=full, aes(x=order,y=coef.mod2,fill=group,group=group)) +
geom_rect(aes(xmin = 4.5, xmax = 8.5,ymin = -Inf, ymax = Inf),fill="gray90", alpha = 0.4) +
geom_bar(stat="identity",position="dodge")+
geom_bar(stat="identity",position="dodge",colour="black",show_guide=FALSE)+
theme_bw()+
scale_y_continuous("Model Coefficient",limits=c(-.12,.15))+
theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill="transparent", colour=NA),
    axis.text.x=element_text(angle=45,hjust=1))+
scale_fill_manual(values=CL,name="",labels=c("Alaska rockfish","California rockfish","Alaska greenling","California greenling"))+
ggtitle("SST")+ 
geom_hline(yintercept=0)+
annotate("text",x=2.5,y=-.12,label="2 yrs prior")+
annotate("text",x=6.5,y=-.12,label="1 yr prior")+
annotate("text",x=10.5,y=-.12,label="growth yr")+
theme(legend.position="bottom", legend.direction="horizontal")+
scale_x_discrete("Season", breaks=1:12,labels=rep(c("Fall","Winter","Spring","Summer"),times=3))+
geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5),linetype="dotted")
plot


#export plot as png
png("SST_quarterEx.png",width=7,height=6,units="in", res=300, bg="transparent")
plot
dev.off()

#export as eps if needed
postscript("SST_quarterEx.eps",width = 7.0, height = 6.0,
           horizontal = FALSE, onefile = FALSE)
plot
dev.off()



#example of some diagnostics
mod<-Summer0  #place to switch model name
print(summary(mod))
require(lattice)
xyplot(resid(mod)~annulus|ID,data=BR4,layout=c(3,4)) #plot of individual residuals
plot(augPred(mod, primary= ~logannulus, level = 2),layout=c(3,4)) #plot of individual fits by individual
xyplot(resid(mod)~Yearoto|Site,data=BR2) #residuals by year and location


#effect plot example
#must rerun the model with an interaction using * (rather than /) for the effect package to work
Summer0.<-lme(logwidth ~logannulus + Current*SummerSST0 ,
random=~1|Site/ID, correlation=corAR1(form=~annulus|Site/ID), BR4, method="ML") 
plot(effect("Current*SummerSST0",Summer0.))
#keep in mind that the 'rug' shown is for all samples, not just those in each current system




