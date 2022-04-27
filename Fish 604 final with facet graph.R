#####################################################################

# FISH 604: Modern Applied Statistics for Fisheries
# Lab 7:  Model fitting
# Kevin McNeel

#####################################################################

require(lme4)
require(nlme)
require(mgcv)
require(plyr)
require(effects)
library(lattice)
library(lmtest)
library(AICcmodavg)
library(ggplot2)
######################## Import and Scale

##Import Ringwidths, Sum data, and Raw data
ShortRF = read.csv(file.choose()) 
Gak1 = read.csv(file.choose()) #for plotting

###################Unused but useful###################
##ShortRF = read.table("clipboard", header=T, sep="\t")
##Environ = read.table("clipboard", header=T, sep="\t")
##PDO= read.table("clipboard", header=T, sep="\t") #for plotting
##########Dataframe cleaning
#Environ<-unique(Environ) ##Remove Duplicates
#Environ<-na.omit(Environ) ##Remove NAs
#write.csv(Gak1, file = file.choose(), na="")##export dataframe as csv
#ShortRF<-edit(ShortRF) ##Directly edit dataframe
##############################################3#######

##Sumarize raw Gak1 for graphing
Gak1Tempsum<-aggregate( Temp.C~Month+Depth.m, Gak1, mean )#Summarize monthly Temp
Gak1Saltsum<-aggregate( Sal.ppt~Month+Depth.m, Gak1, mean )
#PDOsum<-aggregate(NPGO.index~YEAR, PDO, mean )#Summarize monthly Temp

#######Unused but useful: Summarize Ringwidth data and merge width and temp data
#summary(ShortRF)
#head(ShortRF)
#Short1<-ShortRF ##store a reserve 
#ShortRF<-merge(x=ShortRF,y=Environ,by.x=c("YEAR"),by.y=c("YEAR"), all.x=TRUE)
#ShortRF<-merge(x=ShortRF,y=PDOsum,by.x=c("YEAR"),by.y=c("YEAR"), all.x=TRUE)
########
catch=Environ
catch[,2:7]=scale(catch[,2:7])

N_columns = length(names(catch))
names(catch[,2:N_columns])                                                                                                                                                        #Shows starting and ending column names
                                                                                                                                                                                                                                                #,# refers to column of data you want to start normalizing   
x.norm = matrix(scale(catch[,2:N_columns]), ncol=dim(catch[,2:N_columns])[2])     #normalizes
colnames(x.norm) = names(catch[,2:N_columns])                                                                                                              #names columns with original names
catch[,2:N_columns] = x.norm                                                                                                                      #puts in normalized columns in new matrix
x.norm                                                                                                                                                                                                                 #make sure you call new data matrix in models                                                                                                                                                                                                          
                                                                                                                                                                                                                                                #make sure all columns you want to normalize are
                                                                                                                                                                                                                         #numeric or won't run; to make numeric go to Excel
?scale_colour_manual()                                                                                                                                                                                                                                                #file and formate data to number

##Plot monthly temp means by depth
e1<-ggplot(Environ, aes(x=Month,y=Value,colour=factor(Depth),group=Depth))+
	geom_line(size=1)+
	theme(panel.background = element_rect(fill = "NA", color = "gray50"),
		panel.grid.major = element_line(color = "gray90"),
		panel.grid.minor = element_blank(),
		strip.background = element_rect(colour="NA", fill="NA", size=1.5),
		axis.title.y = element_blank())+
	scale_x_continuous(breaks=seq(0, 12, 1))+
	scale_colour_brewer(name  ="Depth",palette = "Spectral")+
	facet_grid(VAR ~ .,scales="free_y",switch = "y")
	

##Plot monthly sal means by depth
e2<-ggplot(Gak1Saltsum, aes(x=Month,y=Sal.ppt,colour=factor(Depth.m),group=Depth.m))+
	geom_line(size=1)+
	ggtitle("(b)")+
	theme(panel.background = element_blank(),panel.grid.major = element_line(color = "gray90"))+
	scale_x_continuous(breaks=seq(0, 12, 1))+
	scale_colour_grey(name  ="Depth m")
##Plot individual growth
ggplot(Environ, aes(x=YEAR,y=Sablefish))+
gad<-ggplot()+
	geom_line(data=Environ,aes(x=YEAR,y=(P.cod)),size=1,linetype=1)+
	geom_line(data=Environ,aes(x=YEAR,y=(Pollock)),size=1,linetype=2)+
	xlab("Year")+
	ggtitle("(a)")+
	theme(panel.background = element_blank())+
	ylab("Round Weight (lbs)")+
	theme(axis.title.y=element_text(vjust=1))+
	scale_x_continuous(breaks=c(2006,2008,2010,2012,2014))
sab<-ggplot()+
	geom_line(data=Environ,aes(x=YEAR,y=(Sablefish)),size=1,linetype=1)+
	xlab("Year")+
	ggtitle("(b)")+
	theme(panel.background = element_blank(),
		axis.title.y=element_blank(),
		legend.title=element_blank(),
		legend.position="none")+
	scale_x_continuous(breaks=c(2006,2008,2010,2012,2014))
roc<-ggplot()+
	geom_line(data=Environ,aes(x=YEAR,y=Shortraker.rf),size=1,linetype=1)+
	geom_line(data=Environ,aes(x=YEAR,y=Rougheye.rf),size=1,linetype=2)+
	xlab("Year")+
	ggtitle("(c)")+
	theme(panel.background = element_blank(),axis.title.y=element_blank(),legend.title=element_blank())+
	scale_x_continuous(breaks=c(2006,2008,2010,2012,2014))
multiplot( e1,e2, cols=2)

#Scale
Short<-ShortRF #backup
ShortRF[,-c(1,2,3,4)]<-scale(ShortRF[,-c(1,2,3,4)])#scales environmental data
summary(Short)
summary(ShortRF)#compare dataframes

#new variables to check width~age transformation
ShortRF$logAGE<-log(ShortRF$AGE)
ShortRF$logINC<-log(ShortRF$INC)

#####base models without any environmental effects
summary(ShortRF$INC)

base<-lme(INC~AGE,random=~AGE|ID, ShortRF, method="ML") #Age not significant
summary(base)
ar(resid(base))
lnbase<-lme(logINC~AGE,random=~AGE|ID, ShortRF, method="ML")
summary(lnbase) #relationship not significant
plot(lnbase,main="ln(Y)")

lnlnbase<-lme(logINC~logAGE,random=~logAGE|ID, ShortRF, method="ML")
summary(lnlnbase) #not significant
plot(lnlnbase,main="ln:ln") 
plot(Shortraker.rf~YEAR,data=Environ)
lnbasecor<-lme(logINC~AGE,
random=~AGE|ID, correlation=corAR1(form=~AGE), ShortRF, method="ML")
#Error running the model

ShortRF$base.resid<-residuals(base)#residuals for plotting

#Remove NAs
ShortRF2<-na.omit(ShortRF)

#####models with environmental effects

##Check transformation and autocorrelation
#corARMA
#correlation=corAR1(form=~AGE|ID)

AVE0<-lme(INC ~AGE + AVE0,random=~AGE|ID, 
ShortRF2, na.action = na.exclude, method="ML") 
summary(AVE0)
plot(ranef(AVE0))
plot(AVE0,main="Linear")#appears normal
r.auto=resid(AVE0)
pacf(r.auto, lag.max=10)
(arma<-round(ar(resid(AVE0))$ar,2))

AVE0<-lme(INC ~AGE + AVE0,random=~AGE|ID, correlation=corARMA(p=12,form=~AGE|ID,value=c(0.09,0.16,-0.021,-0.11,0.14,-0.16,-0.09,-0.07,0.01,-0.10,-0.09,-0.11),fixed=TRUE),
ShortRF2, na.action = na.exclude, method="ML") 
summary(MAX250)
plot(ranef(AVE0arma))
plot(AVE0arma,main="Linear")#appears normal
r.auto=resid(AVE0arma)
pacf(r.auto, lag.max=10)
ar(resid(AVE01))


logAVE0<-lme(logINC ~AGE + AVE0,random=~AGE|ID, correlation=corAR1(form=~AGE|ID),
ShortRF2, na.action = na.exclude, method="ML") 
summary(logAVE0)
plot(logAVE0,main="ln(y)")#appears not normal
?corAR1

AVE01<-lme(INC ~AGE + AVE0,random=~AGE|ID, ShortRF2, method="ML")
n=length(Short
r=resid(AVE01)
ar(r)
plot(r~ShortRF2$YEAR)
abline(h=0)
dwtest(AVE01)
acf(r, lag.max=10)
pacf(r, lag.max=10)
lrAVE0<-lrtest(AVE01,AVE0)	###Including autocorrelation was significant

AVElm<-lme(INC ~AGE + AVE0,random=~1|ID, correlation=corAR1(form=~1|ID),
ShortRF2, na.action = na.exclude, method="ML") 
summary(AVElm)
lrtest(AVE0,AVE0arma)##random effect on age significant
 
##basecor from above with reduced dataframe with no NAs for AIC compare
base<-lme(INC~AGE,
random=~AGE|ID, correlation=corAR1(form=~AGE|ID), ShortRF2, method="ML")
baseresid<-resid(base)

AVE250<-lme(INC ~AGE+AVE250,
random=~AGE|ID, correlation=corAR1(form=~AGE|ID), ShortRF2, method="ML")

AVETOTAL<-lme(INC ~AGE+AVETOTAL,
random=~AGE|ID, correlation=corAR1(form=~AGE|ID), ShortRF2, method="ML") 

MAX0<-lme(INC ~AGE+MAX0,
random=~AGE|ID, correlation=corAR1(form=~AGE|ID), ShortRF2, method="ML") 

MAX250<-lme(INC ~AGE+MAX250,
random=~AGE|ID, correlation=corAR1(form=~AGE|ID), ShortRF2, method="ML") 

MAXTOTAL<-lme(INC ~AGE+MAXTOTAL,
random=~AGE|ID, correlation=corAR1(form=~AGE|ID), ShortRF2, method="ML")

MIN0<-lme(INC ~AGE+MIN0,
random=~AGE|ID, correlation=corAR1(form=~AGE|ID), ShortRF2, method="ML")

MIN250<-lme(INC ~AGE+MIN250,
random=~AGE|ID, correlation=corAR1(form=~AGE|ID), ShortRF2, method="ML") 

MINTOTAL<-lme(INC ~AGE+MINTOTAL,
random=~AGE|ID, correlation=corAR1(form=~AGE|ID), ShortRF2, method="ML") 


AVE250<-lme(INC ~AGE+AVE250, random=~AGE|ID, correlation=corARMA(p=10,form=~AGE|ID,value=c(0.16,0.17,-0.06,-0.11,0.16,-0.14,-0.08,-0.04,-0.02,-0.15),fixed=TRUE), ShortRF2, method="ML")

AVETOTAL<-lme(INC ~AGE+AVETOTAL, random=~AGE|ID,  correlation=corARMA(p=7,form=~AGE|ID,value=c(0.12,0.22,-0.02,-0.08,0.15,-0.16,-0.11),fixed=TRUE),ShortRF2, method="ML") 
AVETOTAL1<-lme(INC ~AGE+AVETOTAL, random=~1|ID,  correlation=corARMA(p=7,form=~AGE|ID,value=c(0.12,0.22,-0.02,-0.08,0.15,-0.16,-0.11),fixed=TRUE),ShortRF2, method="ML") 
lrtest(AVETOTAL1,AVETOTAL)

MIN0<-lme(INC ~AGE+MIN0,random=~AGE|ID, correlation=corARMA(p=10,form=~AGE|ID,value=c(0.15,0.12,-0.02,-0.06,0.13,-0.15,-0.06,-0.04,-0.03,-0.18),fixed=TRUE), ShortRF2, method="ML") 

MIN250<-lme(INC ~AGE+MIN250,random=~AGE|ID, correlation=corARMA(p=10,form=~AGE|ID,value=c(0.16,0.17,-0.05,-0.11,0.16,-0.15,-0.08,-0.04,-0.02,-0.17),fixed=TRUE),ShortRF2, method="ML") 

MINTOTAL<-lme(INC ~AGE+MINTOTAL,random=~AGE|ID, correlation=corARMA(p=10,form=~AGE|ID,value=c(0.16,0.12,-0.02,-0.05,0.14,-0.15,-0.06,-0.04,-0.03,-0.17),fixed=TRUE),  ShortRF2, method="ML")

MAX0<-lme(INC ~AGE+MAX0,random=~AGE|ID,correlation=corARMA(p=7,form=~AGE|ID,value=c(0.14,0.18,0,-0.08,0.16,-0.15,-0.12),fixed=TRUE),  ShortRF2, method="ML")

MAX250<-lme(INC ~AGE+MAX250,random=~AGE|ID,correlation=corARMA(p=6,form=~AGE|ID,value=c(0.19,0.17,-0.02,-0.09,0.11,-0.14),fixed=TRUE),  ShortRF2, method="ML") 
MAX2501<-lme(INC ~AGE+MAX250,random=~1|ID,correlation=corARMA(p=6,form=~AGE|ID,value=c(0.19,0.17,-0.02,-0.09,0.11,-0.14),fixed=TRUE),  ShortRF2, method="ML") 
lrtest(MAX250,MAX2501)
MAXTOTAL<-lme(INC ~AGE+MAXTOTAL,random=~AGE|ID,correlation=corARMA(p=7,form=~AGE|ID,value=c(0.14,0.18,0, -0.08,0.16,-0.15,-0.12),fixed=TRUE),  ShortRF2, method="ML") 

TOTAL<-lme(INC ~AGE+AVE0+AVE250,
random=~AGE|ID,  ShortRF2, method="ML") 
summary(TOTAL)
library(MASS),
stepAIC(TOTAL,direction="both",steps = 1000)

round(ar(resid(AVE0))$ar,2);round(ar(resid(AVE250))$ar,2);round(ar(resid(AVETOTAL))$ar,2)
round(ar(resid(MIN0))$ar,2);round(ar(resid(MIN250))$ar,2);round(ar(resid(MINTOTAL))$ar,2)
round(ar(resid(MAX0))$ar,2);round(ar(resid(MAX250))$ar,2);round(ar(resid(MAXTOTAL))$ar,2)
pacf(resid(MIN0))
#make a dataframe of the AIC values 
ShortRFaic<-as.data.frame(AIC(base,AVE0,AVE250,AVETOTAL,MIN0,MIN250,MINTOTAL,
MAX0,MAX250,MAXTOTAL))
ShortRFaic$dAIC<-AIC(base)-ShortRFaic$AIC #dAIC shows improvement from the base model (aic improvment of 2 or more)
ShortRFaic$DAIC<-ShortRFaic$AIC-min(ShortRFaic$AIC) #AVETOTAL is the best model in this example and used to calculate DAIC
ShortRFaic<-ShortRFaic[2:10,] #remove the line for the base model, which won't be plotted
ShortRFaic
order<-c(1:9) #needed to plot the results in the right temporal order
ShortRFaic$names<-rownames(ShortRFaic) 
ShortRFaic<-cbind(ShortRFaic,order) #add the order

#extract coefficients from models and add to dataframe
L<-list(AVE0,AVE250,AVETOTAL,MIN0,
MIN250,MINTOTAL,MAX0,MAX250,MAXTOTAL)
ShortRFresults2<-lapply(L,function(x) summary(x)$tTable[2,c(1,4,5)] )  #all the coefficients and pvalues
ShortRFresults2<-as.data.frame(do.call(rbind, ShortRFresults2))  #set up as dataframe
ShortRFar2<-lapply(L,function(x) ar(resid(x))$order )  #all the coefficients and pvalues
ShortRFar2<-as.data.frame(do.call(rbind, ShortRFar2))  #set up as dataframe

ShortRFresults2$order<-rep(1:9,)
ShortRFaicc<-lapply(L,function(x) AICc(x)) 
ShortRFaicc<-as.data.frame(do.call(rbind, ShortRFaicc))
colnames(ShortRFaicc)<-c("AICc") 
ShortRFaicc$DAICc<-ShortRFaicc$AICc-min(ShortRFaicc$AICc)
ShortRFaicc$order<-rep(1:9,)

#merge the dataframes of the AIC values with the coefficient values
ShortRFaic<-merge(ShortRFaic,ShortRFaicc, by="order")
ShortRFaic<-merge(ShortRFaic,ShortRFresults, by="order")

ShortRFaic
ShortRFaic$AR<-lapply(L,function(x) ar(resid(x))$order)  #all the coefficients and pvalues
summary(AVE0)
############
###plot coefficients
#####

full<-ShortRFaic
colnames(full)[c(9,10)]<-c("coef","pvalue") 
full
#this column will only show the coefficient value if there is AIC improvment from base model
full$coef.AICc<-ifelse(full$DAICc>2,full$coef*0,full$coef)
#this column will only show coefficient value if the pvalue of the coefficient is <0.05 and AIC improvement of model >2
full$coef.05<-ifelse(full$pvalue<0.05,full$coef.AICc,full$coef.AICc*0)
#this column will only show coefficient value if the p value of the coefficients is <0.01 and AIC improvement is >2
full$coef.01<-ifelse(full$pvalue<0.01,full$coef.AICc,full$coef.AICc*0)
full

###diagnostics
mod<-AVE0 #place to switch model name
print(summary(mod))
xyplot(resid(mod)~AGE|ID,data=ShortRF2) #plot of individual residuals
plot(augPred(mod, primary= ~AGE)) #plot of individual fits by individual
xyplot(resid(mod)~YEAR,data=ShortRF2,type="b") #residuals by year and location
chronology<-aggregate( scale(resid(mod))~YEAR, ShortRF2, mean )
plot(chronology,type="b")
r<-resid(mod)


#effect plot##Average Temperature Response across other effects
plot(effect("AVE0",AVE0),main=NA)
plot(effect("AVE250",AVE250),main=NA)
plot(effect("AVETOTAL",AVETOTAL),main=NA)
plot(effect("MIN0",MIN0),main=NA)
plot(effect("MIN250",MIN250),main=NA)
plot(effect("MINTOTAL",MINTOTAL),main=NA)
plot(effect("MAX0",MAX0),main=NA)
plot(effect("MAX250",MAX250),main=NA)
plot(effect("MAXTOTAL",MAXTOTAL),main=NA)

#########TO DO#########

AVE02 <- lme(INC ~AGE + AVETOTAL,
random=~AGE|ID, ShortRF2, correl = corARMA( p=4, form =~AGE|ID), method="ML")
names(ShortRF2)
#Test transformation on larger datasets/more measurements

library(corrgram)
corrgram(ShortRF2[,c(3,7,12,9,13,16)], lower.panel=panel.conf, upper.panel=panel.ellipse,
         diag.panel=panel.density)
corrgram(ShortRF2, lower.panel=panel.conf, upper.panel=panel.ellipse,
         diag.panel=panel.density)

##Sumarize yearly residual
Baseresid1<-aggregate(baseresid~YEAR, ShortRF2, mean )

##Plot individual growthcolor="#56B4E9",color="Residual"),color="#999999""#E69F00"
p1<-ggplot()+
	geom_line(data=ShortRF2,aes(x=YEAR,y=MIN0,color="MIN0"),color="grey50",size=1,linetype=6)+
	geom_line(data = Baseresid1, aes(x = YEAR, y = scale(baseresid),color="AveAge|ID Resid"),color="black",size=1.2,linetype=1)+
	xlab("Year")+
	ylab("Scaled Response")+
	theme(legend.title=element_blank())+
      scale_x_continuous(breaks=seq(1965,2015,5))+
	ggtitle("(b)")+
	theme(panel.background = element_blank(),panel.grid.major = element_line(color = "gray90"))

p2<-ggplot()+
	geom_line(data=ShortRF2,aes(x=YEAR,y=AVE0,color="AVE0"),color="grey50",size=1,linetype=6)+
	geom_line(data = Baseresid1, aes(x = YEAR, y = scale(baseresid),color="AveAge|ID Resid"),color="black",size=1.2,linetype=1)+
	xlab("Year")+
	ylab("Scaled Response")+
	theme(legend.title=element_blank())+
	ggtitle("(a)")+
      scale_x_continuous(breaks=seq(1965,2015,5))+
	theme(panel.background = element_blank(),panel.grid.major = element_line(color = "gray90"))


p3<-ggplot()+
	geom_point(data=ShortRF2,aes(x=AGE,y=INC,color=YEAR),size=2)+
  	theme(axis.title.y=element_text(vjust=1))+
	geom_smooth(data=ShortRF2,aes(x=AGE,y=INC),alpha=.2, size=1,col="black",se=F)+
	xlab("Age")+
	theme(panel.background = element_blank(),axis.title.y = element_blank())+
	ggtitle("(b)")+
	scale_colour_gradientn(colours=rainbow(5))

p4<-ggplot()+
	geom_boxplot(data=ShortRF2,aes(x=ID,y=INC),fill="lightblue",size=1)+
	scale_x_discrete(breaks=NULL)+
	xlab("Individual")+
	ylab("Increment Width (mm)")+
	theme(axis.title.y=element_text(vjust=1))+
	theme(legend.title=element_blank())+
	ggtitle("(a)")+
	theme(panel.background = element_blank())


#export plot as png
png("ShortRFdiagnostic.png",width=8,height=4,units="in", res=300, bg="transparent")
par(mfrow=c(1,2));pacf(r,main="(a)");plot(r~fitted(mod),xlab=expression(hat(y)),main="(b)");abline(h=0,lty=2)
dev.off()

png("ShortRFggplot.png",width=8,height=4,units="in", res=300, bg="transparent")
multiplot(p4, p3, cols=2)
dev.off()


png("ShortRFGAK1.png",width=8,height=4,units="in", res=300, bg="transparent")
multiplot(e1, e2, cols=1)
dev.off()

png("ShortRFtempcurve.png",width=8,height=4,units="in", res=300, bg="transparent")
p1
dev.off()

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}