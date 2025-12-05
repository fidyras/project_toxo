
## Location of all of the data (change to your path)
fname <- "~/Princeton Dropbox/C. Jessica Metcalf/Toxoplasmosis/data/"


#### 1. BRING IN THE DATA ##################################################################################################################################################

## A. Data on seroprevalence  ###############
df <- read.csv(paste(fname,"Updated-data-added-Ingrid-Jun28-2025.csv", sep=""))

## Organize names
year <- as.numeric(substring(df$Year.of.data.collection..Final.Year.of.Data.Collection.,1,4))
country <- df$Location..country.
city <- df$Location..town
u.ctry <- unique(country)

## Get the core variables
lower.age <- as.numeric(df$Lower.bound.of.age)
upper.age <- as.numeric(df$Upper.bound.of.age) 
iggpos <- as.numeric(df$Number.of.people.seropositive)
tot <- as.numeric(df$Number.of.people.tested)

#Assume that upper age for missing and putatively reproductive population is 50 
upper.age[is.na(upper.age) & df$Women.of.Reproductive.age==1] <- 50

#Round for the one digitized from Cameroon
iggpos <- round(iggpos)

#Get mid ages
mid.ages <- 0.5*(lower.age+upper.age)






## B. HDI data by year / country ################################################################################################################################################
HDI2 <- read.csv(paste(fname, "human-development-index/human-development-index.csv",sep=""))
HDI2$Entity[HDI2$Entity =="Cote d'Ivoire"] <- "Ivory Coast"
HDI2$Entity[HDI2$Entity =="Sao Tome and Principe"] <- "Democratic Republic of São Tomé and Príncipe"
HDI2$Year <- as.numeric(HDI2$Year)

## Align the full data-frame 
hdi.match <- rep(NA,length(year))
for (j in 1:length(year)) {
   if (is.na(country[j])) next()
   if (is.na(year[j])) next()
  chs <-  HDI2$Year==year[j] & HDI2$Entity==country[j]
  if (sum(chs)>0) hdi.match[j] <- HDI2$Human.Development.Index[chs]
  if (sum(chs)==0) hdi.match[j] <- HDI2$Human.Development.Index[HDI2$Entity==country[j] & HDI2$Year==2022] #if no data available, take the most recent

}

#store the most recent for future projection
recent.hdi <- rep(NA,length(u.ctry)) 
for (j in 1:length(u.ctry)) recent.hdi[j] <-HDI2$Human.Development.Index[HDI2$Year==2022 & HDI2$Entity==country[j]]








## C. Temperature data (Wenchang) ##############################################################################################################################################
temps <- read.csv(paste(fname, "/climate-data-from-Wenchang/era5.skt.daily.Ingrid.1979-2024.degC.csv", sep=""))
country.temps <- city.temps <- rep(NA,ncol(temps)-1)
for (j in 2:ncol(temps)) { country.temps[j-1] <- substring(colnames(temps)[j],1,regexpr("_",colnames(temps)[j])[[1]]-1); 
			    city.temps[j-1] <-  substring(colnames(temps)[j],regexpr("_",colnames(temps)[j])[[1]]+1,nchar(colnames(temps)[j]))}
year.temperature <- as.numeric(substring(temps$time,1,4))
city.temps[city.temps=="SãoTomé"] <- "São Tomé"
city.temps[city.temps=="PenkaMichel"] <- "Penka-Michel"
city.temps[city.temps=="Burch.i"] <- "Burch'i"
city.temps[city.temps=="Jos.North"] <- "Jos–North"

country.temps[country.temps=="BurkinaFaso"] <- "Burkina Faso"
country.temps[country.temps=="DemocraticRepublicofSãoToméandPríncipe"] <- "Democratic Republic of São Tomé and Príncipe"
country.temps[country.temps=="IvoryCoast"] <- "Ivory Coast"


## Find temperature for each year/country/city combo
min.temperature.match <- mean.temperature.match <- median.temperature.match <- max.temperature.match <- rep(NA,length(year)) 
for (j in 1:length(year)) {
	  if (is.na(country[j])) next()
	  row.chs <- which(year.temperature==year[j])
 	  col.chs <- which(city.temps== gsub(" ","",city[j])) 
	  if (length(col.chs)==0) col.chs <-  which(city.temps== gsub("-","",city[j])) 
	  if (length(col.chs)==0) col.chs <- which(country.temps ==country[j]) 

	  #indexing - note, need to add one to the columns to account for the first colum of time. 
	  tmp.here <- as.numeric(c(as.matrix(temps[row.chs,col.chs+1])))	
	
	  mean.temperature.match[j] <- mean(tmp.here)	
	  median.temperature.match[j] <- median(tmp.here)	
	  max.temperature.match[j] <- max(tmp.here)	
	  min.temperature.match[j] <- min(tmp.here)	

}



mean.temperature.match[!is.finite(mean.temperature.match)] <- NA
max.temperature.match[!is.finite(max.temperature.match)] <- NA
median.temperature.match[!is.finite(median.temperature.match)] <- NA
min.temperature.match[!is.finite(min.temperature.match)] <- NA


early.mean.temp <- early.max.temp <- early.min.temp <- rep(NA,length(u.ctry))
recent.mean.temp <- recent.max.temp <- recent.min.temp <- rep(NA,length(u.ctry))
mean.temp.u.ctry <- max.temp.u.ctry <- min.temp.u.ctry <- rep(NA,length(u.ctry))
for (j in 1:length(u.ctry)) { 
	col.chs <- which(country.temps ==u.ctry[j],arr.ind=TRUE)
        row.chs <- which(year.temperature==2024,arr.ind=TRUE)
        row.chs2 <- which(year.temperature==1979,arr.ind=TRUE)

	mean.temp.u.ctry[j] <- mean(as.numeric(c(as.matrix(temps[,col.chs+1]))))	
	max.temp.u.ctry[j] <- max(as.numeric(c(as.matrix(temps[,col.chs+1]))))	
	min.temp.u.ctry[j] <- min(as.numeric(c(as.matrix(temps[,col.chs+1]))))	


	recent.mean.temp [j] <- mean(as.numeric(c(as.matrix(temps[row.chs,col.chs+1]))))	
	recent.max.temp[j] <- max(as.numeric(c(as.matrix(temps[row.chs,col.chs+1]))))	
	recent.min.temp[j] <- min(as.numeric(c(as.matrix(temps[row.chs,col.chs+1]))))	

	early.mean.temp [j] <- mean(as.numeric(c(as.matrix(temps[row.chs2,col.chs+1]))))	
	early.max.temp[j] <- max(as.numeric(c(as.matrix(temps[row.chs2,col.chs+1]))))	
	early.min.temp[j] <- min(as.numeric(c(as.matrix(temps[row.chs2,col.chs+1]))))	


}




## Get most recent scores (for current burden) and order by unique country 
recent.data <- data.frame(country=u.ctry, recent.hdi= recent.hdi,
	recent.mean.temp=recent.mean.temp,recent.max.temp=recent.max.temp,recent.min.temp=recent.min.temp,
	early.mean.temp=early.mean.temp,early.max.temp=early.max.temp,early.min.temp=early.min.temp,  
		mean.temp= mean.temp.u.ctry, max.temp= max.temp.u.ctry, min.temp= min.temp.u.ctry)
recent.data <- recent.data[order(recent.data[,1]),]






## D. Fertility over age ######################
### Data is from UN population projectoins: https://population.un.org/wpp/downloads?folder=Standard%20Projections&group=Fertility
### Births by single age of mother annually in 1000s - now and in future (so presumably correcting for changing populations)


## past and current estimates
estimate.births <- read.csv(paste(fname,"/UNdataFertPop/Estimates.births.by.age.mother.csv",sep=""))
year.estimate <- estimate.births$Year
births.per.age <- estimate.births[,12:ncol(estimate.births)]

## Fiddle names
estimate.births$Region..subregion..country.or.area..[estimate.births $Region..subregion..country.or.area..=="United Republic of Tanzania"] <- "Tanzania"
estimate.births$Region..subregion..country.or.area..[estimate.births $Region..subregion..country.or.area..=="Côte d'Ivoire"] <- "Ivory Coast"
estimate.births$Region..subregion..country.or.area..[estimate.births $Region..subregion..country.or.area..=="Sao Tome and Principe"] <- "Democratic Republic of São Tomé and Príncipe"




## future medium variant 
median.future.births <- read.csv(paste(fname,"/UNdataFertPop/MediumVariantPorjection.births.by.age.mother.csv",sep=""))
year.median <- median.future.births$Year
births.per.age.future <- median.future.births[,12:ncol(median.future.births)]


## Fiddle names
median.future.births$Region..subregion..country.or.area..[median.future.births $Region..subregion..country.or.area..=="United Republic of Tanzania"] <- "Tanzania"
median.future.births$Region..subregion..country.or.area..[median.future.births $Region..subregion..country.or.area..=="Côte d'Ivoire"] <- "Ivory Coast"
median.future.births$Region..subregion..country.or.area..[median.future.births $Region..subregion..country.or.area..=="Sao Tome and Principe"] <- "Democratic Republic of São Tomé and Príncipe"









#### 2. CREATE A DATA FRAME WWITH PUBLICATION, POS, TOT, AGE, TEMPERATURE, HDI ###############################################################################################################

## TODO check if want mid or upper.age
new.data <- data.frame(lower.age=lower.age,upper.age=upper.age,Pos=iggpos,Neg=tot-iggpos,tot=tot,lAge=log(mid.ages),#lAge=log(upper.age),
			year=year, country= country,city=city,hdi= hdi.match,
			mean.temperature=mean.temperature.match, max.temperature= max.temperature.match,
			median.temperature= median.temperature.match, min.temperature= min.temperature.match)


## some convenience vectors
bad <- is.na(new.data$Pos) | is.na(new.data$Neg) | is.na(new.data$lAge) | !is.finite(new.data$lAge)
test.hdi <- seq(0.3,1,length=100)
test.temp <- seq(14,35,length=100)
test.min.temp <- seq(1,26,length=100)
test.max.temp <- seq(20,42,length=100)
u.ctry <- unique(new.data$country[!bad])[order(unique(new.data$country[!bad]))]

## elimatinate the countries that don't make it from recent.data
recent.data <- recent.data[is.element(recent.data[,1],u.ctry),]

## Index the countries that we are focussed on in the fertility data (with the updated u.ctry)
index.u.ctry.in.estimate <- index.u.ctry.in.future <- rep(NA,nrow(estimate.births))
for (j in 1:length(u.ctry)) { 
	index.u.ctry.in.estimate[u.ctry[j]==estimate.births$Region..subregion..country.or.area..] <- j 
	index.u.ctry.in.future[u.ctry[j]==median.future.births$Region..subregion..country.or.area..] <- j 
}
## sanity check
#estimate.births[index.u.ctry.in.estimate==22 & !is.na(index.u.ctry.in.estimate),]








## FiGURES OF THE DATA  #################################

## A. SEROPREV
par(mfrow=c(1,3),bty="n", mar=c(10,5,4,4))
cols <- colorRampPalette(c("blue", "red"))(24)
plot(exp(new.data$lAge[!bad]),new.data$Pos[!bad]/(new.data$tot[!bad]),  pch=19, xlab="Age class midpoint", ylab="IgG Positive (proportion)", col=cols[round(new.data$max.temperature[!bad])-18], cex=new.data$hdi[!bad]*2)
legend("topleft",legend=c(19,30,42), col=cols[c(1,12,24)],pch=15, title="Max T (C)", bty="n")
title("A) Seroprevalence data")

## B. TEMPERATURE
o.temp <- order(recent.data$mean.temp)
plot(c(1:length(u.ctry)), recent.data$mean.temp[o.temp], ylim=range(c(recent.data$mean.temp, recent.data$min.temp, recent.data$max.temp),na.rm=T), xlab="", ylab="Temperature (C)", pch=1, axes=F, type="n")
points(c(1:length(u.ctry))-0.25, recent.data$early.mean.temp[o.temp], pch=19,col="grey")
points(c(1:length(u.ctry))+0.25, recent.data$recent.mean.temp[o.temp], pch=19,col="black")

for (j in  c(1:length(u.ctry))) { 
	#points(c(j,j), c(recent.data$min.temp[o.temp][j], recent.data$max.temp[o.temp][j]),lty=1, type="l", col="black",lwd=0.5)
	points(c(j,j)-0.25, c(recent.data$early.min.temp[o.temp][j], recent.data$early.max.temp[o.temp][j]),lty=1, type="l", col="grey",lwd=0.5)
	points(c(j,j)+0.25, c(recent.data$recent.min.temp[o.temp][j], recent.data$recent.max.temp[o.temp][j]),lty=1, type="l", col="black",lwd=0.5)
	}


axis(1,at=c(1:length(u.ctry)),lab=u.ctry[o.temp],las=2); axis(2)
title("B) Temperature ranges")
#legend("bottomright",legend=c("Average","1979","2023"), pch=c(1,19,19),col=c(1,"grey",1),bty="n")
legend("bottomright",legend=c("1979","2023"), pch=c(19,19),col=c("grey",1),bty="n")


## C. FERTILITY 
j <- 12
births.2023 <- as.numeric(births.per.age[index.u.ctry.in.estimate==j & !is.na(index.u.ctry.in.estimate) & year.estimate==2023,])
births.2050 <- as.numeric(births.per.age.future[index.u.ctry.in.future ==j & !is.na(index.u.ctry.in.future) & year.median ==2050,])
births.2100 <- as.numeric(births.per.age.future[index.u.ctry.in.future ==j & !is.na(index.u.ctry.in.future) & year.median ==2100,])

yupper <- max(c(births.2023, births.2050, births.2100))
plot(15:49, births.2023,xlab="Age", ylab="Births per year age (1000s)", type="l",ylim=c(0,yupper))
points(15:49, births.2050,type="l",lty=2)
points(15:49, births.2100,type="l",lty=3)
title("C) Fertility trajectory example")
#title(u.ctry[j])
legend("topright",legend=c(2023,2050,2100), title=u.ctry[j],lty=c(1,2,3),bty="n")







#### 3. FIT THE MODEL ###########################################################################################################################################################################


require(mgcv)

## need this to add random effect
new.data$country <- as.factor(new.data$country)

## Assume: linear on min and max temperature and HDI; random effect of country 
fit <- gam(cbind(Pos,Neg)~hdi+max.temperature+min.temperature+s(country,bs="re"),offset=lAge,family=binomial(link=cloglog), data=new.data[!bad,])


## Comparison with simpler models
fit.notemp <- gam(cbind(Pos,Neg)~hdi+s(country,bs="re"),offset=lAge,family=binomial(link=cloglog), data=new.data[!bad,])
fit.maxtemp <- gam(cbind(Pos,Neg)~hdi+max.temperature+s(country,bs="re"),offset=lAge,family=binomial(link=cloglog), data=new.data[!bad,])
anova(fit.notemp,fit.maxtemp,fit,test="Chisq") ## does seem to improve


## Check error distribution 
plot(fit) ## reasonable error distribution 


## Plot across the range of variables, others set to mean; and country is Somalia (closest to average random effect); note that lAge is set to 18, but predict does not include the offset
rc.hdi <- predict(fit,newdata=data.frame(hdi=test.hdi, max.temperature =mean(new.data$max.temperature[!bad],na.rm=TRUE), 
		min.temperature =mean(new.data$min.temperature[!bad],na.rm=TRUE),country=as.factor("Somalia"),lAge=log(18)),se.fit = TRUE)
rc.min.temp <- predict(fit,newdata=data.frame(hdi=0.5, min.temperature =test.min.temp,max.temperature=mean(new.data$max.temperature[!bad],na.rm=TRUE),country=as.factor("Somalia"),lAge=log(18)),se.fit = TRUE)
rc.max.temp <- predict(fit,newdata=data.frame(hdi=0.5, max.temperature =test.max.temp, min.temperature =mean(new.data$min.temperature[!bad],na.rm=TRUE),country=as.factor("Somalia"),lAge=log(18)),se.fit = TRUE)
rc.ctry <- predict(fit,newdata=data.frame(hdi=0.5, max.temperature =mean(new.data$max.temperature[!bad],na.rm=TRUE), min.temperature =mean(new.data$min.temperature[!bad],na.rm=TRUE),country=as.factor(u.ctry),lAge=log(18)),se.fit = TRUE)

## convenience to allow toggle between different outputs
conv <- function(x) return(exp(x)) 
#conv <- function(x) return(x)

## Show results
par(mfrow=c(1,4), mar=c(10,5,4,1), bty="l")
plot(test.hdi,conv(rc.hdi$fit), type="l", xlab="Human Development Index", ylab="Force of infection", pch=19, cex.lab=1.5,lwd=2)
points(test.hdi, conv( rc.hdi$fit-1.96*rc.hdi$se.fit), type="l",lty=2)
points(test.hdi, conv( rc.hdi$fit+1.96*rc.hdi$se.fit), type="l",lty=2)
points(new.data$hdi[!bad],rep(min(conv(rc.hdi$fit)),nrow(new.data[!bad,])),pch=3) #tickmarks showing where data is

plot(test.max.temp,conv(rc.max.temp$fit), type="l", xlab="Max temperature (C)", ylab="Force of infection", pch=19, cex.lab=1.5,lwd=2)
points(test.max.temp,conv(rc.max.temp$fit-1.96*rc.max.temp$se.fit), type="l",lty=2)
points(test.max.temp,conv(rc.max.temp$fit+1.96*rc.max.temp$se.fit), type="l",lty=2)
points(new.data$max.temperature[!bad],rep(min(conv(rc.max.temp$fit)),nrow(new.data[!bad,])),pch=3) #tickmarks showing where data is

plot(test.min.temp,conv(rc.min.temp$fit), type="l", xlab="Min temperature (C)", ylab="Force of infection", pch=19, cex.lab=1.5,lwd=2)
points(test.min.temp,conv(rc.min.temp$fit-1.96*rc.min.temp$se.fit), type="l",lty=2)
points(test.min.temp,conv(rc.min.temp$fit+1.96*rc.min.temp$se.fit), type="l",lty=2)
points(new.data$min.temperature[!bad],rep(min(conv(rc.min.temp$fit)),nrow(new.data[!bad,])),pch=3) #tickmarks showing where data is

order.c <- order(rc.ctry$fit)
plot(1:length(u.ctry),conv(rc.ctry$fit[order.c]), type="p", xlab="", ylab="Force of infection", pch=19, cex.lab=1.5,lwd=2, axes=FALSE)
for (j in 1:length(u.ctry)) points(c(j,j),conv(rc.ctry$fit[order.c][j]+c(1.96,-1.96)*rc.ctry$se.fit[order.c][j]), type="l",lty=2)
axis(1,at=1:length(u.ctry),lab=u.ctry[order.c],las=2); axis(2)


## cbind(u.ctry,coef(fit)[5:length(coef(fit))])[order(coef(fit)[5:length(coef(fit))]),] #sanity check this matches what you get with predict for country level






#### 4. ESTIMATE THE BURDEN ###########################################################################################################################################################################
## One issue is that no babies are projected for Sao Tome and Principe....


model.now <-  gam(cbind(Pos,Neg)~hdi+max.temperature+min.temperature+s(country,bs="re"),offset=lAge,family=binomial(link=cloglog), data=new.data[!bad,])


##orders and checks: 
recent.data <- recent.data[order(recent.data[,1]),]
#cbind(recent.data[,1],u.ctry,rownames(ranef(model.now)$country))


## use recent.data which contains most recent HDI and temperature
burden.now <- burden.future.baseline <- burden.future.hotter <- burden.longfuture.baseline <- burden.longfuture.hotter <- rep(NA,length(u.ctry))
burden.now.rel <- burden.future.baseline.rel <- burden.future.hotter.rel <- burden.longfuture.baseline.rel <- burden.longfuture.hotter.rel <- rep(NA,length(u.ctry))

par(mfrow=c(1,4), bty="l")
yupper <- 60 ##upper bound for plots of countries births per age class

## Loop over 
for (j in 1:length(u.ctry)) { 

	## relies on order u.ctry being same as in model.now and in recent.data; and there obly being 4 fixed effects (intercept, slope1,2,3); can check order of coefficients in model.now by doing levels(new.data$country)	
	country.effect <- coef(model.now)[5:length(coef(model.now))][j]

	## get current age profile risk
	risk.now <- exp(coef(model.now)[1]+coef(model.now)[2]*recent.data$recent.hdi[j]+
					coef(model.now)[3]*recent.data$recent.max.temp[j]+coef(model.now)[4]*recent.data$recent.min.temp[j]+ country.effect)
	staying.susceptible <- exp(-risk.now*(15:49))
	first.infection.at.age <- risk.now* staying.susceptible

	## get same with 3 degree increase
	risk.now.hot <- exp(coef(model.now)[1]+coef(model.now)[2]*recent.data$recent.hdi[j]+
					coef(model.now)[3]*(recent.data$recent.max.temp[j]+3)+coef(model.now)[4]*recent.data$recent.min.temp[j]+ country.effect)
	staying.susceptible.hot <- exp(-risk.now.hot*(15:49))
	first.infection.at.age.hot <- risk.now.hot* staying.susceptible.hot

		
	## get births in the future
	births.2023 <- as.numeric(births.per.age[index.u.ctry.in.estimate==j & !is.na(index.u.ctry.in.estimate) & year.estimate==2023,])
	births.2050 <- as.numeric(births.per.age.future[index.u.ctry.in.future ==j & !is.na(index.u.ctry.in.future) & year.median ==2050,])
	births.2100 <- as.numeric(births.per.age.future[index.u.ctry.in.future ==j & !is.na(index.u.ctry.in.future) & year.median ==2100,])

	## get burdens assuming only demography changes
	burden.now[j] <- sum(first.infection.at.age*births.2023)
	burden.future.baseline[j] <- sum(first.infection.at.age*births.2050)
	burden.future.hotter[j] <- sum(first.infection.at.age.hot*births.2050)

	burden.longfuture.baseline[j] <- sum(first.infection.at.age*births.2100)
	burden.longfuture.hotter[j] <- sum(first.infection.at.age.hot*births.2100)


	## get burdens assuming only demography changes
	burden.now.rel[j] <- sum(first.infection.at.age*births.2023)/sum(births.2023)
	burden.future.baseline.rel[j] <- sum(first.infection.at.age*births.2050)/sum(births.2050)
	burden.future.hotter.rel[j] <- sum(first.infection.at.age.hot*births.2050)/sum(births.2050)

	burden.longfuture.baseline.rel[j] <- sum(first.infection.at.age*births.2100)/sum(births.2100)
	burden.longfuture.hotter.rel[j] <- sum(first.infection.at.age.hot*births.2100)/sum(births.2100)


	## Contrast fertility patterns of the extremes (demography and one temperature; not Eswatini, just too too small to be interesting looking)
	if (j==12 | j==15 | j==21 | j==22 ) {

		yupper <- max(c(births.2023, births.2050, births.2100))

 		plot(15:49, births.2023,xlab="Age", ylab="Births per year age (1000s)", type="l",ylim=c(0,yupper))
		points(15:49, births.2050,type="l",lty=2)
		points(15:49, births.2100,type="l",lty=3)
		title(u.ctry[j])

		points(15:49, yupper*first.infection.at.age/max(first.infection.at.age), type="l", col="grey",lwd=2)
		points(15:49, yupper*first.infection.at.age.hot/max(first.infection.at.age.hot), type="l", col="coral",lwd=2)

		#susceptible is more interpretable?
		#points(15:49, yupper*(1-staying.susceptible), type="l", col="grey",lwd=2)
		#points(15:49, yupper*(1-staying.susceptible.hot), type="l", col="coral",lwd=2)

		## what I really want is the average age of first infection
		
		axis(4, at=c(0,0.25,0.5,0.75,1)*yupper, lab=c(0,0.25,0.5,0.75,1),lwd=2)
	}


}
legend("topright",legend=c(2023,2050,2100), lty=c(1,2,3),bty="n")



par(mfrow=c(1,2))
plot(burden.now,burden.future.baseline, xlab="Current burden (1000s)", ylab="Future burden (1000s)",pch=19,log="xy")
points(burden.now,burden.longfuture.baseline,pch=1) 
points(burden.now,burden.future.hotter,pch=19,col="red") 
points(burden.now,burden.longfuture.hotter,pch=1,col="red") 
abline(0,1)
title("Number of births affected") ## TODO CHECK


plot(burden.now.rel,burden.future.baseline.rel, xlab="Current burden (1000s)", ylab="Future burden (1000s)",pch=19,log="xy", ylim=range(c(burden.future.baseline.rel, burden.future.hotter.rel),na.rm=T))
points(burden.now.rel,burden.longfuture.baseline.rel,pch=1) 
points(burden.now.rel,burden.future.hotter.rel,pch=19,col="red") 
points(burden.now.rel,burden.longfuture.hotter.rel,pch=1,col="red") 
abline(0,1)
title("Proportion of births affected") ## CHECK





o1<- order(burden.future.baseline/burden.now)

par(mfrow=c(1,1), mar=c(5,10,2,2))
plot(c(burden.future.baseline/burden.now)[o1],1:length(o1), ylab="", xlab="Relative future burden",pch=19,xlim=c(0.5,2), axes=FALSE)
axis(2, at=1:length(o1),lab=u.ctry[o1],las=2);axis(1)

points(c(burden.longfuture.baseline/burden.now)[o1],1:length(o1),pch=1) 
points(c(burden.future.hotter/burden.now)[o1],1:length(o1),pch=19,col="red") 
points(c(burden.longfuture.hotter/burden.now)[o1],1:length(o1),pch=1,col="red") 
abline(v=1)





o1<- order(burden.future.baseline/burden.now)

par(mfrow=c(1,1), mar=c(10,5,2,2))
plot(1:length(o1),c(burden.future.baseline/burden.now)[o1], xlab="", ylab="Relative future burden",pch=19,ylim=c(0.5,2.5), axes=FALSE)
axis(1, at=1:length(o1),lab=u.ctry[o1],las=2);axis(2)

points(1:length(o1),c(burden.longfuture.baseline/burden.now)[o1],pch=1) 
points(1:length(o1),c(burden.future.hotter/burden.now)[o1],pch=19,col="red") 
points(1:length(o1),c(burden.longfuture.hotter/burden.now)[o1],pch=1,col="red") 
abline(h=1)
title("Relative number of births affected") ## TODO CHECK




##tunisia and morocco demography improves; but tunisia temperature big effect makes worse
##tz and ivory coast, demography makes worse; ivory coast tmperarure  no effect, tz, makes better


### IT SHOULD INDEED BE THAT COUNTRIES ON LEFT HAND SIDE EXPERIENCE AN INCREASE (LOW BURDEN, PUSH HIGHER) AND RIGHT A DECERASE, OR SOME SUCH SWITCH 
### TRIPLE CHECK HOW THE FOI IS WORKING - DO SANITY CHECKS ON E.G., BURDEN AND MEDIAN AGE BIRTH

## Could play games like find the temperature that would result in the same risk as the demography


### THE DUBEY ET AL. TELLS US THAT FROM 35 degree celsius, IT SHOULD START GOING DOWN. 
##TODO refit using simple lmer to avoid issue with predicting countrie s



