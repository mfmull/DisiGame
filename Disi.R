spy = 60*60*24*365 #(*seconds per year*);
Tt = 55;#Planning horizon
src='CONF'#CONF UNCONF

#Empirical consumptions
qJemp = 184000000./spy;
qJDrink = 120000000/spy #// N;
qJPipe = 100000000/spy #// N;
qJDrinkhist = 20000000/spy #// N;
qJAg = qJemp - qJDrink;
qSemp = 1032000000/spy;
#(*hJemp=130;*)


#yearly cost per flow rate per m head: depends on unit cost of fuel.
betJ = (1*1000*9.81*1)/0.23*0.9/(36.4e6)*spy ;#http://www.indexmundi.com/facts/jordan/pump-price-for-diesel-fuel Ave 1996-2012
betS = (1*1000*9.81*1)/0.23*0.07/(36.4e6)*spy ;#Efficiency of Diesel pump:http://www.waterandenergyprogress.org/library/Kranz11a.pdf
betD=(1*1000*9.81*1)/0.69*0.15/3.6e6*spy;#Efficiency of Electrical pump (Disi pump paper); electricity price:http://www.nepco.com.jo/en/electricity_tariff_en.aspx

#Historical Use
use=read.csv("WaterUse.csv")


#hydrogeol
load(paste("PHI2015.rdata",sep=''))
PHI=-PHI
PHI=subset(PHI,select=-c(y,yo))[-1,]
PHI$lag=0:(nrow(PHI)-1)
phidS=function(t) PHI$SauDom[which(PHI$lag==t)]
phida=function(t) PHI$AgDom[which(PHI$lag==t)]
phidd=function(t) PHI$DomDom[which(PHI$lag==t)]
phiaa=function(t) PHI$AgAg[which(PHI$lag==t)]
phiaS=function(t) PHI$AgSau[which(PHI$lag==t)]
phiSS=function(t) PHI$SauSau[which(PHI$lag==t)]
SS<-SD<-SA<-AA<-DD<-AD<-NULL
for(i in 1:nrow(PHI)){
  for(j in 1:nrow(PHI)){
    if(i-j<0){
      SS=c(SS,0);SD=c(SD,0);SA=c(SA,0);AA=c(AA,0);DD=c(DD,0);AD=c(AD,0)
    }else if(i-j==0){
      SS=c(SS,phiSS(i-j)/2);SD=c(SD,phidS(i-j)/2);SA=c(SA,phiaS(i-j)/2);AA=c(AA,phiaa(i-j)/2);DD=c(DD,phidd(i-j)/2);AD=c(AD,phida(i-j)/2)
    }else{
      SS=c(SS,phiSS(i-j));SD=c(SD,phidS(i-j));SA=c(SA,phiaS(i-j));AA=c(AA,phiaa(i-j));DD=c(DD,phidd(i-j));AD=c(AD,phida(i-j))}}}
SS=t(matrix(SS,nrow=nrow(PHI)));SD=t(matrix(SD,nrow=nrow(PHI)));SA=t(matrix(SA,nrow=nrow(PHI)));AA=t(matrix(AA,nrow=nrow(PHI)));DD=t(matrix(DD,nrow=nrow(PHI)));AD=t(matrix(AD,nrow=nrow(PHI)))

IC=data.frame(Ag = 103,Dom = 140,Sau=13)





lift.ts=function(schedule, IC){
  t=nrow(schedule)
  DDi=DD[1:t,1:t];ADi=AD[1:t,1:t];SDi=SD[1:t,1:t];AAi=AA[1:t,1:t];SAi=SA[1:t,1:t];SSi=SS[1:t,1:t]
  schedule=schedule*1e6/spy
  Dom=IC$Dom+DDi%*%schedule$Dom+ADi%*%schedule$Ag+SDi%*%schedule$Sau
  Ag=IC$Ag+ADi%*%schedule$Dom+AAi%*%schedule$Ag+SAi%*%schedule$Sau
  Sau=IC$Sau+SDi%*%schedule$Dom+SAi%*%schedule$Ag+SSi%*%schedule$Sau
  return(data.frame(Ag,Dom,Sau))
}

lft=lift.ts(use,IC)
h=lft[nrow(lft),]
q=use[nrow(use),]*1e6/spy
alphJa=betJ*(h$Ag+0.5*q$Dom*phida(0)+q$Ag*phiaa(0)+0.5*q$Sau*phiaS(0)) #ommitted 0.5 comes from the derivative in the first order condition!
alphSa=betS*(h$Sau+0.5*q$Dom*phidS(0)+0.5*q$Ag*phiaS(0)+q$Sau*phiSS(0))




#lift function
lift=function(schedule,IC){
  #requires consumption schedule and initial lift conditions; both data frames with a, d, S
  #returns lift conditions one period after the last consumption of the schedule
  qJahEMP=function(t) return(schedule$Ag[t]*1e6/spy)
  qJdhEMP=function(t) return(schedule$Dom[t]*1e6/spy)
  qShEMP=function(t) return(schedule$Sau[t]*1e6/spy)
  t=1:nrow(schedule) 
  Ag=IC$Ag+sum(unlist(sapply(t,qJdhEMP))*unlist(sapply(rev(t),phida)))+sum(unlist(sapply(t,qJahEMP))*unlist(sapply(rev(t),phiaa)))+sum(unlist(sapply(t,qShEMP))*unlist(sapply(rev(t),phiaS)))
  Sau=IC$Sau+sum(unlist(sapply(t,qJdhEMP))*unlist(sapply(rev(t),phidS)))+sum(unlist(sapply(t,qJahEMP))*unlist(sapply(rev(t),phiaS)))+sum(unlist(sapply(t,qShEMP))*unlist(sapply(rev(t),phiSS)))
  Dom=IC$Dom+sum(unlist(sapply(t,qJdhEMP))*unlist(sapply(rev(t),phidd)))+sum(unlist(sapply(t,qJahEMP))*unlist(sapply(rev(t),phida)))+sum(unlist(sapply(t,qShEMP))*unlist(sapply(rev(t),phidS)))
  return(data.frame(Ag,Sau,Dom))
}


data.frame(matrix(unlist(t(sapply(1:nrow(use), function(t) return(lift(use[1:t,],IC))))),ncol=3))


#Spontaneous drawdown: Preexisting conditions and drinking water.
