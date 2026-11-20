#Predict arrival timing at delta entry in forecast year. Based on forecasted outmigration timing at each RST from
#inseason model and forecasted travel time from CJS model
print("Arrival Timing")

#Trib_Abun[3]=mean(Trib_Abun[1:2]);Trib_Abun[4]=mean(Trib_Abun[1:2])
pwks=9:45#for plotting Nov02-Jul12
Plot_Type="Timing"
#Plot_Type="TT"
#Plot_Type="Surv"

predRST_stats=array(dim=c(Ntribs+1,Ncalwks,3))
predDelta_stats=predRST_stats
predTis_stats=predRST_stats;predKni_stats=predRST_stats

Nsamps=1000

#if you want to look at
#load(fn_surv)#look at intercepts
#dsv2=as.data.frame(fit,pars=c("T_bReach","T_bTrib"));summary(dsv2)
#browser()

#Plot relationships between flow and TT############
if(Ntribs==4){
  par(mfrow=c(2,2),xaxs="i",mai=c(.8,.7,.2,.2), omi=c(0.25,0.25,0.5,0.25),cex.axis=1.1,cex.lab=1.3,cex.main=1.4)
} else {
  par(mfrow=c(2,4),xaxs="i",mai=c(.8,.7,.2,.2), omi=c(0.25,0.25,0.5,0.25),cex.axis=1.1,cex.lab=1.3,cex.main=1.4)
}
isz=max(which(Xsz<=stm$musz))

for(itrib in 1:Ntribs){
  
  if(RSTInTrib[itrib]==0){ #mainstem
    xv=Xvec/1000
  } else { #tributary
    xv=XvecT[,SurvRstIndex[itrib]]/1000
  }
  
  plot_TT=matrix(nrow=NsX,ncol=3)
  plot_TT_Tis=plot_TT;plot_TT_Kni=plot_TT;
  for(ix in 1:NsX){
    
    if(RSTInTrib[itrib]==0){ #mainstem
      plot_TT[ix,]=quantile(stm$TTForecastSz_rst[,isz,ix,SurvRstIndex[itrib]],probs=CI)
      if(itrib<=4){
        plot_TT_Tis[ix,]=quantile(stm$TTForecastSz_rst_tis[,isz,ix,SurvRstIndex[itrib]],probs=CI)
        plot_TT_Kni[ix,]=quantile(stm$TTForecastSz_rst_kni[,isz,ix,SurvRstIndex[itrib]],probs=CI)
      }
      
    } else {
      plot_TT[ix,]=quantile(stm$TTForecastTSz_rst[,isz,ix,SurvRstIndex[itrib]],probs=CI)
    }
  }
  
  ymax=max(plot_TT)
  plot(xv,plot_TT[,2],type='l',col="red",bty='n',xlab="Peak flow during release (kcfs)",ylab="",ylim=c(0,ymax),main=TribName[itrib])
  polygon(x = c(rev(xv), xv), y = c(rev(plot_TT[,1]), plot_TT[,3]), col = alpha("grey", 0.4), border = NA) 
  lines(x=c(SurvCov[itrib]/1000,SurvCov[itrib]/1000),y=c(0,ymax),lty=2)
  
  if(itrib<=4){
    lines(xv,plot_TT_Tis[,2],col="blue")
    lines(xv,plot_TT_Kni[,2],col="yellow")
    points(x=sac_Q_range,y=c(0.01,0.01),pch=19,col='blue')
    if(itrib==1) legend("bottomleft",legend=c("Tisdale","Knights","Delta"),lty=rep(1,3),lwd=rep(2,3),col=c("blue","yellow","red"),bty='n')
  } else if(itrib==5) {
    points(x=butte_Q_range,y=c(0.01,0.01),pch=19,col='blue')
  } else if(itrib==7) {
    points(x=feath_Q_range,y=c(0.01,0.01),pch=19,col='blue')
  }
}
mtext("Travel Time (days)",side=2, las=3,line=0,outer=T,cex=1.3,font=2)



## Plot relationship between fish size and TT
if(Ntribs==4){
  par(mfrow=c(2,2),xaxs="i",mai=c(.8,.7,.2,.2), omi=c(0.25,0.25,0.5,0.25),cex.axis=1.1,cex.lab=1.3,cex.main=1.4)
} else {
  par(mfrow=c(2,4),xaxs="i",mai=c(.8,.7,.2,.2), omi=c(0.25,0.25,0.5,0.25),cex.axis=1.1,cex.lab=1.3,cex.main=1.4)
}
for(itrib in 1:Ntribs){
  
  plot_TT=matrix(nrow=Nsz,ncol=3);plot_TT_Tis=plot_TT;plot_TT_Kni=plot_TT
  for(isz in 1:Nsz){
    
    if(RSTInTrib[itrib]==0){ #mainstem
      #ix=max(which(Xvec<=stm$mux)) 
      ix=which.min(abs(Xvec-stm$mux)) 
      plot_TT[isz,]=quantile(stm$TTForecastSz_rst[,isz,ix,SurvRstIndex[itrib]],probs=CI)
      if(itrib<=4){
        plot_TT_Tis[isz,]=quantile(stm$TTForecastSz_rst_tis[,isz,ix,SurvRstIndex[itrib]],probs=CI)
        plot_TT_Kni[isz,]=quantile(stm$TTForecastSz_rst_kni[,isz,ix,SurvRstIndex[itrib]],probs=CI)
      }
      
    } else {
      if(SurvRstIndex[itrib]==1){
        muxT=muxB
      } else {
        muxT=muxF
      }
      ix=max(which(XvecT[,SurvRstIndex[itrib]]<=muxT))
      plot_TT[isz,]=quantile(stm$TTForecastTSz_rst[,isz,ix,SurvRstIndex[itrib]],probs=CI)
    }
  }#isz
  
  ymax=max(plot_TT)
  plot(Xsz,plot_TT[,2],type='l',col="red",bty='n',xlab="Forklength (mm)",ylab="",ylim=c(0,ymax),main=TribName[itrib])
  polygon(x = c(rev(Xsz), Xsz), y = c(rev(plot_TT[,1]), plot_TT[,3]), col = alpha("grey", 0.4), border = NA)   
  
  irecs=which(Xsz<=Fl_mu[itrib])#For x-axis position, find the Xsz index closes to FL_mu, the pRun-weighted mean size 
  if(length(irecs)==0){isz_mu=1} else{isz_mu=max(irecs)}
  lines(x=c(Xsz[isz_mu],Xsz[isz_mu]),y=c(0,ymax),lty=2)
  
  if(itrib<=4){
    lines(Xsz,plot_TT_Tis[,2],col="blue")
    lines(Xsz,plot_TT_Kni[,2],col="yellow")
    points(x=sac_sz_range,y=c(0.01,0.01),pch=19,col='blue')
    if(itrib==1) legend("bottomleft",legend=c("Tisdale","Knights","Delta"),lty=rep(1,3),lwd=rep(2,3),col=c("blue","yellow","red"),bty='n')
  } else if(itrib==5) {
    points(x=butte_sz_range,y=c(0.01,0.01),pch=19,col='blue')
  } else if(itrib==7) {
    points(x=feath_sz_range,y=c(0.01,0.01),pch=19,col='blue')
  }
  
}#itrib
mtext("Travel Time (days)",side=2, las=3,line=0,outer=T,cex=1.3,font=2)



AggRST=matrix(data=0,nrow=Nsamps,ncol=Ncalwks); AggDelta=AggRST
AggRST_temp=AggRST;AggDelta_temp=AggRST
AggTis=AggRST;AggTis_temp=AggDelta_temp
AggKni=AggRST;AggKni_temp=AggDelta_temp

TT_stats=array(dim=c(Ntribs,Ncalwks,3))

rst_delta_surv=array(dim=c(Nsamps,Ntribs,Ncalwks));rd_surv_stats=array(dim=c(Ntribs,Ncalwks,3))
rst_Tis_surv=array(dim=c(Nsamps,4,Ncalwks));rst_Kni_surv=rst_Tis_surv

par(mfrow=c(2,4),xaxs="i",mai=c(.8,.7,.2,.2), omi=c(0.25,0.25,0.5,0.25),cex.axis=1.1,cex.lab=1.3,cex.main=1.4)

for(itrib in 1:Ntribs){ ##############################################################################
  
  fnm=paste0(InseasDir,DoTrib[itrib],"_null.Rdata")
  load(file=fnm) #loads fit object from inseason model
  dis=as.data.frame(fit, pars = c("For_cp"))#grab forecast variable
  
  for(iwk in 1:Ncalwks){
    
    cpRST_temp=dis[,iwk];grecs=which(is.finite(cpRST_temp)==T)
    cpRST=sample(x=cpRST_temp[grecs],size=Nsamps,replace=F)
    
    predRST_stats[itrib,iwk,]=as.double(quantile(cpRST,probs=CI,na.rm=T)) #;predRST_stats[itrib,iwk,2]=mean(cpRST,na.rm=T)
    
    #find the closest size in Xsz for TT forecast given median size in this week from plad (Sz)
    #This hets the best isz index for Xsz used in CJS model for this itrib-iwk
    usewk=iwk;if(iwk>Nwks) usewk=Nwks  #calendar week can exceed modelled week in BT-SPAS/PLAD. Assume size after last modelled week = last size in PLAD week
    #irecs=which(Xsz<=Sz[itrib,usewk])#Sz populated in DS_Surv.R
    #if(length(irecs)==0){ 
    #  isz=1 #Sz is smaller than lowest value of Xsz so set to lowest index (unlikely)
    #} else{
    #  isz=max(irecs) #For Sz>max(Xsz) this will set index to last element for Xsz
    #}
    isz=which.min(abs(Xsz-Sz[itrib,usewk])) 
    
    #Find the closest index for TT flow covariate in forecast
    if(RSTInTrib[itrib]==0){ #mainstem
      #useXvec=max(which(Xvec<=SurvCov[itrib])) #find the index closest to the forecasted maxflow a
      useXvec=which.min(abs(Xvec-SurvCov[itrib])) 
      
      px=stm$TTForecastSz_rst[,isz,useXvec,SurvRstIndex[itrib]]
      
      if(itrib<=4){
        TT_days_Tis=sample(x=stm$TTForecastSz_rst_tis[,isz,useXvec,SurvRstIndex[itrib]],size=Nsamps,replace=F)
        TT_Tis=round(TT_days_Tis/7,digits=0); Arrive_wk_Tis=iwk+TT_Tis 
        irecs=which(Arrive_wk_Tis>Ncalwks);Arrive_wk_Tis[irecs]=Ncalwks
        
        TT_days_Kni=sample(x=stm$TTForecastSz_rst_kni[,isz,useXvec,SurvRstIndex[itrib]],size=Nsamps,replace=F)
        TT_Kni=round(TT_days_Kni/7,digits=0); Arrive_wk_Kni=iwk+TT_Kni
        irecs=which(Arrive_wk_Kni>Ncalwks);Arrive_wk_Kni[irecs]=Ncalwks
        
      }
      
    } else {
      #useXvec=max(which(XvecT[,SurvRstIndex[itrib]]<=SurvCov[itrib]))
      useXvec=which.min(abs(XvecT[,SurvRstIndex[itrib]]-SurvCov[itrib])) 
      px=stm$TTForecastTSz_rst[,isz,useXvec,SurvRstIndex[itrib]]
    }
    
    #icol=which(names(dsv)==Vnm)
    #TT_days=sample(x=dsv[,icol],size=Nsamps,replace=F)
    TT_days=sample(x=px,size=Nsamps,replace=F)
    TT_stats[itrib,iwk,]=as.double(quantile(TT_days,probs=CI))
    
    #The week of arrival in Delta is just week at RST + travel time from RST to Delta    
    TT=round(TT_days/7,digits=0)#convert from TT in units of days (CJS model) to TT in units of weeks
    Arrive_wk=iwk+TT #have to round as it points to a weekly element in tpA below
    irecs=which(Arrive_wk>Ncalwks);Arrive_wk[irecs]=Ncalwks
    
    if(iwk==1){
      pD=matrix(data=0,nrow=Nsamps,ncol=Ncalwks)#non-cumulative proportion at RST by week
      pA=matrix(data=0,nrow=Nsamps,ncol=Ncalwks)#non-cumulative proportion arriving in Delta by week
      N_pA=pA; cpA=pA
      
      pA_Tis=pA;N_pA_Tis=pA; cpA_Tis=pA
      pA_Kni=pA;N_pA_Kni=pA; cpA_Kni=pA
      
      pD[,iwk]=cpRST#non-cumulative probability of proportion past RST this week. For first week it is the same as cumulative probability
      
    } else {
      pD[,iwk]=cpRST-cpRST_last #cum prob for current week - cu prob from previous week = proportion in current week
    }
    cpRST_last=cpRST
    
    #Sample for the poterior distribution of RST-Delta survival rates for itrib for week wik. This depends on fish size for this week
    rst_delta_surv[,itrib,iwk] = sample(x=SurvPost[,isz,itrib],size=Nsamps,replace=F)
    rd_surv_stats[itrib,iwk,1:3]=quantile(rst_delta_surv[,itrib,iwk],probs=CI)
    
    #Combined RST arrival probabilities for each iwk into predicted arrival weeks. This creates the x-axis shift
    #the probabilities being accumulated are probability of passing RST * probability of surviving to Delta. Thus survival rate is reflected
    #in the contribution. Multiplying by Trib_Abun predicts the number arriving in Delta by week, needed for trib-aggregated computation
    for(isim in 1:Nsamps)N_pA[isim,Arrive_wk[isim]]=N_pA[isim,Arrive_wk[isim]]+Trib_Abun[itrib]*pD[isim,iwk]*mean(rst_delta_surv[,itrib,iwk])
    
    if(itrib<=4){
      rst_Tis_surv[,itrib,iwk] = sample(x=SurvPost_Tis[,isz,itrib],size=Nsamps,replace=F)
      rst_Kni_surv[,itrib,iwk] = sample(x=SurvPost_Kni[,isz,itrib],size=Nsamps,replace=F)
      
      for(isim in 1:Nsamps){
        N_pA_Tis[isim,Arrive_wk_Tis[isim]]=N_pA_Tis[isim,Arrive_wk_Tis[isim]]+Trib_Abun[itrib]*pD[isim,iwk]*mean(rst_Tis_surv[,itrib,iwk])
        N_pA_Kni[isim,Arrive_wk_Kni[isim]]=N_pA_Kni[isim,Arrive_wk_Kni[isim]]+Trib_Abun[itrib]*pD[isim,iwk]*mean(rst_Kni_surv[,itrib,iwk])
      }
    }
  }#next iwk
  
  for(isim in 1:Nsamps){
    for(iwk in 1:Ncalwks){
      #Standardize pA (proportion arriving by week at Delta) and then create cumulative probabilities
      pA[isim,iwk]=N_pA[isim,iwk]/sum(N_pA[isim,1:Ncalwks])  #so the sum of pAs across 1:Nwks=1 for each isim
      cpA[isim,iwk]=sum(pA[isim,1:iwk])#create cummulative distribution
      
      #Computations for across-tributary aggregated RST and Delta arrival timings weighted by Trib_wgt     
      #for AggDelta aso weight also weight by tributary-weekly survival rate
      AggRST_temp[isim,iwk]=AggRST_temp[isim,iwk]+Trib_Abun[itrib]*pD[isim,iwk]
      AggDelta_temp[isim,iwk]=AggDelta_temp[isim,iwk]+N_pA[isim,iwk]
      
      if(itrib<=4){
        pA_Tis[isim,iwk]=N_pA_Tis[isim,iwk]/sum(N_pA_Tis[isim,1:Ncalwks]) 
        cpA_Tis[isim,iwk]=sum(pA_Tis[isim,1:iwk])
        AggTis_temp[isim,iwk]=AggTis_temp[isim,iwk]+N_pA_Tis[isim,iwk]
        
        pA_Kni[isim,iwk]=N_pA_Kni[isim,iwk]/sum(N_pA_Kni[isim,1:Ncalwks]) 
        cpA_Kni[isim,iwk]=sum(pA_Kni[isim,1:iwk])
        AggKni_temp[isim,iwk]=AggKni_temp[isim,iwk]+N_pA_Kni[isim,iwk]
      }
      
    }
  }
  
  for(iwk in 1:Ncalwks) predDelta_stats[itrib,iwk,]=as.double(quantile(cpA[,iwk],probs=CI))    
  if(itrib<=4){
    for(iwk in 1:Ncalwks){
      predTis_stats[itrib,iwk,]=as.double(quantile(cpA_Tis[,iwk],probs=CI))    
      predKni_stats[itrib,iwk,]=as.double(quantile(cpA_Kni[,iwk],probs=CI))    
    }
  }
  
  
  if(Plot_Type=="Timing"){
    
    plot(pwks,predRST_stats[itrib,pwks,2],xlim=c(min(pwks),max(pwks)),axes=F,lwd=2,type='l',bty='n',main=paste0(TribName[itrib]," (",round(Trib_Abun[itrib],digits=0),")"),col="black",xlab="",ylab="",ylim=c(0,1))
    axis(2);uwks=seq(min(pwks),max(pwks), by=4)
    #axis(1,at=pwks,labels=rep("",length(pwks)))
    axis(1,at=uwks,labels=calwk[uwks],las=3,cex.lab=0.5)
    polygon(x = c(rev(pwks), pwks), y = c(rev(predRST_stats[itrib,pwks,1]), predRST_stats[itrib,pwks,3]), col = alpha("grey", 0.4), border = NA)
    lines(pwks,predDelta_stats[itrib,pwks,2],lwd=2,col="red")
    polygon(x = c(rev(pwks), pwks), y = c(rev(predDelta_stats[itrib,pwks,1]), predDelta_stats[itrib,pwks,3]), col = alpha("red", 0.15), border = NA)
    
    lines(pwks,predTis_stats[itrib,pwks,2],lwd=2,col="blue")
    lines(pwks,predKni_stats[itrib,pwks,2],lwd=2,col="yellow")
    if(itrib==1) legend("bottomright",legend=c("RST","Tisdale","Knights","Delta"),lwd=rep(2,4),col=c("black","blue","yellow","red"),bty='n')
    
    #week when median cumulative proportion is 0.5 by interpolation
    wkt=max(which(predDelta_stats[itrib,,2]<0.5))+1
    xslope=predDelta_stats[itrib,wkt+1,2]-predDelta_stats[itrib,wkt,2]
    med_wk_arr=wkt+(0.5-predDelta_stats[itrib,wkt,2])*1/xslope
    lines(x=c(0,med_wk_arr),y=c(0.5,0.5),col="red",lty=3);lines(x=c(med_wk_arr,med_wk_arr),y=c(0,0.5),col="red",lty=3)
    
    wkt=max(which(predRST_stats[itrib,,2]<0.5))+1
    xslope=predRST_stats[itrib,wkt+1,2]-predRST_stats[itrib,wkt,2]
    med_wk_dep=wkt+(0.5-predRST_stats[itrib,wkt,2])*1/xslope
    
    print(c(DoTrib[itrib],calwk[round(med_wk_dep,digits=0)],calwk[round(med_wk_arr,digits=0)],round(med_wk_arr,digits=2)))
    
  } else if (Plot_Type=="TT"){
    
    plot(pwks,TT_stats[itrib,pwks,2],xlim=c(min(pwks),max(pwks)),axes=F,lwd=2,type='l',bty='n',main=TribName[itrib],col="black",xlab="",ylab="",ylim=c(0,max(TT_stats[itrib,,])))
    axis(2);uwks=seq(min(pwks),max(pwks), by=4);  axis(1,at=uwks,labels=calwk[uwks],las=3,cex.lab=0.5)
    polygon(x = c(rev(pwks), pwks), y = c(rev(TT_stats[itrib,pwks,1]), TT_stats[itrib,pwks,3]), col = alpha("grey", 0.4), border = NA)
    
  } else if (Plot_Type=="Surv"){
    
    plot(pwks, rd_surv_stats[itrib,pwks,2],xlim=c(min(pwks),max(pwks)),axes=F,lwd=2,type='l',bty='n',main=TribName[itrib],col="black",xlab="",ylab="",ylim=c(0,max(rd_surv_stats[itrib,,])))
    axis(2);uwks=seq(min(pwks),max(pwks), by=4);  axis(1,at=uwks,labels=calwk[uwks],las=3,cex.lab=0.5)
    polygon(x = c(rev(pwks), pwks), y = c(rev(rd_surv_stats[itrib,pwks,1]), rd_surv_stats[itrib,pwks,3]), col = alpha("grey", 0.4), border = NA)
    
  }
  
}#next itrib ##############################################################


#Calculate tributary-aggregated RST and Delta entry timing curves, create statistics, and plot
for(isim in 1:Nsamps){
  for(iwk in 1:Ncalwks){
    AggRST[isim,iwk]=sum(AggRST_temp[isim,1:iwk]/sum(AggRST_temp[isim,1:Ncalwks]))
    AggDelta[isim,iwk]=sum(AggDelta_temp[isim,1:iwk]/sum(AggDelta_temp[isim,1:Ncalwks]))
    
    AggTis[isim,iwk]=sum(AggTis_temp[isim,1:iwk]/sum(AggTis_temp[isim,1:Ncalwks]))
    AggKni[isim,iwk]=sum(AggKni_temp[isim,1:iwk]/sum(AggKni_temp[isim,1:Ncalwks]))
  }
}

itrib=Ntribs+1
for(iwk in 1:Ncalwks) {
  predRST_stats[itrib,iwk,]=as.double(quantile(AggRST[,iwk],probs=CI))   
  predDelta_stats[itrib,iwk,]=as.double(quantile(AggDelta[,iwk],probs=CI))   
  
  predTis_stats[itrib,iwk,]=as.double(quantile(AggTis[,iwk],probs=CI))   
  predKni_stats[itrib,iwk,]=as.double(quantile(AggKni[,iwk],probs=CI))   
}


if(Plot_Type=="Timing"){
  plot(pwks,predRST_stats[itrib,pwks,2],xlim=c(min(pwks),max(pwks)),axes=F,lwd=2,type='l',bty='n',main="Aggregate",col="black",xlab="",ylab="",ylim=c(0,1))
  uwks=seq(min(pwks),max(pwks), by=4)
  axis(2)
  axis(1,at=pwks,labels=rep("",length(pwks)))
  axis(1,at=uwks,labels=calwk[uwks],las=3,cex.lab=0.5)
  polygon(x = c(rev(pwks), pwks), y = c(rev(predRST_stats[itrib,pwks,1]), predRST_stats[itrib,pwks,3]), col = alpha("grey", 0.4), border = NA)
  lines(pwks,predDelta_stats[itrib,pwks,2],lwd=2,col="red")
  polygon(x = c(rev(pwks), pwks), y = c(rev(predDelta_stats[itrib,pwks,1]), predDelta_stats[itrib,pwks,3]), col = alpha("red", 0.15), border = NA)
  
  if(Ntribs==4){
    lines(pwks,predTis_stats[itrib,pwks,2],lwd=2,col="blue")
    lines(pwks,predKni_stats[itrib,pwks,2],lwd=2,col="yellow")
  }
  
  wkt=max(which(predDelta_stats[itrib,,2]<0.5))+1
  xslope=predDelta_stats[itrib,wkt+1,2]-predDelta_stats[itrib,wkt,2]
  med_wk_arr=wkt+(0.5-predDelta_stats[itrib,wkt,2])*1/xslope
  lines(x=c(0,med_wk_arr),y=c(0.5,0.5),col="red",lty=2);lines(x=c(med_wk_arr,med_wk_arr),y=c(0,0.5),col="red",lty=2)
  
  wkt=max(which(predRST_stats[itrib,,2]<0.5))+1
  xslope=predRST_stats[itrib,wkt+1,2]-predRST_stats[itrib,wkt,2]
  med_wk_dep=wkt+(0.5-predRST_stats[itrib,wkt,2])*1/xslope
  
  print(c("Aggregate",calwk[round(med_wk_dep,digits=0)],calwk[round(med_wk_arr,digits=0)],round(med_wk_arr,digits=2)))
  
  mtext("End Date of Model Week",side=1, las=1,line=-1,outer=T,cex=1.3,font=2)
  mtext("Cumulative Proportion of Outmigrant Run",side=2, las=3,line=0,outer=T,cex=1.3,font=2)
  
  
  #Compare Delta arrival timing for aggregate with Tisdale and KL
  #if(Ntribs==4){#Only valid if Butte, Feather/Yuba not included in forecast
  par(mfrow=c(2,1),xaxs="i",mai=c(.8,.7,.5,.5), omi=c(0.25,0.25,0.5,0.25),cex.axis=0.9)
  main_nm=c("tisdale","knights landing")
  for(i in 1:2){
    
    pstats=switch(i,predTis_stats,predKni_stats)
    
    plot(pwks,pstats[itrib,pwks,2],xlim=c(min(pwks),max(pwks)),axes=F,lwd=2,main=main_nm[i],type='l',bty='n',col="red",xlab="",ylab="",ylim=c(0,1))
    uwks=seq(min(pwks),max(pwks), by=4)
    axis(2)
    axis(1,at=pwks,labels=rep("",length(pwks)))
    axis(1,at=uwks,labels=calwk[uwks],las=3,cex.lab=0.5)
    
    polygon(x = c(rev(pwks), pwks), y = c(rev(pstats[itrib,pwks,1]), pstats[itrib,pwks,3]), col = alpha("red", 0.15), border = NA)
    
    wkt=max(which(predDelta_stats[itrib,,2]<0.5))+1
    xslope=pstats[itrib,wkt+1,2]-pstats[itrib,wkt,2]
    med_wk_arr=wkt+(0.5-pstats[itrib,wkt,2])*1/xslope
    lines(x=c(0,med_wk_arr),y=c(0.5,0.5),col="red",lty=2);lines(x=c(med_wk_arr,med_wk_arr),y=c(0,0.5),col="red",lty=2)
    
    wkt=max(which(predRST_stats[itrib,,2]<0.5))+1
    xslope=predRST_stats[itrib,wkt+1,2]-predRST_stats[itrib,wkt,2]
    med_wk_dep=wkt+(0.5-predRST_stats[itrib,wkt,2])*1/xslope
    print(c("Aggregate",calwk[round(med_wk_dep,digits=0)],calwk[round(med_wk_arr,digits=0)],round(med_wk_arr,digits=2)))
    
    fn_main=paste0(InseasDir,main_nm[i],"_null.Rdata")
    load(fn_main)
    dmain=as.data.frame(fit, pars = c("For_cp"))#grab forecast variable
    main_stats=matrix(nrow=Ncalwks,ncol=3)
    for(iwk in 1:Ncalwks) main_stats[iwk,]=as.double(quantile(dmain[,iwk],probs=CI,na.rm=T))
    lines(pwks,main_stats[pwks,2],lwd=2,col="black")
    polygon(x = c(rev(pwks), pwks), y = c(rev(main_stats[pwks,1]), main_stats[pwks,3]), col = alpha("gray", 0.4), border = NA)
    
    wkt=max(which(main_stats[,2]<0.5))+1
    xslope=main_stats[wkt+1,2]-main_stats[wkt,2]
    med_wk=wkt+(0.5-main_stats[wkt,2])*1/xslope
    lines(x=c(0,med_wk),y=c(0.5,0.5),col="black",lty=2);lines(x=c(med_wk,med_wk),y=c(0,0.5),col="black",lty=2)
    print(c(main_nm[i],round(med_wk,digits=2),calwk[round(med_wk,digits=0)]))
    
    legend("bottomright",legend=c(paste0("JPE @ ",main_nm[i]),paste0(main_nm[i]," RST")),lwd=c(2,2),col=c("red","black"),bty='n')
  }
  
  mtext("End Date of Model Week",side=1, las=1,line=0,outer=T,cex=1.3,font=2)
  mtext("Cumulative Proportion of Outmigrant Run",side=2, las=3,line=0,outer=T,cex=1.3,font=2)
  #}
} else if(Plot_Type=="TT"){
  mtext("End Date of Model Week",side=1, las=1,line=-1,outer=T,cex=1.3,font=2)
  mtext("Travel Time (days)",side=2, las=3,line=0,outer=T,cex=1.3,font=2)
  
} else if(Plot_Type=="Surv"){
  mtext("End Date of Model Week",side=1, las=1,line=-1,outer=T,cex=1.3,font=2)
  mtext("RST-Delta Survival Rate",side=2, las=3,line=0,outer=T,cex=1.3,font=2)
  
}