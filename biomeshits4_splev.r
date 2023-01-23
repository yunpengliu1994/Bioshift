rm(list=ls())
## 计算每个物种的温度和降水生态位 ----	
	clim=read.csv("vars.csv")[,c("ADCODE99","Lon","Lat","bio01","bio05","bio06","bio12","bio16","bio17")]#Worldclim from luoao
	library(data.table)
	dis=as.data.frame(fread("Spdat_an_isrm_splev.csv"))#种级物种分布，详见biomeshigts3.r
	rownames(dis)=dis[,1]
	dis=dis[,-1]
	
	library(parallel)
	no_cores <- detectCores() - 1
	mycl <- makePSOCKcluster(no_cores);
	xx0<-parLapply(cl=mycl, X=1:dim(dis)[2],function(i,dis,clim){#ref Liu et al. NEE 2020 & Jezkova T, Wiens JJ. 2016 Rates of change in climatic niches ...
		x=dis[,i]
		adcd=rownames(dis)[which(as.logical(x))]	
		n<-length(adcd)
		data1=clim[clim$ADCODE99%in%adcd,]
		
		# climate mean
		AGG <- sapply(data1[,c("Lon","Lat","bio01","bio05","bio06","bio12","bio16","bio17")],function(x)
			c(mean=mean(x),median=median(x),sd=ifelse(length(x)==1,0,sd(x)),quant=ifelse(length(x)<=10,max(x),quantile(x,0.90)),quant=ifelse(length(x)<=10,min(x),quantile(x,0.10)),r=max(x)-min(x)))
		x1<-data.frame(t(matrix(AGG)))
		names(x1)<-paste(rep(names(data1)[-1],each=6),
				c(".mean",".median",".sd",".90",".10",".r"),sep="")
				
		# niche breadth
		# WLNBT:The within-locality niche breadth for temperature (WLNBT) is the difference between the maximum temperature of the warmest month (Bio5) 
		    # and the minimum temperature of the coldest month (Bio6) for each locality.
		# sppWLNBT: The mean WLNBT across localities for each species was calculated. 
		# SNBT: The species niche breadth for temperature (SNBT) is the difference between the maximum Bio5 and the minimum Bio6 across all the species’ sampled localities.
		# The ratioT is the mean WLNBT divided by SNBT for each species, which estimate how much within-locality niche breadth
			# contributes to species niche breadth.			
		WLNBT <- (data1$bio05-data1$bio06)	
		SNBT <- (max(data1$bio05)-min(data1$bio06))
		ratioT <- WLNBT/SNBT
		sppWLNBT <- mean(WLNBT)
		sppratioT <-mean(ratioT)
		NBVT <- ifelse(n==1,0,sd(WLNBT))
		NPVT <- ifelse(n==1,0,sd((data1$bio05+data1$bio06)/2))

		WLNBP <- data1$bio16-data1$bio17 		
		SNBP <- max(data1$bio16)-min(data1$bio17) 
		ratioP <- WLNBP/SNBP
		sppWLNBP <- mean(WLNBP)
		sppratioP <-mean(ratioP)
		NBVP <- ifelse(n==1,0,sd(WLNBP))
		NPVP <- ifelse(n==1,0,sd((data1$bio16+data1$bio17)/2))
		
		x2<-cbind(SNBT,sppWLNBT,sppratioT,NBVT,NPVT,
				SNBP,sppWLNBP,sppratioP,NBVP,NPVP)
		x3<-data.frame(x2)
		x4<-cbind(n,x1,x3)
		rownames(x4)<-colnames(dis)[i]
		return(x4)	
	},dis,clim)
	
	xx1<-do.call(rbind,xx0)
	stopCluster(mycl)
	write.csv(xx1,"Climate_niche.csv")
	
	xx0<-parLapply(cl=mycl, X=1:dim(dis)[2],function(i,dis,clim){#ref Liu et al. NEE 2020 & Jezkova T, Wiens JJ. 2016 Rates of change in climatic niches ...
		x=dis[,i]
		adcd=rownames(dis)[which(as.logical(x))]	
		n<-length(adcd)
		data1=clim[clim$ADCODE99%in%adcd,]
		
		# climate mean
		AGG <- sapply(data1[,c("bio01","bio05","bio06","bio12","bio16","bio17")],function(x)
			c(quant=ifelse(length(x)<=10,max(x),quantile(x,0.90)),quant=ifelse(length(x)<=10,min(x),quantile(x,0.10))))
		x1<-data.frame(t(matrix(AGG)))
		names(x1)<-paste(rep(c("bio01","bio05","bio06","bio12","bio16","bio17"),each=2),
				c(".90",".10"),sep="")
		
		rownames(x1)<-colnames(dis)[i]
		return(x1)	
	},dis,clim)
	
	xx1<-do.call(rbind,xx0)
	clim.niche=as.data.frame(fread("Climate_niche.csv"))
	rownames(clim.niche)=clim.niche$V1
	clim.niche=clim.niche[,-1]
	clim.niche2=data.frame(xx1[match(rownames(clim.niche),rownames(xx1)),],
		clim.niche[,which(!colnames(clim.niche)%in%colnames(xx1))])
	write.csv(clim.niche2,"Climate_niche.csv")
	
	biodis=as.data.frame(fread("Splev.new2.biomat.csv"))[,-1]	
	rownames(biodis)=biodis$genus	
	biodis2=data.frame(biodis[,which(!colnames(biodis)%in%colnames(xx1))],xx1[match(rownames(biodis),rownames(xx1)),])
	write.csv(biodis2,"Splev.new2.biomat.csv")
	
## Estimate absolute rates of niche evolution for each species ---
	library(data.table)
	clim.niche=as.data.frame(fread("Climate_niche.csv"))
	vlist=c(paste(rep(c("bio01","bio05","bio06","bio12","bio16","bio17"),each=3),c(".mean",".90",".10"),sep=""),"SNBT","sppWLNBT","SNBP","sppWLNBP")
	status=clim.niche[,vlist]
	rownames(status)=clim.niche$V1	
	
	library(parallel)
	no_cores <- detectCores() - 1
	mycl <- makePSOCKcluster(no_cores); 	
	biomat0=parLapply(cl=mycl, X=1:100,function(i,status,vlist){
		require(ape)
		tre0=read.tree(paste("FromAo_Smith_100treesV2/",i,".tre",sep=""))			
		#tre0=get(load(paste("trees/",i,".Rdata",sep="")))#100 radom trees of polytomy resolved tree from Smith et al.
		tre=drop.tip(tre0, tip=subset(tre0$tip.label,!tre0$tip.label%in%rownames(status)))		
		status2=as.matrix(status[tre$tip.label,])
		ans=do.call(cbind,lapply(1:length(vlist),function(j){
			stat=status2[,vlist[j]]			
			ans=castor::asr_squared_change_parsimony(tre,stat,weighted = TRUE,check_input = TRUE)$ance
			return(ans)
		}))			
		colnames(ans)=paste("ans",vlist,sep=".")
		re=unique(tre$edge[,1])
		re=re[order(re,decreasing=F)]
		rownames(ans)=re
		
		pos=which(tre$edge[,2]<=Ntip(tre))
		genus=data.frame(X1=tre$edge[pos,1],X2=tre$edge[pos,2],genus=tre$tip.label[tre$edge[pos,2]],age=tre$edge.length[pos])		
		biomat=cbind(genus,status2[match(genus$genus,rownames(status2)),],ans[match(genus$X1,rownames(ans)),])				
		biomat2=biomat[order(biomat$genus),]
		return(biomat2)	
	},status,vlist)
	stopCluster(mycl)
	
	vlist.ans=paste("ans",vlist,sep=".")
	bioans=do.call(cbind,lapply(1:length(vlist.ans),function(j,biomat0,vlist.ans){	
		re=do.call(cbind,lapply(1:100,function(i,j,biomat0,vlist.ans){			
			return(biomat0[[i]][,vlist.ans[j]])					
			},j,biomat0,vlist.ans))
		fin=apply(re,1,mean)
		return(fin)
	},biomat0,vlist.ans))
	colnames(bioans)=vlist.ans	
	biomat=cbind(biomat0[[1]][,c("genus","age",vlist)],bioans)	
	
	rate=do.call(cbind,lapply(1:length(vlist),function(i){
		rate=(biomat[,vlist[i]]-biomat[,vlist.ans[i]])/biomat$age
		return(rate)
	}))
	colnames(rate)=paste("rate",vlist,sep=".")
	rownames(rate)=biomat$genus
	
	biomat2=data.frame(biomat,rangesize=clim.niche[match(biomat$genus,clim.niche$V1),"n"],rate)	
	#write.csv(biomat2,"Splev.new2.biomat.csv")	#只有末端节点溶解了
	write.csv(biomat2,"Splev.new2.biomat.V2.csv")	#所有节点都溶解了
	
	##estimate niche evl rate use mocecular smith tree
	#Smith phylogenetic tree
	# ref:Smith, S. A., and J. W. Brown. 2018. Constructing a broadly inclusive seed plant phylogeny. American Journal of Botany 105(3): 1–13.
	# download from: https://github.com/FePhyFoFum/big_seed_plant_trees
	library(ape)
	tre.sm0=read.tree("GBOTB.tre")#mol tree	
	tre=drop.tip(tre.sm0, tip=subset(tre.sm0$tip.label,!tre.sm0$tip.label%in%rownames(status)))		
	status2=as.matrix(status[tre$tip.label,])
	ans=do.call(cbind,lapply(1:length(vlist),function(j){
		stat=status2[,vlist[j]]			
		ans=castor::asr_squared_change_parsimony(tre,stat,weighted = TRUE,check_input = TRUE)$ance
		return(ans)
	}))			
	colnames(ans)=paste("ans",vlist,sep=".")
	re=unique(tre$edge[,1])
	re=re[order(re,decreasing=F)]
	rownames(ans)=re
		
	pos=which(tre$edge[,2]<=Ntip(tre))
	genus=data.frame(X1=tre$edge[pos,1],X2=tre$edge[pos,2],genus=tre$tip.label[tre$edge[pos,2]],age=tre$edge.length[pos])		
	biomat=cbind(genus,status2[match(genus$genus,rownames(status2)),],ans[match(genus$X1,rownames(ans)),])	
	
	rate=do.call(cbind,lapply(1:length(vlist),function(i){
		rate=(biomat[,vlist[i]]-biomat[,paste("ans",vlist[i],sep=".")])/biomat$age
		return(rate)
	}))
	colnames(rate)=paste("rate",vlist,sep=".")
	rownames(rate)=biomat$genus
	
	biomat2=data.frame(biomat,rangesize=clim.niche[match(biomat$genus,clim.niche$V1),"n"],rate)
	biomat2=biomat2[order(biomat2$genus),]
	write.csv(biomat2,"Splev.new2.biomat.V2.smith.csv")#estimate niche evl rate use mocecular smith tree
	
##evaluate different evl Ref.Liu et al. NEE,2020,http://www.phytools.org/Cordoba2017/
	library(geiger)
	library(phytools)
	library(data.table)	
	clim.niche=as.data.frame(fread("Climate_niche.csv"))
	vlist=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")
	status=clim.niche[,vlist]
	rownames(status)=clim.niche$V1
	
	tre0=read.tree("FromAo_Smith_100treesV2/1.tre")		
	tre=drop.tip(tre0, tip=subset(tre0$tip.label,!tre0$tip.label%in%rownames(status)))		
	status2=as.matrix(status[tre$tip.label,])	

	evl.model=function(vlist.t,tre,status2){
		# calculate models of evolution
		WNfit1 <-fitContinuous (tre,status2[,vlist.t], model="white")
		BMfit1 <-fitContinuous (tre,status2[,vlist.t], model="BM")
		OUfit1 <-fitContinuous (tre,status2[,vlist.t], model="OU")
		LAfit1 <-fitContinuous (tre,status2[,vlist.t], model="lambda")
		EBfit1 <-fitContinuous (tre,status2[,vlist.t], model="EB")#maybe true		
		# compare individual models 
		AIC1 <- c (WNfit1$opt$aic, BMfit1$opt$aic,OUfit1$opt$aic,LAfit1$opt$aic,EBfit1$opt$aic) 
		AIC1 <-data.frame (AIC1)
		AIC1 <- t (AIC1)
		colnames (AIC1) <- c("WN", "BM","OU","LA","EB")
		return(data.frame(niche=vlist.t,AIC1))
	}		
	AIC=do.call(rbind,lapply(vlist,evl.model,tre,status2))
	write.csv(AIC,"AIC.csv")
	
	#using LA model to caculate niche evol rate
	library(parallel)
	no_cores <- detectCores() - 1
	mycl <- makePSOCKcluster(no_cores); 	
	biomat0=parLapply(cl=mycl, X=1:100,function(i,status,vlist){
		require(ape)
		require(geiger)
		require(phytools)
		tre0=read.tree(paste("FromAo_Smith_100treesV2/",i,".tre",sep=""))		
		tre=drop.tip(tre0, tip=subset(tre0$tip.label,!tre0$tip.label%in%rownames(status)))		
		status2=as.matrix(status[tre$tip.label,])
		ans=do.call(cbind,lapply(1:length(vlist),function(j){
			stat=status2[,vlist[j]]
			LAfit1 <-fitContinuous (tre,stat, model="lambda")
			# rescale the LA tree	
			LA1tree<-rescale(tre,model='lambda',lambda=LAfit1$opt$lambda)	
			ans=castor::asr_squared_change_parsimony(LA1tree,stat,weighted = TRUE,check_input = TRUE)$ance
			return(ans)
		}))			
		colnames(ans)=paste("ans",vlist,sep=".")
		re=unique(tre$edge[,1])
		re=re[order(re,decreasing=F)]
		rownames(ans)=re
		
		pos=which(tre$edge[,2]<=Ntip(tre))
		genus=data.frame(X1=tre$edge[pos,1],X2=tre$edge[pos,2],genus=tre$tip.label[tre$edge[pos,2]],age=tre$edge.length[pos])		
		biomat=cbind(genus,status2[match(genus$genus,rownames(status2)),],ans[match(genus$X1,rownames(ans)),])				
		biomat2=biomat[order(biomat$genus),]
		return(biomat2)	
	},status,vlist)
	stopCluster(mycl)
	
	vlist.ans=paste("ans",vlist,sep=".")
	bioans=do.call(cbind,lapply(1:length(vlist.ans),function(j,biomat0,vlist.ans){	
		re=do.call(cbind,lapply(1:100,function(i,j,biomat0,vlist.ans){			
			return(biomat0[[i]][,vlist.ans[j]])					
			},j,biomat0,vlist.ans))
		fin=apply(re,1,mean)
		return(fin)
	},biomat0,vlist.ans))
	colnames(bioans)=vlist.ans	
	biomat=cbind(biomat0[[1]][,c("genus","age",vlist)],bioans)	
	
	rate=do.call(cbind,lapply(1:length(vlist),function(i){
		rate=(biomat[,vlist[i]]-biomat[,vlist.ans[i]])/biomat$age
		return(rate)
	}))
	colnames(rate)=paste("rate",vlist,sep=".")
	rownames(rate)=biomat$genus
	
	biomat2=data.frame(biomat,rangesize=clim.niche[match(biomat$genus,clim.niche$V1),"n"],rate)		
	write.csv(biomat2,"Splev.new2.biomat.V2.LA.csv")	#所有节点都溶解了
	
## 比较不同类群的生态位进化速率 ----	
##1.age #rangesize
	library(ggpubr);library(scales)
	plotn=lapply(1:length(niche.rt),function(i){
		biodisp=biodis[,c("age","rangesize",niche.rt[i])]
		colnames(biodisp)[3]="vars"
		biodisp[,3]=abs(biodisp[,3])
		p <-ggplot(data=biodisp,aes(x=age,y=vars)) + 			
			scale_y_continuous(trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+
			#scale_x_continuous(trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+
			geom_point(size=0.6,colour='gray') +
			geom_density2d(show.legend=FALSE)+
			geom_smooth(method = "gam",colour="red",fill="darkgray",size=1)		
		 p+theme(axis.title.x = element_blank(),
			axis.title.y = element_text(size=15),
			axis.text.x  = element_text(size=15),
			axis.text.y  = element_text(size=15),
			panel.grid.major =element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black"),
			legend.position=c(0.8,0.2))+
			labs(y=niche.rt[i])
	})
	figure=ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],plotn[[5]],plotn[[6]],
	nrow=2,ncol =3,widths=rep(2,3),heights=rep(2,3),hjust=0,align="hv",common.legend=TRUE)	
	annotate_figure(figure, bottom = text_grob("age (Ma)", size=10))

##2.比较温带和热带类群;木本和草本 
	library(data.table)
	plant.order=as.data.frame(fread("GyAn.csv"))
	plant.fam=unique(as.data.frame(fread("SpLevDis2.csv"))[,c("Family_E","Genus_E","Species_E1")])
	taxonmic=cbind(Order=plant.order[match(plant.fam$Genus_E,plant.order$Genus),"Order"],plant.fam)
	write.csv(taxonmic,"taxonmic.csv")
	#科和目的信息
	taxonmic=as.data.frame(fread("taxonmic.csv",header=T))[,-1]	
	biodis0=as.data.frame(fread("Splev.new2.biomat.V2.csv"))[,-1]	
	rownames(biodis0)=biodis0$genus
	biodis=cbind(taxonmic[match(biodis0$genus,taxonmic$Species_E1),1:3],biodis0)
	# #划分热带和温带
	# geoname=read.csv("Geo-isrm.csv")
	# dis=as.data.frame(fread("SpLevDis2.csv"))
	# dis.lat=cbind(dis,geoname[match(dis$Adcode99,geoname$ADCODE99),c("Lon","Lat","continent")])
	# sp.lat=tapply(dis.lat$Lat,dis.lat$Species_E1,mean)
	# biodis2=cbind(biodis,lat=sp.lat[biodis$genus],to.tropic=ifelse(abs(sp.lat[biodis$genus])<=23.43677,"tro","tmp"))
	# #祖先重建
	# require(parallel)
	# sp2bio=biodis2$to.tropic
	# names(sp2bio)=biodis2$genus
	# no_cores <- detectCores() - 1
	# mycl <- makePSOCKcluster(no_cores);
	# node2bio0=lapply(1:100,function(i,sp2bio,mycl){
			# require(ape)
			# require(castor)
			# require(parallel)
			# tre_ori=get(load(paste("trees/",i,".Rdata",sep="")))			
			# d=subset(tre_ori$tip.label,!tre_ori$tip.label%in%names(sp2bio))
			# tre <- drop.tip(tre_ori, tip=d)		
			# status=sp2bio[tre$tip.label]		
			# stus0=map_to_state_space(status)
			# stus=stus0$mapped_states
			# names(stus)=names(status)		
			# ans=asr_max_parsimony(tre,stus)	
			# node2bio=ans$ancestral_likelihoods
			# colnames(node2bio)=names(stus0$name2index)
			# re=unique(tre$edge[,1])
			# re=re[order(re,decreasing=F)]
			# rownames(node2bio)=re
			# pos=which(tre$edge[,2]<=Ntip(tre))
			# node2bio.t=node2bio[match(unique(tre$edge[pos,1]),rownames(node2bio)),]
			# f2=function(x){
				# x=sort(x,decreasing=T)
				# ans.stat=ifelse(max(x)>=0.95,names(x[1]),NA)
			# }
			# node2bio2=parApply(cl=mycl,node2bio.t,1,f2)		
			# node2bio3=node2bio2[match(tre$edge[pos,1],names(node2bio2))]
			# names(node2bio3)=tre$tip.label[tre$edge[pos,2]]
			# return(node2bio3)
		# },sp2bio,mycl)
	# stopCluster(mycl)	
	
	# re=do.call(cbind,node2bio0)	
	# mf=function(x){
		# v=x[!is.na(x)]
		# uniqv=unique(v)
		# uniqv[which.max(tabulate(match(v,uniqv)))]	
	# }
	# node2bio=apply(re,1,mf)		
	# biodis3=cbind(biodis2,from.tropic=node2bio[biodis$genus])
	#lifeform
	life.form=as.data.frame(fread("GrowthForm_Angiosperms.csv"))
	library(stringr)
	Species_E1=str_squish(chartr(" × ","   ",as.character(life.form$Accepted_SpName1)))
	Species_E1=chartr(" ","_",Species_E1)	
	life.form2=data.frame(Species_E1=Species_E1,lifeform=ifelse(life.form$W==0,"H","W"))	
	biodis2=cbind(lifefrom=life.form2[match(biodis$genus,life.form2$Species_E1),-1],biodis)
	
	#APGIV clade的信息，同时依据apGIV更新了目	
	apg=read.csv("APGIV-taxonmic.csv")
	biodis3=cbind(apg[match(biodis2$Family_E,apg$Family_E),],biodis2[,-c(2,3)])		
	write.csv(biodis3,"Splev.new3.biomat.V2.csv")

	#t检验表 nich.rate.stat.csv
	library(data.table)
	biodis=as.data.frame(fread("Splev.new3.biomat.V2.csv"))[,-1]	
	niche=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")
	niche.rt=c("age","rangesize",paste("rate",niche,sep="."))	
	rate.ttest=function(measurevar,groupvar,biodis){
		require(Rmisc)
		#groupvar="to.tropic"
		#measurevar=niche.rt[3]
		biodisp=na.omit(biodis[,c(groupvar,"genus",measurevar)])		
		biodisp[,3]=abs(biodisp[,3])			
		tmp=summarySE(biodisp[,c(groupvar,measurevar)], measurevar=measurevar, groupvars=groupvar)
		test=function(biodisp){
			test=t.test(biodisp[,measurevar]~biodisp[,groupvar])
			P=test$p.value
			Px=ifelse(P <= 0.0001,"****",
							ifelse(P>0.0001&P<=0.001,"***",
							ifelse(P>0.001&P<=0.01,"**",
							ifelse(P>0.01&P<=0.05,"*",
							#ifelse(P>0.05&P<0.1,"",
							ifelse(P>0.05,"ns",P)))))
			t.test.p=paste(round(P,2),Px,sep="")
			return(t.test.p)		
		}		
		t.test.p=test(biodisp)
		re=cbind(niche=measurevar,tmp[,c(groupvar,"N")],niche.rate=paste(round(tmp[,measurevar],2),"±",round(tmp$se,2),sep=""),t.test.p)	
		return(re)
	}	
	#do.call(rbind,lapply(niche.rt,rate.ttest,"to.tropic",biodis))
	do.call(rbind,lapply(niche.rt,rate.ttest,"lifefrom",biodis))
	
	#作图
	library(ggplot2)
	biodis2=cbind(biodis,group=paste(biodis$from.tropic,biodis$to.tropic,sep="-"))		
	ggplot(na.omit(biodis2[,c("genus","rate.bio01.mean","lifefrom")]), aes(x=log(abs(rate.bio01.mean)), fill=lifefrom)) + geom_density(alpha=.3)
	ggplot(na.omit(biodis2[,c("genus","rangesize","rate.bio06.mean","from.tropic","to.tropic","group")]), aes(x=log(rangesize), fill=to.tropic)) + geom_density(alpha=.3)
	
	biodis3=biodis2[biodis2$group%in%c("tmp-tmp","tro-tmp"),]
	do.call(rbind,lapply(niche.rt,rate.ttest,"group",biodis3))
	ggplot(na.omit(biodis3[,c("genus","rate.bio12.mean","from.tropic","to.tropic","group")]), aes(x=log(abs(rate.bio12.mean)), fill=group)) + geom_density(alpha=.3)
	
## 进化速率相关性 ----
	library(data.table)	
	biodis=as.data.frame(fread("Splev.new3.biomat.V2.csv"))[,-1]	
	rownames(biodis)=biodis$genus
	var=paste(rep(c("bio01","bio05","bio06","bio12","bio16","bio17"),each=3),c("mean","90","10"),sep=".")
	niche.tp=c("MAT","MTWM","MTCM","MAP","MPWQ","MPDQ")
	niche=abs(biodis[,var])
	colnames(niche)=paste(rep(niche.tp,each=3),c("mean","90","10"),sep=".")
	library(ggcorrplot)	
	library(ggpubr)	
	cor.p=function(niche,tp){
		cord=niche[,c(paste(tp,"mean",sep="."),paste(tp,"90",sep="."),paste(tp,"10",sep="."))]
		corr <- round(cor(cord), 2)
		p.mat <- cor_pmat(cord)#sig
		ggcorrplot(corr, hc.order = FALSE, type = "lower", p.mat = p.mat,outline.color = "gray",insig = "blank",lab = TRUE)		
	}
	plotn=vector("list",6)	
	for (i in 1:6) plotn[[i]]=cor.p(niche,niche.tp[i])
	ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],plotn[[5]],plotn[[6]],
			nrow=2,ncol = 3,widths=c(2,2),heights=c(2,2,2),common.legend=TRUE,hjust=0,align="hv")

## 科的谱系信号 --- 
	library(data.table)
	taxonmic=unique(as.data.frame(fread("taxonmic.csv",header=T))[,c("Family_E","Genus_E","Species_E1")])	
	biodis0=as.data.frame(fread("Splev.new2.biomat.csv"))[,-1]	
	rownames(biodis0)=biodis0$genus
	biodis=cbind(taxonmic[match(biodis0$genus,taxonmic$Species_E1),1:2],biodis0)
	library(geiger)	
	tre=read.tree("ALLMB.tre")#与物种分布数据匹配1313807条，得到235249物种	
	faml=tapply(biodis$genus,biodis$Family_E,length)
	famlist=names(faml[faml>2])
	biodis.f=lapply(1:length(famlist),function(i){
		biodis[biodis$Family_E%in%famlist[i],]
	})
	#vlist0=paste(rep(c("bio01","bio05","bio06","bio12","bio16","bio17"),each=3),c(".mean",".90",".10"),sep="")
	vlist0=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")
	vlist=c(vlist0,paste("rate",vlist0,sep="."))
	signal=function(biodis.f){		
		#biodis.f=biodis[biodis$Family_E%in%"Metteniusaceae",]
		biodata=abs(biodis.f[,vlist.t])
		names(biodata)=biodis.f$genus		
		tre.t=ape::drop.tip(tre, tip=subset(tre$tip.label,!tre$tip.label%in%names(biodata)))
		treeWData <- try(geiger::treedata(tre.t, biodata, sort = T),silent=TRUE)
		temp <- try(phytools::phylosig(treeWData$phy, treeWData$data, 
										 method = "lambda", test = T),silent=TRUE);#Calculate Pagel's lambda
		Lambda=ifelse('try-error' %in% class(temp),NA,temp$lambda);
		Lambda.p=ifelse('try-error' %in% class(temp),NA,temp$P);
		temp2 <- try(phytools::phylosig(treeWData$phy, treeWData$data, 
										 method = "K", test = T),silent=TRUE); #Calculate Blomberg's k
		K <- ifelse('try-error' %in% class(temp2),NA,temp2$K);
		K.p <- ifelse('try-error' %in% class(temp2),NA,temp2$P);				
		sig=cbind(Lambda,Lambda.p,K,K.p)
			# re=rbind(re,sig)
			# }			
		return(sig)
	}
	require(parallel)	
	mycl <- makeCluster(120);	
	physig.all=famlist		
	for (j in 1:length(vlist)){
		vlist.t=vlist[j]
		clusterExport(cl = mycl, varlist = c("vlist.t","tre"))
		physig=parLapply(mycl,biodis.f,signal)
		physig.fam=do.call(rbind,physig)
		colnames(physig.fam)=paste(colnames(physig.fam),vlist.t,sep="_")
		physig.all=data.frame(physig.all,physig.fam)
	}	
	stopCluster(mycl)
	write.csv(physig.all,"physig.all.csv")
	date()	
	
## 进化速率格局图 ----
	library(data.table)	
	#biodis=as.data.frame(fread("Splev.new3.biomat.csv"))[,-1]
	biodis=as.data.frame(fread("Splev.new2.biomat.V2.smith.csv"))[,-1] #estimate niche evl rate use mocecular smith tree
	biodis=as.data.frame(fread("Splev.new3.biomat.V2.csv"))[,-1]#estimate niche evl rate BM
	biodis=as.data.frame(fread("Splev.new2.biomat.V2.LA.csv"))[,-1]#estimate niche evl rate LA
	rownames(biodis)=biodis$genus			
	dis=unique(as.data.frame(fread("SpLevDis2.csv"))[,c("Adcode99","Species_E1")])	
	niche=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")
	niche.rt=paste("rate",niche,sep=".")
		
	library(rgdal);library(ggpubr);library(dplyr)	
	map <- readOGR(dsn = "county_islandrm.dbf",stringsAsFactors=FALSE);map1<-fortify(map)
	cordacd<-as.data.frame(map@data)%>%mutate(ADCODE99=as.numeric(ADCODE99));
	cordacd$id<-as.character(0:(length(unique(map1$id))-1))
	
	# library(sp)
	# library(maptools)
	# shape<-readShapeSpatial("county_islandrm.shp")
	# shape2<-readShapeSpatial("coast_world/coast_world.shp")
	# library(foreign)
	# a=read.dbf("coast_world/Output.dbf")
	# shape3=subset(shape2,shape2@data$File_ID%in%a$File_ID)	
	# #raster file
	# library(raster)
	# data.test <- raster("extend.tree/rasters/0.00_temperature.asc",native=TRUE)
	# coord=coordinates(data.test)
	# rownames(coord)=paste(coord[,1],coord[,2],sep="_")	
	# cordacd <- read.csv("POINT2.csv")
	# cordacd=subset(cordacd,cordacd$ADCODE99>0)
	# rownames(cordacd)=paste(cordacd[,1],cordacd[,2],sep="_")	
	
	#plot	
	rate.p=function(biodis,rate,ratetype,dis,cordacd,map1,main){
		#rate=niche.rt[1];ratetype="inc";main="MAT"
		biodis=biodis[biodis$age>=1,]#182277 out of 231567	
		bios0=biodis[,rate]
		names(bios0)=rownames(biodis)
		if (rate%in%c("lifefrom","to.tropic","genus")&ratetype%in%"sprich"){				
			bios=bios0
			rate.dis=na.omit(cbind(dis,rate=bios[match(dis$Species_E1,names(bios))]))
			if (rate%in%"genus") {nich.rate=tapply(rate.dis$rate,rate.dis$Adcode99,length)} else {
				pro.f=function(x) length(x[x%in%"W"|x%in%"tro"])/(length(x[x%in%"W"|x%in%"tro"])+length(x[x%in%"H"|x%in%"tmp"]))
				nich.rate=tapply(rate.dis$rate,rate.dis$Adcode99,pro.f)
			}									
		} else{			
			if(ratetype=="inc") bios=bios0[bios0>0]#current>ans
			if(ratetype=="dec")	bios=abs(bios0[bios0<0])#current<ans
			if(ratetype=="all")	bios=abs(bios0)				
			rate.dis=na.omit(cbind(dis,rate=bios[match(dis$Species_E1,names(bios))]))		
			nich.rate=tapply(rate.dis$rate,rate.dis$Adcode99,mean)
		}	
		range01 <- function(x){(x-min(x))/(max(x)-min(x))}
			
		mdat=cbind(cordacd,nich.rate=nich.rate[match(cordacd$ADCODE99,names(nich.rate))])
		map2<-inner_join(map1,mdat)
			
		mdat.p=mdat[,c("Lat","nich.rate")]
		mdat.p$nich.rate2=183+range01(mdat.p$nich.rate)*40
		n=ifelse(rate%in%c("rate.bio01.mean","rate.bio05.mean","rate.bio06.mean"),1,0)	
		p<-ggplot()+
		geom_rect(aes(ymin=min(map1$lat)-5,ymax=max(map1$lat)+5,
					  xmin=min(map1$long),xmax=max(map1$long)),fill='#F2F7FB',col='gray75')+##E5F1FB
		geom_rect(aes(ymin=min(map1$lat)-5,ymax=max(map1$lat)+5,
					  xmin=max(map1$long),xmax=max(map1$long)+45),fill='transparent',col='gray75')+
		geom_polygon(data=map2,aes(x=long,y=lat,group=group,fill=nich.rate),color="transparent")+ 
		scale_fill_gradient2(midpoint=mean(nich.rate),low="blue",mid="yellow",high="red",name=main)+
		geom_point(data=mdat.p,aes(y=Lat,x=nich.rate2,fill=nich.rate),shape=21,size=1.5,pch=1,show.legend = F)+
		scale_y_continuous(expand = c(0,0),breaks = seq(-60,80,20),
						   labels = paste0(seq(-60,80,20),"°"),position = "right")+
		scale_x_continuous(expand = c(0,0),breaks = 183+c(14,28),
						   labels = c(round((max(nich.rate)-min(nich.rate))/3,n),round(2*(max(nich.rate)-min(nich.rate))/3,n)),
						   position = 'bottom')+
		theme(#aspect.ratio=(diff(range(test_df$y))+5)/(diff(range(test_df$x))+10*max(nich.rate)+2),#+75*max(srlist$P)
			  #text = element_text(size =10),
			  legend.background=element_rect(fill='transparent'),
			  panel.border = element_blank(),
			  axis.title = element_blank(),
			  axis.text = element_text(face="bold",size =10),
			  axis.text.x = element_text(angle=30),
			  axis.ticks = element_blank(),
			  panel.grid.minor = element_blank(),
			  legend.text=element_text(face="bold",size=10),
			  legend.title=element_text(face="bold",size=12),
			  legend.position = c(0.09,0.38)) 	
		return(p)
	}
	
	#all
	#bio01
	p1=rate.p(biodis,niche.rt[1],ratetype="all",dis,cordacd,map1,main="MAT")	
	p2=rate.p(biodis,niche.rt[1],ratetype="inc",dis,cordacd,map1,main="MAT_inc")
	p3=rate.p(biodis,niche.rt[1],ratetype="dec",dis,cordacd,map1,main="MAT_dec")
	#bio12
	p4=rate.p(biodis,niche.rt[4],ratetype="all",dis,cordacd,map1,main="MAP")
	p5=rate.p(biodis,niche.rt[4],ratetype="inc",dis,cordacd,map1,main="MAP_inc")
	p6=rate.p(biodis,niche.rt[4],ratetype="dec",dis,cordacd,map1,main="MAP_dec")
	ggarrange(p1,p4,p2,p5,p3,p6,
		nrow=3,ncol =2,widths=c(1,1,1),heights=c(1,1),labels=c("a","b","c","d","e","f"),font.label = list(size = 20))

	
	#bio05 Max Temperature of Warmest Month
	p1=rate.p(biodis,niche.rt[2],ratetype="all",dis,cordacd,map1,main="MTWM")
	#bio16 Precipitation of Wettest Quarter
	p2=rate.p(biodis,niche.rt[5],ratetype="all",dis,cordacd,map1,main="MPWQ")
	#bio06
	p3=rate.p(biodis,niche.rt[3],ratetype="all",dis,cordacd,map1,main="MTCM")
	#bio17
	p4=rate.p(biodis,niche.rt[6],ratetype="all",dis,cordacd,map1,main="MPDQ")
	ggarrange(p1,p2,p3,p4,
		nrow=2,ncol =2,widths=c(1,1),heights=c(1,1),labels=c("a","b","c","d"),font.label = list(size = 20))

	
	##plot niche rate of different lifeform
	#lifeform
	#190492 out of 231567 (82.26% with lifeform)
	biodis.h=biodis[biodis$lifefrom%in%"H",]#104353 out of 231567; 82152 out of 182325 
	biodis.w=biodis[biodis$lifefrom%in%"W",]#86139 out of 231567; 67987 out of 182325
	#herb
	#bio01
	p1=rate.p(biodis.h,niche.rt[1],ratetype="all",dis,cordacd,map1,main="MAT")
	#bio12
	p3=rate.p(biodis.h,niche.rt[4],ratetype="all",dis,cordacd,map1,main="MAP")
	#woody
	#bio01
	p2=rate.p(biodis.w,niche.rt[1],ratetype="all",dis,cordacd,map1,main="MAT")
	#bio12
	p4=rate.p(biodis.w,niche.rt[4],ratetype="all",dis,cordacd,map1,main="MAP")
	ggarrange(p1,p2,p3,p4,
		nrow=2,ncol =2,widths=c(1,1),heights=c(1,1),labels=c("a","b","c","d"),label.x=0.03,font.label = list(size = 20))

	## plot niche width and pos
	p1=rate.p(biodis,"SNBT",ratetype="all",dis,cordacd,map1,main="NichWid\n(Thermal)")
	p2=rate.p(biodis,"SNBP",ratetype="all",dis,cordacd,map1,main="NichWid\n(Hydrologic)")
	p3=rate.p(biodis,"bio01.mean",ratetype="all",dis,cordacd,map1,main="NichMean\n(Thermal)")
	p4=rate.p(biodis,"bio12.mean",ratetype="all",dis,cordacd,map1,main="NichMean\n(Hydrologic)")
	ggarrange(p1,p2,p3,p4,
		nrow=2,ncol =2,widths=c(1,1),heights=c(1,1),labels=c("a","b","c","d"),label.x=0.03,font.label = list(size = 20))

	##plot climate ---
	library(sp)
	library(maptools)
	shape<-readShapeSpatial("county_islandrm.shp")
	shape2<-readShapeSpatial("coast_world/coast_world.shp")
	library(foreign)
	a=read.dbf("coast_world/Output.dbf")
	shape3=subset(shape2,shape2@data$File_ID%in%a$File_ID)
	
	getColor = function(mapdata, provname, provcol, othercol){
			f = function(x, y) ifelse(x %in% y, which(y == x), 0);
			colIndex = sapply(mapdata@data$ADCODE99, f, provname);
			fg = c(othercol, provcol)[colIndex + 1];
			return(fg);
		}
	
	clim=read.csv("vars.csv")
	#bio04 = Temperature Seasonality (standard deviation ×100)
	#bio07 = Temperature Annual Range (BIO5-BIO6)
	#bio15 = Precipitation Seasonality (Coefficient of Variation)
	shape@data=cbind(shape@data,clim[match(shape@data$ADCODE99,clim$ADCODE99),c("bio04","bio07","bio15")])	
	plot.clim=function(shape,shape3,clim.var,main,lab){
		pop0=shape@data[,c("ADCODE99",clim.var)]				
		pop=pop0[pop0[,2]!=0,2];provname=as.character(pop0[pop0[,2]!=0,1])				
		col=data.frame(cwe=pop[order(pop)],col=colorRampPalette(c("palegreen3","sandybrown"))(length(pop)))
		provcol=as.character(col[match(pop,col$cwe),2])					
		
		width <- 60; height <- 30
		windows(width=width, height=height)	
		par(fig=c(0,1,0,1),new=FALSE)
		plot(shape, col = getColor(shape, provname, provcol, "gray"),border = getColor(shape, provname, provcol, "gray"))
		plot(shape3,col="black",add=T)
		box()
		mtext(lab,side=3,cex=2.5,col="black",adj=0.01,line=-2)
		
		#legend
		par(fig=c(0.075,0.2,0.05,0.55),new=TRUE)				
		barplot(as.matrix(rep(1,length(pop))),col=as.character(col$col),horiz=F,axes=F,border = NA)
		axis(2,c(seq(1,length(pop),round(length(pop)/4)),length(pop)),signif(sort(pop)[c(seq(1,length(pop),round(length(pop)/4)),length(pop))],2),cex.axis=1.3)
		title(main = list(main, cex = 1.5))
		box()
	}
	
	plot.clim(shape,shape3,clim.var="TSN",main="TSN","a")	
	plot.clim(shape,shape3,clim.var="ART",main="ART","b")
	plot.clim(shape,shape3,clim.var="PSN",main="PSN","c")
	
##不同类群进化速率的比较 ---
	#不同分类群genus richness
	clades=tapply(biodis$genus,biodis$clades,length)
					# ANA                 asterids     asterids-campanulids 
                     # 123                     8953                    10173 
        # asterids-lamiids           Caryophyllales            Chloranthanae 
                   # 31553                     7050                       65 
           # core eudicots              Dilleniales                 eudicots 
                      # 60                      342                     4648 
# eudicots-Ceratophyllanae               Magnoliids         mono-commelinids 
                       # 6                     7085                    17749 
                # monocots                   rosids               rosids-COM 
                   # 29753                      741                      810 
           # rosids-fabids       rosids-fabids-Nfix           rosids-malvids 
                     # 210                    39116                    20245 
              # Santalales            superasterids              superrosids 
                    # 1730                        3                     1862 				
	rate.p(biodis[biodis$clades%in%names(clades)[6],],niche.rt[1],ratetype="all",dis,shape,main=names(clades)[6])
	
	#不同clades进化速率的比较 bar plot with error bar	  
	library(ggplot2)
	#options("install.lock"=FALSE)
	library(Rmisc)	
	library(scales)
	library(data.table)	
	biodis=as.data.frame(fread("Splev.new3.biomat.V2.csv"))[,-1]	
	rownames(biodis)=biodis$genus	
	# vlist=paste("rate",c("bio01","bio05","bio06","bio12","bio16","bio17"),"mean",sep=".")
	# label.list=	c("MAT","MTWM","MTCM","MAP","MPWQ","MPDQ")
	vlist=paste("rate",c("bio01","bio12"),"mean",sep=".")
	label.list=	c("MAT","MAP")
	bardat=do.call(rbind,lapply(1:length(vlist),function(i){
		a=biodis[,c("genus",vlist[i],"clades")]	
		a_inc=subset(a,a[,vlist[i]]>0)
		a_dec=subset(a,a[,vlist[i]]<0)
		a[,vlist[i]]=abs(a[,vlist[i]])
		a_dec[,vlist[i]]=abs(a_dec[,vlist[i]])
		tmp=summarySE(a, measurevar=vlist[i], groupvars="clades")
		tmp_inc=summarySE(a_inc, measurevar=vlist[i], groupvars="clades")
		tmp_dec=summarySE(a_dec, measurevar=vlist[i], groupvars="clades")
		niche.stat=rbind(cbind(nichetype=label.list[i],tmp),cbind(nichetype=paste(label.list[i],"inc",sep="_"),tmp_inc),
			cbind(nichetype=paste(label.list[i],"dec",sep="_"),tmp_dec))
		niche.stat[is.na(niche.stat)]=0
		colnames(niche.stat)[4]="mean"		
		return(niche.stat)
	}))
 	
	clade.n=data.frame(
		clade=c("ANA","Chloranthanae","Magnoliids","mono-commelinids","monocots",
			"eudicots-Ceratophyllanae","eudicots","core eudicots","Dilleniales","superasterids","Caryophyllales","Santalales","asterids","asterids-campanulids",
			"asterids-lamiids","superrosids","rosids","rosids-malvids","rosids-fabids","rosids-COM","rosids-fabids-Nfix"),
		clade.n=1:21,
		labels=c("ANA", "Chl", "Mgn","MonC","Mon","DicC","Dic","DicCor","Dill","SpAst","Caryo","Santa","Ast","AstC","AstL",
			"SpRos","Ros","RosM","Fab","COM","Nfix"),
		col=c("#B51E90", "#E13A94", "#F588B8","#FEF201","brown","#8C52A2","#B384BA","#DEC3DE","pink","#00BA38","#9EBC90","green","#2ABAA0","#03B4F0","#D6F1FE",
			"black","darkgray","lightgray","#EC7862","#F48F74","#FBAA72"),id=rep(1,21)		
	)
scales::show_col(c("#B51E90", "#E13A94", "#F588B8","#FEF201","brown","#8C52A2","#B384BA","#DEC3DE","pink","#00BA38","#9EBC90","green","#2ABAA0","#03B4F0","#D6F1FE",
			"black","darkgray","lightgray","#EC7862","#F48F74","#FBAA72"))
	
	#plot the legend in a phylogeny
	## clade 和目画在谱系树上 ref:https://www.blog4xiang.world/posts/fb57811f.html
	library(ape)
	tre_ori=read.tree("ALLMB.tre")	
	apg=tapply(biodis$genus,biodis$clades,sample,1)
	tre=drop.tip(tre_ori, tip=subset(tre_ori$tip.label,!tre_ori$tip.label%in%apg))
	tre$tip.label=names(apg[match(tre$tip.label,apg)])		
	pos=which(tre$edge[,2]<=Ntip(tre))
	apg.edge=data.frame(node=tre$edge[pos,2],clade=tre$tip.label[tre$edge[pos,2]])		
	tre$tip.label=clade.n[match(tre$tip.label,clade.n$clade),"labels"]
	# if (!requireNamespace("BiocManager", quietly = TRUE))
    # install.packages("BiocManager")
	# BiocManager::install("ggtree")
	plot(tre,show.tip.label=T,edge.width=rep(3,Ntip(tre)),cex.main = 1.2, label.offset =10, no.margin = F, cex = 1)
	tiplabels(text=rep("",Ntip(tre)),tip=seq(1,Ntip(tre),1),cex=1,frame="circle",bg=clade.n[match(tre$tip.label,clade.n$labels),"col"])
	axisPhylo()
	
	bardat$clade.abbr=clade.n[match(bardat$clades,clade.n$clade),"labels"]
	bardat$clade.abbr=factor(bardat$clade.abbr,levels=clade.n[match(tre$tip.label,clade.n$labels),"labels"])	
	bardat$nichetype=factor(bardat$nichetype,levels=c("MAT","MAT_inc","MAT_dec","MAP","MAP_inc","MAP_dec"))	
	#scanter plot show nich evl rate ~ clade sprich
	p.scan=ggplot(bardat, aes(x=N, y=mean,color=clade.abbr),show.legend=FALSE) + 	
	scale_x_continuous(trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+
	scale_y_continuous(trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+
	geom_errorbar(aes(ymin=mean-se, ymax=mean+se),size=3,width=0.3,show.legend=FALSE)+
	scale_color_manual(values=clade.n[match(tre$tip.label,clade.n$labels),"col"])+		
	geom_point(aes(fill=clade.abbr),size=6,shape=21,color="black",show.legend=FALSE) + 
	scale_fill_manual(values=clade.n[match(tre$tip.label,clade.n$labels),"col"])+	
		xlab("Species richness") +ylab("Niche evolution rate") +theme_bw()+facet_wrap( ~ nichetype,nrow = 2,scales='free_y')+
	ggrepel::geom_text_repel(aes(x=N, y=mean,label=clade.abbr,fontface = "bold.italic"),size=6,color="black",show.legend=FALSE)+
	theme(axis.text = element_text(size=20,color='black'),		
		axis.title.x = element_text(size=25,color='black',angle=0),
		axis.title.y = element_text(size=25,color='black',angle=90),
		#panel.grid.major.y=element_blank(),
		panel.grid.minor.y=element_blank(),
		axis.line.y=element_line(linetype=1,color='black'),
		axis.line.x=element_line(linetype=1,color='black'),
		axis.ticks = element_line(linetype=2,color='black'),
		panel.grid=element_line(linetype=2,color='grey'),
		strip.text = element_text(size = 25),# 设置分面的字字体大小、颜色、背景、ref:https://www.jianshu.com/p/7c917d821641
		strip.background = element_rect(fill = "white", colour = "black")) 
	
	#bar plot with error bar	
	p.rate=ggplot(bardat, aes(x=clade.abbr, y=mean,color=clade.abbr),show.legend=FALSE) + 			
	scale_y_continuous(trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+	
	geom_errorbar(aes(ymin=mean-se, ymax=mean+se),size=1.8,width=1,show.legend=FALSE)+
	scale_color_manual(values=clade.n[match(tre$tip.label,clade.n$labels),"col"])+		
	geom_point(aes(fill=clade.abbr),size=5,shape=21,color="black",show.legend=FALSE) + 
	scale_fill_manual(values=clade.n[match(tre$tip.label,clade.n$labels),"col"])+	
	geom_vline(xintercept=3.5,col="red",size=1,linetype="longdash",alpha=0.5)+	
	geom_vline(xintercept=5.5,col="red",size=1,linetype="longdash",alpha=0.5)+	
	geom_vline(xintercept=9.5,col="red",size=1,linetype="longdash",alpha=0.5)+
	geom_vline(xintercept=15.5,col="red",size=1,linetype="longdash",alpha=0.5)+	
	ylab("Niche evolution rate") +theme_bw()+coord_flip()+facet_wrap( ~ nichetype,ncol = 6,scales='free_x')+
	#theme(strip.text= element_blank())+#theme_minimal()+
	theme(axis.text = element_text(size=25,color='black'),		
		axis.title.x = element_text(size=25,color='black',angle=0),
		axis.title.y = element_blank(),
		#panel.grid.major.y=element_blank(),
		panel.grid.minor.y=element_blank(),
		axis.line.y=element_line(linetype=1,color='black'),
		axis.line.x=element_line(linetype=1,color='black'),
		axis.ticks = element_line(linetype=2,color='black'),
		panel.grid=element_line(linetype=2,color='grey'),
		strip.text = element_text(size = 25),# 设置分面的字字体大小、颜色、背景、ref:https://www.jianshu.com/p/7c917d821641
		strip.background = element_rect(fill = "white", colour = "black")) 
	
	#不同目进化速率与原始程度的比较	
	library(data.table)	
	biodis=as.data.frame(fread("Splev.new3.biomat.V2.csv"))[,-1]	
#match species with order
	tpl=read.csv("Genus_synonym.csv")#from Zhiheng 2022.11.1
	biodis$Order=tpl[match(biodis$Genus_E,tpl$Genus_Accepted_E1),"Order_Soltis_Genbank"]
	unmat=subset(biodis,is.na(biodis$Order))
	unmat$Order=tpl[match(unmat$Genus_E,tpl$Genus_E1),"Order_Soltis_Genbank"]
	mat=subset(biodis,!is.na(biodis$Order))
	biodis=rbind(mat,unmat)
	write.csv(biodis, "Splev.new3.biomat.V2.csv")
	
	library(data.table)	
	biodis=as.data.frame(fread("Splev.new3.biomat.V2.csv"))[,-1]	
	rownames(biodis)=biodis$genus		
	# niche=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")
	# niche.rt=paste("rate",niche,sep=".")
	# vlist=niche.rt
	vlist=c("rate.bio01.mean","rate.bio12.mean")
	bardat=do.call(rbind,lapply(1:length(vlist),function(i){
		require(Rmisc)
		a=biodis[,c("genus",vlist[i],"Order")]		
		a.inc=subset(a,a[,vlist[i]]>0)
		a.dec=subset(a,a[,vlist[i]]<0)
		a[,vlist[i]]=abs(a[,vlist[i]])			
		tmp=na.omit(summarySE(a, measurevar=vlist[i], groupvars="Order"))
		colnames(tmp)[3]="mean"	
		tmp.inc=na.omit(summarySE(a.inc, measurevar=vlist[i], groupvars="Order"))
		colnames(tmp.inc)[3]="mean"	
		tmp.dec=na.omit(summarySE(a.dec, measurevar=vlist[i], groupvars="Order"))
		colnames(tmp.dec)[3]="mean"	
		tmp.dec$mean=abs(tmp.dec$mean)
		re=rbind(cbind(nichetype=vlist[i],tmp),
			cbind(nichetype=paste(vlist[i],"inc",sep="."),tmp.inc),
			cbind(nichetype=paste(vlist[i],"dec",sep="."),tmp.dec))
		return(re)
	}))
	
	root.dist <- function(tree, method = c("n.node","branch.length"), type=c("tip","node","both"))
	{
	if (length(method) > 1) method <- method[1]
	if (length(type) > 1) type <- type[1]
	
	if (class(tree) != "phylo") stop("The tree is not a phylo tree")
	
	## find the root
	root.label <- unique(tree$edge[,1][!tree$edge[,1] %in% tree$edge[,2]])
	
	tree.edge <- rbind(tree$edge,c(0,root.label))
	ii <- order(tree.edge[,2])
	tree.edge <- tree.edge[ii,]
	
	## the results of root distance
	N.tip <- Ntip(tree)
	N.node <- tree$Nnode
	N.edge <- Nedge(tree)
	if (type == "tip") {
		N.ii <- 1:N.tip
		rd <- numeric(N.tip)
		names(rd) <- tree$tip.label
		}
	else if (type =="node") {
		N.ii <- N.tip + (1:N.node)
		rd <- numeric(N.node)
		names(rd) <- tree$node.label
		}
	else if (type == "both") {
		N.ii <- 1:(N.tip + N.node)
		rd <- numeric(N.tip + N.node)
		names(rd) <- c(tree$tip.label, tree$node.label)
		}
	
	## count the number of nodes from root to tips
	if (method == "n.node")
		{
		for (i in 1:length(N.ii))
			{
			tip <- descendant <- N.ii[i]
			for (j in 1:N.node)
				{
				if (tree.edge[descendant,1] == 0)
					{
					rd[i] <- 0
					break
					}
				if (tree.edge[descendant,1] == root.label) 
					{
					rd[i] <- j
					break
					}
				ancestor <- tree.edge[descendant,1]
				descendant <- tree.edge[ancestor,2]
				}
			}
		}
	else if (method == "branch.length")
		{
		if (is.null(tree$edge.length)) edge.length <- rep(1, times = N.tip + N.node)
		else	edge.length <- c(tree$edge.length, 0)[ii]
		
		for (i in 1:length(N.ii))
			{
			tip <- descendant <- N.ii[i]
			for (j in 1:N.node)
				{
				rd[i] <- rd[i] + edge.length[descendant]
				if (tree.edge[descendant,1] == root.label | tree.edge[descendant,1] == 0) 
					{
					#rd[i] <- rd[i] + edge.length[root.label]
					break
					}
				ancestor <- tree.edge[descendant,1]
				descendant <- tree.edge[ancestor,2]
				}
			}		
		}
	return(rd)
	}
	ss.glm <- function(r.glm)
					{
					r.ss <- summary(r.glm)
					rsq <- 100*(r.ss$null.deviance-r.ss$deviance)/r.ss$null.deviance
					adj.rsq <- 100*(1-(r.ss$deviance/r.ss$df.residual)/(r.ss$null.deviance/r.ss$df.null))
					f.stat <- ((r.ss$null.deviance-r.ss$deviance)/(r.ss$df.null-
							r.ss$df.residual))/(r.ss$deviance/r.ss$df.residual)
					p <- pf(f.stat, r.ss$df.null-r.ss$df.residual, r.ss$df.residual, lower.tail=FALSE)
					return(c(r2=rsq,adj.r2=adj.rsq,p=p))
					}
	library(ape)
	tre0=read.tree("ALLMB.tre")	
	taxa=tapply(biodis$genus,biodis$Order,sample,1)
	tre=drop.tip(tre0, tip=subset(tre0$tip.label,!tre0$tip.label%in%taxa))
	root=root.dist(tre,"n.node","tip")
	root2=data.frame(order=names(taxa),rootdis=root[taxa])
	rate.rootdis=cbind(bardat,rootdis=root2[match(bardat$Order,root2$order),"rootdis"])
	#去掉极端值Metteniusales（47个物种)
	rate.rootdis2=rate.rootdis[!rate.rootdis$Order%in%"Metteniusales",]
	lab=data.frame(niche.type=unique(rate.rootdis2$nichetype),label=c("MAT","MAT_inc","MAT_dec","MAP","MAP_inc","MAP_dec"))
	rate.rootdis2$nichetype=lab[match(rate.rootdis2$nichetype,lab$niche.type),"label"]
	rate.rootdis2$nichetype=factor(rate.rootdis2$nichetype,levels=c("MAT","MAT_inc","MAT_dec","MAP","MAP_inc","MAP_dec"))
	#bar plot with error bar
	require(ggpubr);require(scales)
	the=theme(axis.text.x = element_text(size=18,color='black'),	
		axis.text.y = element_text(size=18,color='black',angle=60),	
		axis.title.x = element_text(size=28,color='black',angle=0),
		axis.title.y = element_text(size=28,color='black',angle=90),
		#panel.grid.major.y=element_blank(),
		panel.grid.minor.y=element_blank(),
		axis.line.y=element_line(linetype=1,color='black'),
		axis.line.x=element_line(linetype=1,color='black'),
		axis.ticks = element_line(linetype=2,color='black'),
		panel.grid=element_line(linetype=2,color='grey'),
		strip.text = element_text(size = 25),# 设置分面的字字体大小、颜色、背景、ref:https://www.jianshu.com/p/7c917d821641
		strip.background = element_rect(fill = "white", colour = "black")) 
	
	ggplot(rate.rootdis2, aes(x=rootdis, y=mean),show.legend=FALSE) + 	
	#scale_x_continuous(trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+
	geom_errorbar(aes(ymin=mean-se, ymax=mean+se),size=0.8,width=.4)+	
	geom_point(size=2,show.legend=FALSE) + 
	#scale_y_continuous(trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+
	xlab("Degree of divergences") +ylab("Niche evolutionary rate")+
	geom_smooth(method = "gam",colour="red",fill="#00BD5F",size=1)+
	theme_bw()+facet_wrap( ~ nichetype,nrow = 2,scales='free_y')+the
	
	ggplot(rate.rootdis2, aes(x=N, y=mean),show.legend=FALSE) + 	
	#scale_x_continuous(trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+
	geom_errorbar(aes(ymin=mean-se, ymax=mean+se),size=0.8,width=.4)+	
	geom_point(size=2,show.legend=FALSE) + 
	#scale_y_continuous(trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+
	xlab("Species richness within each order") +ylab("Niche evolutionary rate")+
	geom_smooth(method = "gam",colour="red",fill="#00BD5F",size=1)+
	theme_bw()+facet_wrap( ~ nichetype,nrow = 2,scales='free')+the
	
	
##on grid cell level
	library(data.table)
	dis=as.data.frame(fread("SpLevDis2.csv"))[,-1]	
	biodis=as.data.frame(fread("Splev.new3.biomat.V2.csv"))[,-1]	
	rownames(biodis)=biodis$genus
	biodis=biodis[biodis$age>=1,]#182277 out of 231567
	
	# load("biostat.fam.Rdata")
	# biostat.fam=biostat.fam[biostat.fam$age>=1,]
	niche=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")
	#niche.lmd=paste("lmd",niche,sep=".")
	#niche.K=paste("K",niche,sep=".")		
	niche.rt=paste("rate",niche,sep=".")	
	dis0=cbind(dis,#biostat.fam[match(dis$Family_E,rownames(biostat.fam)),c(niche.K,"r0_mag","r0.9_mag","r0_bar","r0.9_bar")],
		biodis[match(dis$Species_E1,rownames(biodis)),
			c(niche.rt,"clades","lifefrom","rangesize","age","SNBT","sppWLNBT","SNBP","sppWLNBP",niche)])
	dis2=subset(dis0,!is.na(dis0$rate.bio01.mean))
	
	dis.w=subset(dis2,dis2$lifefrom%in%"W")
	dis.h=subset(dis2,dis2$lifefrom%in%"H")
	dis.ast=subset(dis2,dis2$clades%in%c("superasterids","asterids","asterids-campanulids","asterids-lamiids"))
	dis.mono=subset(dis2,dis2$clades%in%c("mono-commelinids","monocots"))
	dis.basal=subset(dis2,dis2$clades%in%c("ANA"))
	dis.Nfix=subset(dis2,dis2$clades%in%"rosids-fabids-Nfix")
	dis.mal=subset(dis2,dis2$clades%in%"rosids-malvids")
	
	get.dat.grid=function(dis2){
		mf=function(x) length(x[x%in%"W"])/(length(x[x%in%"W"])+length(x[x%in%"H"]))	
		woody.pro=tapply(dis2$lifefrom,dis2$Adcode99,mf)		
		mf2=function(x) mean(na.omit(x))
		mf3=function(x) mean(log(abs(na.omit(x))))
		mf3.dec=function(x) mean(log(abs(na.omit(x[x<0]))))
		mf3.inc=function(x) mean(log(abs(na.omit(x[x>0]))))	
		#species richness weighted mean phylosig and diverficate rate of each family	
		dat.grid=cbind(woody.pro,
			do.call(cbind,lapply(c("rangesize","age","SNBT","sppWLNBT","SNBP","sppWLNBP",niche),
				function(i){tapply(dis2[,i],dis2[,"Adcode99"],mf2)})),
			do.call(cbind,lapply(niche.rt,function(i){tapply(dis2[,i],dis2[,"Adcode99"],mf3)})),
			do.call(cbind,lapply(niche.rt,function(i){tapply(dis2[,i],dis2[,"Adcode99"],mf3.dec)})),
			do.call(cbind,lapply(niche.rt,function(i){tapply(dis2[,i],dis2[,"Adcode99"],mf3.inc)})))
		colnames(dat.grid)=c("WoodyPro","Range","Age","SNBT","sppWLNBT","SNBP","sppWLNBP",
			niche,niche.rt,paste(niche.rt,"dec",sep="_"),paste(niche.rt,"inc",sep="_"))
		dat.grid=data.frame(adcode=rownames(dat.grid),dat.grid)
		return(dat.grid)
	}		
	dat.grid=get.dat.grid(dis2)
	dat.grid=get.dat.grid(dis.w)
	dat.grid=get.dat.grid(dis.h)
	dat.grid=get.dat.grid(dis.ast)
	dat.grid=get.dat.grid(dis.mono)
	dat.grid=get.dat.grid(dis.basal)
	dat.grid=get.dat.grid(dis.Nfix)
	dat.grid=get.dat.grid(dis.mal)
	
	clim=read.csv("vars.csv")
	#bio04 = Temperature Seasonality (standard deviation ×100)
	#bio07 = Temperature Annual Range (BIO5-BIO6)
	#bio15 = Precipitation Seasonality (Coefficient of Variation)
	dat.grid2=cbind(dat.grid,clim[match(dat.grid$adcode,clim$ADCODE99),c("bio04","bio07","bio15")])
	
	#scanter plot grid-level
	getscater.grid=function(x.list,y.list,xlab,ylab,dat.grid){
		ss.glm <- function(r.glm){
			r.ss <- summary(r.glm)
			rsq <- 100*(r.ss$null.deviance-r.ss$deviance)/r.ss$null.deviance
			adj.rsq <- 100*(1-(r.ss$deviance/r.ss$df.residual)/(r.ss$null.deviance/r.ss$df.null))
			f.stat <- ((r.ss$null.deviance-r.ss$deviance)/(r.ss$df.null-
				r.ss$df.residual))/(r.ss$deviance/r.ss$df.residual)
			p <- pf(f.stat, r.ss$df.null-r.ss$df.residual, r.ss$df.residual, lower.tail=FALSE)
			return(c(r2=rsq,adj.r2=adj.rsq,p=p))
		}		
		temp.x<-as.data.frame(dat.grid[,x.list])				
		temp.y<-as.data.frame(dat.grid[,y.list])			
		temp=cbind(temp.x,temp.y)	
		colnames(temp)=c(x.list,y.list)
		width <- 70; height <- 40
		windows(width=width, height=height)
		n.col <- length(x.list); n.row <- length(y.list)
		mylayout <- layout(matrix(1:((2+n.col)*(2+n.row)), ncol=2+n.col, byrow=T), 
			width=c(0.7*width/(2+n.col), rep(width/(2+n.col),times = n.col), 0.3*width/(2+n.col)),
			height=c(0.3*height/(2+n.row), rep(height/(2+n.row),times=n.row), 0.7*height/(2+n.row)))
		layout.show((2+n.col)*(2+n.row))			
		glmr=c()
		for (i in 1:length(y.list)) {		
			glmr.t2=c()
			if (i == 1) {
				for (j in 1:(2+n.col)) {par(mar=c(0,0,0,0)); plot.new()}
				}	
				for (xi in 1:length(x.list)) {
					temp2=na.omit(temp[,c(x.list[xi],y.list[i])])				
					x <- temp2[,1]
					y <- temp2[,2]			
					if (xi == 1) {par(mar=c(0,0,0,0)); plot.new()}
					par(mar=c(0.5,0.5,0,0), cex=1.5, cex.axis=1.3, cex.lab=1.5, mgp=c(2.4,0.3,0), tck=-0.04)
					plot(y~x, data=temp2, pch=19, xlab="", ylab="", cex=0.5, col='grey',axes=F); box()
					m <- glm(y~x,data=temp2)
					glmr.t=data.frame(x=xlab[xi],y=ylab[i],adj.r2=round(ss.glm(m)[2],1),p=round(ss.glm(m)[3],2))
					a <- seq(min(x), max(x), (max(x)-min(x))/10)
					if (ss.glm(m)[3] < 0.05)  lines(x=a, y=predict(m, newdata=list(x=a)), col="black", lty=1, lwd=0.8)
					if (xi == 1) {
						axis(side=2, label=T)
						mtext(side=2, text=ylab[i], line=1.6, cex=1.5, las=0)
					}
					if (i == length(y.list)) {
						axis(side=1, lwd = 0.5)
						mtext(side=1, text=xlab[xi], line=1.6, cex=1.5, las=0)
					}
					if (xi == length(x.list)) {par(mar=c(0,0,0,0)); plot.new()}
					glmr.t2=rbind(glmr.t2,glmr.t)	
				}
			glmr=rbind(glmr,glmr.t2)		
			}
			return(glmr)
		}
			
	x.list.t=c("SNBT","bio01.mean","bio04","WoodyPro")
	x.list.p=c("SNBP","bio12.mean","bio15","WoodyPro")
	xlab=c("NichWid","NichMean","ClimVar","WoodyPro")
	
	y.list.t=c("rate.bio01.mean","rate.bio01.mean_inc","rate.bio01.mean_dec")
	y.list.p=c("rate.bio12.mean","rate.bio12.mean_inc","rate.bio12.mean_dec")
	ylab.t=c("MAT","MAT_inc","MAT_dec")
	ylab.p=c("MAP","MAP_inc","MAP_dec")	
	
	getscater.grid(x.list.t,y.list.t,xlab,ylab.t,dat.grid2)
	getscater.grid(x.list.p,y.list.p,xlab,ylab.p,dat.grid2)
	
	y.list.div=c("r0_mag","r0.9_mag","r0_bar","r0.9_bar")
	getscater.grid(c(y.list.t,y.list.p),y.list.div,c(ylab.t,ylab.p),y.list.div,dat.grid,"Divsificate rate")
	
	#scanter plot species level
	niche=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")	
	niche.rt=paste("rate",niche,sep=".")	
	biodis.plt=biodis
	biodis.t.dec=subset(biodis.plt,biodis.plt$rate.bio01.mean<0)
	biodis.t.inc=subset(biodis.plt,biodis.plt$rate.bio01.mean>0)
	biodis.p.dec=subset(biodis.plt,biodis.plt$rate.bio12.mean<0)
	biodis.p.inc=subset(biodis.plt,biodis.plt$rate.bio12.mean>0)
	biodis.plt[,niche.rt]=log(abs(biodis.plt[,niche.rt]))
	biodis.t.dec[,niche.rt]=log(abs(biodis.t.dec[,niche.rt]))
	biodis.t.inc[,niche.rt]=log(abs(biodis.t.inc[,niche.rt]))
	biodis.p.dec[,niche.rt]=log(abs(biodis.p.dec[,niche.rt]))
	biodis.p.inc[,niche.rt]=log(abs(biodis.p.inc[,niche.rt]))
	getscater.grid(x.list.t,"rate.bio01.mean",xlab,"MAT",biodis.plt,"MAT rate")
	getscater.grid(x.list.t,"rate.bio01.mean",xlab,"MAT_dec",biodis.t.dec,"MAT dec rate")
	getscater.grid(x.list.t,"rate.bio01.mean",xlab,"MAT_inc",biodis.t.inc,"MAT inc rate")
	getscater.grid(x.list.p,"rate.bio12.mean",xlab,"MAP",biodis.plt,"MAP rate")
	getscater.grid(x.list.p,"rate.bio12.mean",xlab,"MAP_dec",biodis.p.dec,"MAP dec rate")
	getscater.grid(x.list.p,"rate.bio12.mean",xlab,"MAP_inc",biodis.p.inc,"MAP inc rate")
	#hp
	library(ggplot2)
	hp.dat=c()
	y.var=c(y.list.t,y.list.p)
	y.var.lab=c(ylab.t,ylab.p)
	for(i in 1:length(y.var)){
		if (i<=3) hp.data=dat.grid2[,c("WoodyPro","Range","SNBT","bio01.mean","bio04",y.var[i])] else 
			hp.data=dat.grid2[,c("WoodyPro","Range","SNBP","bio12.mean","bio15",y.var[i])]
		colnames(hp.data)=c("WoodyPro","Range","NichWid","NichMean","ClimVar","NichRate")
		hp=hier.part::hier.part(hp.data$NichRate,hp.data[,c("ClimVar","NichWid","NichMean","WoodyPro")],gof = "Rsqu",barplot = FALSE)$IJ
		hp.dat.t=rbind(data.frame(vars=rownames(hp),rsq=hp$I*100,rsq.type="Independent",hp.type=y.var.lab[i]),data.frame(vars=rownames(hp),rsq=hp$J*100,rsq.type="Joint",hp.type=y.var.lab[i]))		
		hp.dat=rbind(hp.dat,hp.dat.t)
	}
	hp.dat$hp.type=factor(hp.dat$hp.type,levels=y.var.lab)
	ggplot(hp.dat,aes(x=vars,y=rsq,fill=rsq.type))+ geom_bar(stat="identity",position = "stack")+ 
	geom_hline(yintercept=0,col="black",size=1,linetype="longdash",alpha=0.5)+
	facet_wrap(~hp.type)+theme_bw()+ylab("Rsq of hierarchical partitioning") +theme_bw()+
	theme(axis.text.x = element_text(size=12,color='black',angle=15),
		axis.text.y = element_text(size=12,color='black'),
		axis.title.y = element_text(size=15,color='black'),
		axis.title.x = element_blank(),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),		
		strip.text = element_text(size = 15),# 设置分面的字字体大小、颜色、背景、ref:https://www.jianshu.com/p/7c917d821641
		strip.background = element_rect(fill = "white", colour = "black")) 
	
	#残差与网格内气候变异的关系
	clim.var=read.csv("datall_clim.csv")
	hp.data=dat.grid2[,c("WoodyPro","Range","SNBT","SNBP","bio01.mean","bio12.mean","bio15","bio04","rate.bio01.mean","rate.bio12.mean")]
	hp.data2=cbind(hp.data,clim.var[match(rownames(hp.data),clim.var$ADCODE99),])
	M.t <- lm(rate.bio01.mean ~ bio04 + SNBT + bio01.mean + WoodyPro,data=hp.data2)
	E.t <- resid(M.t)
	
	M.p <- lm(rate.bio12.mean ~ bio15 + SNBP + bio12.mean + WoodyPro,data=hp.data2)
	E.p <- resid(M.p)

	par(mfrow = c(2,2),mar = c(5,5,2,2))
	plot(E.t ~ hp.data2$bio01_RANGE,xlab = "Range within GSU(MAT)",ylab = "Residuals",cex.lab = 1.5)
	abline(h = 0, lty = 2,col="red")
	plot(E.t ~ hp.data2$bio04_RANGE,xlab = "Range within GSU(TSN)",ylab = "Residuals",cex.lab = 1.5)
	abline(h = 0, lty = 2,col="red")
	plot(E.p ~ hp.data2$bio12_RANGE,xlab = "Range within GSU(MAP)",ylab = "Residuals",cex.lab = 1.5)
	abline(h = 0, lty = 2,col="red")
	plot(E.p ~ hp.data2$bio15_RANGE,xlab = "Range within GSU(PSN)",ylab = "Residuals",cex.lab = 1.5)
	abline(h = 0, lty = 2,col="red")

	#SEM
	library(lavaan)
	library(semPlot)
	library(magrittr)
	library(dplyr)
	par(mfrow = c(2,3),mar=c(0,0,0.4,0),oma=c(0,0,0,0))	
	#AIC=c()		
	for(i in 1:length(y.var)){
		if (i<=3) sem.data=dat.grid2[,c("WoodyPro","SNBT","bio01.mean","bio04",y.var[i])] else 
			sem.data=dat.grid2[,c("WoodyPro","SNBP","bio12.mean","bio15",y.var[i])]
		colnames(sem.data)=c("WoodyPro","NichWid","NichMean","ClimVar","NichRate")
		myModel1="NichRate~NichWid+ClimVar+	NichMean +	WoodyPro
			NichWid~ClimVar
			WoodyPro~ClimVar
			NichMean~~ClimVar"			
		fit1 <- sem(myModel1, data = scale(sem.data))
		fit=fit1
		# myModel2="NichRate~WoodyPro+NichMean		
			# NichWid~NichRate+WoodyPro+NichMean			
			# NichMean~~WoodyPro"			
		# fit2 <- sem(myModel2, data = scale(sem.data))			
		
		# #AIC值越小的模型表明越有可能准确地预测新数据
		# if (AIC(fit1)<AIC(fit2)) {fit=fit1;best="wid2rate"} else {fit=fit2;best="rate2wid"}
		# AIC.t=data.frame(SEM1=AIC(fit1),SEM2=AIC(fit2),Best=best)
		# AIC=rbind(AIC,AIC.t)
		
		#display only significant path [R: lavaan, semPlot]
		#ref:https://stackoverflow.com/questions/51270032/how-can-i-display-only-significant-path-lines-on-a-path-diagram-r-lavaan-sem
		pvalue_cutoff <- 0.05
		obj <- semPlot:::semPlotModel(fit)
		# save a copy of the original, so we can compare it later and be sure we removed only what we intended to remove
		# original_Pars <- obj@Pars
		check_Pars <- obj@Pars %>% filter(!(edge %in% c("int","<->") | lhs == rhs)) # this is the list of paramater to sift thru
		keep_Pars <- obj@Pars %>% filter(edge %in% c("int","<->") | lhs == rhs) # this is the list of paramater to keep asis
		test_against <- standardizedSolution(fit) %>% filter(pvalue < pvalue_cutoff, rhs != lhs)
		test_against_rev <- test_against %>% rename(rhs2 = lhs,lhs = rhs) %>% rename(rhs = rhs2)
		checked_Pars <- check_Pars %>% semi_join(test_against, by = c("lhs", "rhs")) %>% bind_rows(check_Pars %>% semi_join(test_against_rev, by = c("lhs", "rhs")))
		obj@Pars <- keep_Pars %>% bind_rows(checked_Pars)
		#let's verify by looking at the list of the edges we removed from the object
		# anti_join(original_Pars,obj@Pars)
		semPaths(obj,what = "std",layout = "groups",residuals=FALSE,edge.label.cex =5,label.cex =10,asize=15,borders=FALSE,nCharNodes=0,fade = F)#standardized parameter estimate in edge labels.
		#title(main=y.var[i])
	}	
	
	#每个地理单元内显著高于/低于平均生态位进化速率(ln)的类群
	adc=unique(dis2$Adcode99)
	q.cal2=function(i,dis2,niche.rt,type="all"){
		dis3=subset(dis2,dis2$Adcode99==i)
		if (type%in%"all") dis.t=subset(dis3,!is.na(dis3[,niche.rt]))
		if (type%in%"inc") dis.t=subset(dis3,dis3[,niche.rt]>0)
		if (type%in%"dec") dis.t=subset(dis3,dis3[,niche.rt]<0)
		#dis.t[,1]=abs(dis.t[,1])		
		if (length(grep("rate",niche.rt))==0&length(grep("K",niche.rt))==0) dis.t[,niche.rt]=dis.t[,niche.rt] else dis.t[,niche.rt]=log(abs(dis.t[,niche.rt]))
		#近似正态分布求95%区间的方法https://www.jianshu.com/p/cb53a7dc00e3,
		# -先求取左侧部分的sd,但是要补足右侧对称的数据
		# -同样求右侧部分的sd,同时补足左侧对称的数据
		# -用最高密度值时max_gc值加减1.96倍左右侧的sd
		conf_func = function(data,col){
			set = density(data[,col])
			max_midu = set$x[which.max(set$y)]
			right = data[data[,col] > max_midu,col]
			left = data[data[,col] < max_midu,col]
			max = max_midu + sd(c(right,max_midu*2-right))*1.96
			min = max_midu - sd(c(left,max_midu*2-left))*1.96			
			return(c(min,max))
		}
		conf.t=conf_func(dis.t,niche.rt)
		conf.t[is.na(conf.t)]=min(dis.t[,niche.rt])
		conf.cal=function(mu,conf.t){			
			if(mu<conf.t[1]) return(mu-conf.t[1]) else {
				if (mu>conf.t[2])return (mu-conf.t[2]) else return (0)
			}
		}		
		clade.conf=tapply(dis.t[,niche.rt],dis.t$Species_E1,conf.cal,conf.t)		
		re=data.frame(id=paste(i,names(clade.conf),sep="_"),conf=clade.conf)
		colnames(re)[2]=paste("conf",niche.rt,sep=".")				
		return(re)
	}
	#每个地理单元内，生态位进化速率最快的类群
	#每个地理单元内平均生态位宽度、pos和physig, woody.pro最高和最低的类群
	#ret=c(); for (i in 1:length(adc)) tmp=q.cal2(adc[i],dis2,"SNBT");ret=rbind(tmp,re)
	dis4=data.frame(id=paste(dis2$Adcode99,dis2$Species_E1,sep="_"),dis2[,c("Adcode99","Species_E1","clades","lifefrom","rate.bio01.mean","rate.bio12.mean","SNBT","SNBP","bio01.mean","bio12.mean")])
	rate.MAT.all=do.call(rbind,lapply(adc,q.cal2,dis2,niche.rt="rate.bio01.mean"))
	rate.MAP.all=do.call(rbind,lapply(adc,q.cal2,dis2,niche.rt="rate.bio12.mean"))
	rate.MAT.inc=do.call(rbind,lapply(adc,q.cal2,dis2,niche.rt="rate.bio01.mean",type="inc"))
	rate.MAP.inc=do.call(rbind,lapply(adc,q.cal2,dis2,niche.rt="rate.bio12.mean",type="inc"))
	rate.MAT.dec=do.call(rbind,lapply(adc,q.cal2,dis2,niche.rt="rate.bio01.mean",type="dec"))
	rate.MAP.dec=do.call(rbind,lapply(adc,q.cal2,dis2,niche.rt="rate.bio12.mean",type="dec"))
	Rate=cbind(dis4,rate.MAT.all[match(dis4$id,rate.MAT.all$id),2],rate.MAP.all[match(dis4$id,rate.MAP.all$id),2],rate.MAT.inc[match(dis4$id,rate.MAT.inc$id),2],rate.MAP.inc[match(dis4$id,rate.MAP.inc$id),2],rate.MAT.dec[match(dis4$id,rate.MAT.dec$id),2],rate.MAP.dec[match(dis4$id,rate.MAP.dec$id),2])
	vars=c("MAT","MAP","MAT_inc","MAP_inc","MAT_dec","MAP_dec")
	colnames(Rate)=c(colnames(dis4),vars)
	save(Rate, file="Rate.Rdata")
	
	#match species with order
	load("Rate.Rdata")
	Rate$Genus=do.call(rbind,lapply(strsplit(rownames(Rate),"_"),function(i){return(i[1])}))
	tpl=read.csv("Genus_synonym.csv")#from Zhiheng 2022.11.1
	Rate$Order=tpl[match(Rate$Genus,tpl$Genus_Accepted_E1),"Order_Soltis_Genbank"]
	unmat=subset(Rate,is.na(Rate$Order))
	unmat$Order=tpl[match(unmat$Genus,tpl$Genus_E1),"Order_Soltis_Genbank"]
	mat=subset(Rate,!is.na(Rate$Order))
	Rate=rbind(mat,unmat)
	save(Rate, file="Rate.Rdata")
	
	load("Rate.Rdata")
	#dis.t=Rate;var.t=vars[11];niche.rt="K.bio01.mean"
	plot.pre=function(var.t,Rate){
		dis.t=Rate[!is.na(Rate[,var.t]),]
		mf.cal=function(dis.t,var.t,max=TRUE){
			if (max==TRUE) {
				rate.dat=dis.t[dis.t[,var.t]>0,c(var.t,"Adcode99","Species_E1","lifefrom","Order")]
			}else {
				rate.dat=dis.t[dis.t[,var.t]<0,c(var.t,"Adcode99","Species_E1","lifefrom","Order")]
			}			
			adc=unique(rate.dat$Adcode99)		
			clade.cal=function(i,dat) {
				dat.cld=dat[dat$Adcode99==i,];
				clade.cout=tapply(dat.cld$Species_E1,dat.cld$Order,length)
				clade.t=names(clade.cout[which(clade.cout==max(clade.cout,na.rm=T))])
				#若并列最多/最少，取中值最大/最小的
				clade=clade.t[which.max(abs(tapply(dat.cld[dat.cld$Order%in%clade.t,var.t],dat.cld[dat.cld$Order%in%clade.t,"Order"],median,na.rm=T)))]
				return(data.frame(Adcode99=i,clade))
			}
			clade.all =do.call(rbind,lapply(adc,clade.cal,rate.dat))
			clade.all[is.na(clade.all)]=0
			return(clade.all)
		}
		dmax=mf.cal(dis.t,var.t)
		#dmin=mf.cal(dis.t,var.t,max=FALSE)
		adc2=unique(dis.t$Adcode99)
		re=data.frame(Adcode99=adc2,dmax[match(adc2,dmax$Adcode99),-1])#,dmin[match(adc2,dmin$Adcode99),-1])
		colnames(re)=c("Adcode99","clade.max")#,"clade.min")
		return(re)
	}	
	tag=list()
	for (i in 1:length(vars)) tag[[i]]=plot.pre(vars[i],Rate)
	names(tag)=vars
	
	## draw plot	
	q.plot=function(rate.dat,data.test,clade.n,cordacd){		
		datap=cbind(rate.dat,clade.n[match(rate.dat[,"clade.max"],clade.n$clade),-1])		
		cordacd2=cbind(cordacd,datap[match(cordacd$ADCODE99,datap$Adcode99),])
		pC.t4=cbind(coord,cordacd2[match(rownames(coord),rownames(cordacd2)),])		
		data.test[]=pC.t4$clade.n
		test_spdf <- as(data.test, "SpatialPixelsDataFrame")
		test_df <- as.data.frame(test_spdf)
		colnames(test_df) <- c("value", "x", "y")
		test_df2=cbind(test_df,clade.n[match(test_df$value,clade.n$clade.n),])
		col.tmp=unique(test_df2[,c("col","clade","clade.n")])
		col.tmp=col.tmp[order(col.tmp$clade.n),]
		ggplot()+geom_raster(data = test_df2 , aes(x = x, y = y,fill = factor(value))) + 
				scale_fill_manual(name = "Order",values = as.character(col.tmp$col),labels = col.tmp$clade)+
				#labs(title = paste(title,var,sep="_"))+
				theme_bw()+	
				theme(axis.title = element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
				panel.grid=element_blank(),panel.background = element_rect(fill = 'darkgray', colour = 'black'),
				legend.text=element_text(face="bold.italic",size=16),
				legend.title=element_text(face="bold",size=18)
		)			
						
	}		
	library(raster)	
	library(ggpubr)
	data.test <- raster("extend.tree/rasters/0.00_temperature.asc",native=TRUE)
	coord=coordinates(data.test)
	rownames(coord)=paste(coord[,1],coord[,2],sep="_")
	cordacd <- read.csv("POINT2.csv")
	cordacd=subset(cordacd,cordacd$ADCODE99>0)
	rownames(cordacd)=paste(cordacd[,1],cordacd[,2],sep="_")
	# clade.n=data.frame(
		# clade=c("ANA","Chloranthanae","Magnoliids","mono-commelinids","monocots",
			# "eudicots-Ceratophyllanae","eudicots","core eudicots","Dilleniales","superasterids","Caryophyllales","Santalales","asterids","asterids-campanulids",
			# "asterids-lamiids","superrosids","rosids","rosids-malvids","rosids-fabids","rosids-COM","rosids-fabids-Nfix"),
		# clade.n=1:21,
		# labels=c("ANA", "Chl", "Mgn","MonC","Mon","DicC","Dic","DicCor","Dill","SpAst","Caryo","Santa","Ast","AstC","AstL",
			# "SpRos","Ros","RosM","Fab","COM","Nfix"),
		# col=c("#B51E90", "#E13A94", "#F588B8","#FEF201","brown","#8C52A2","#B384BA","#DEC3DE","pink","#00BA38","#9EBC90","green","#2ABAA0","#03B4F0","#D6F1FE",
			# "black","darkgray","lightgray","#EC7862","#F48F74","#FBAA72"),id=rep(1,21)		
	# )
	clade.n=data.frame(
		clade=na.omit(unique(do.call(rbind,tag)$clade.max)),
		clade.n=1:16,
		col=c("#B51E90", "#F588B8","#FEF201","brown","#8C52A2","#B384BA","pink","#00BA38","green","#2ABAA0","#03B4F0","#D6F1FE",
			"black","#EC7862","#F48F74","#FBAA72"))	
	 	
	#每个地理单元内平均生态位进化速率(ln)最高和最低的类群
	plotn=vector("list",6)	
	for (i in 1:6) plotn[[i]]=q.plot(tag[[i]],data.test,clade.n,cordacd)
		
	windows(80,50)
	ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],plotn[[5]],plotn[[6]],
			nrow=3,ncol = 2,widths=c(2,2,2),heights=c(2,2),common.legend=TRUE,legend="right",
			hjust=-0.5,align="hv",labels=c("a","b","c","d","e","f"),font.label = list(size = 25))
		
##种级净分化速率与生态位进化速率的关系（净分化速率来自蒋可）
	library(data.table)	
	biodis=as.data.frame(fread("Splev.new3.biomat.csv"))[,-1]	
	rownames(biodis)=biodis$genus
	div.rate=as.data.frame(fread("tipratesmith.csv"))
	# div.rate=as.data.frame(fread("final_maga_pgls.csv"))#family level magallon
	# div.rate=as.data.frame(fread("final_bar_pgls.csv"))#family level Barbara
	vlist=paste("rate",c("bio01","bio05","bio06","bio12","bio16","bio17"),"mean",sep=".")
	comb=na.omit(cbind(div.rate,log(abs(biodis[match(div.rate$species,biodis$genus),vlist]))))
		
	xlab.t=c("MAT","MTWM","MTCM")
	xlab.p=c("MAP","MPWQ","MPDQ")	
	getscater.grid(vlist[1:3],"ndr",xlab.t,"div",comb,"MAT rate")
	getscater.grid(vlist[4:6],"ndr",xlab.p,"div",comb,"MAP rate")
	
	y.list.div=c("r0_mag","r0.9_mag","r0_bar","r0.9_bar")
	getscater.grid(c(y.list.t,y.list.p),y.list.div,c(ylab.t,ylab.p),y.list.div,dat.grid,"Divsificate rate")
	
	
###### code graveyard #######
## niche evl rate of tropic and temperte
#tropic
	biodis2=cbind(biodis,group=paste(biodis$from.tropic,biodis$to.tropic,sep="-"))	
	rate.p(biodis,"to.tropic",ratetype="sprich",dis,shape,main="TropicPro")	
	tapply(biodis2$genus,biodis2$group,length)
	# NA-tmp  NA-tro tmp-tmp tmp-tro tro-tmp tro-tro 
    # 714     736   74720    8929   12649  133819
	# sprich ~ age
	ca.sprich=function(type,biodis2,thes){
		biodis.t=biodis2[biodis2$group%in%type,c("age","genus")]		
		p=cbind(age.cat=round(biodis.t$age/thes)*thes,biodis.t)
		sprich=tapply(p$genus,p$age.cat,length)
		p1=cbind(age.cat=round(biodis2$age/thes)*thes,biodis2)
		sprich1=tapply(p1$genus,p1$age.cat,length)
		pro=sprich/sprich1[names(sprich)]
		return(data.frame(group=type,time=as.numeric(names(sprich)),sprich=sprich,pro=pro))
	}
	group=c("tmp-tmp","tmp-tro","tro-tmp","tro-tro")
	sprich=do.call(rbind,lapply(group,ca.sprich,biodis2,10))
	ggplot(data = sprich, mapping = aes(x = time, y = sprich, colour = group))+geom_line(size=1)	
	ggplot(data = sprich, mapping = aes(x = time, y = pro, colour = group)) + geom_line(size=1)	
	#tmp
	rate.p(biodis2[biodis2$to.tropic%in%"tmp",],niche.rt[1],ratetype="all",dis,shape,main="tmp")
	rate.p(biodis2[biodis2$to.tropic%in%"tmp",],niche.rt[4],ratetype="all",dis,shape,main="tmp")
	#tro
	rate.p(biodis2[biodis2$to.tropic%in%"tro",],niche.rt[1],ratetype="all",dis,shape,main="tro")
	rate.p(biodis2[biodis2$to.tropic%in%"tro",],niche.rt[4],ratetype="all",dis,shape,main="tro")
	
	#tmp-tmp
	rate.p(biodis2[biodis2$group%in%"tmp-tmp",],"genus",ratetype="sprich",dis,shape,main="tmp-tmp")
	rate.p(biodis2[biodis2$group%in%"tmp-tmp",],niche.rt[1],ratetype="all",dis,shape,main="tmp-tmp")
	rate.p(biodis2[biodis2$group%in%"tmp-tmp",],niche.rt[4],ratetype="all",dis,shape,main="tmp-tmp")	
	#tmp-tro
	rate.p(biodis2[biodis2$group%in%"tmp-tro",],"genus",ratetype="sprich",dis,shape,main="tmp-tro")
	rate.p(biodis2[biodis2$group%in%"tmp-tro",],niche.rt[1],ratetype="all",dis,shape,main="tmp-tro")
	rate.p(biodis2[biodis2$group%in%"tmp-tro",],niche.rt[4],ratetype="all",dis,shape,main="tmp-tro")
	#tro-tmp
	rate.p(biodis2[biodis2$group%in%"tro-tmp",],"genus",ratetype="sprich",dis,shape,main="tro-tmp")
	rate.p(biodis2[biodis2$group%in%"tro-tmp",],niche.rt[1],ratetype="all",dis,shape,main="tro-tmp")
	rate.p(biodis2[biodis2$group%in%"tro-tmp",],niche.rt[4],ratetype="all",dis,shape,main="tro-tmp")
	#tro-tro
	rate.p(biodis2[biodis2$group%in%"tro-tro",],"genus",ratetype="sprich",dis,shape,main="tro-tro")
	rate.p(biodis2[biodis2$group%in%"tro-tro",],niche.rt[1],ratetype="all",dis,shape,main="tro-tro")
	rate.p(biodis2[biodis2$group%in%"tro-tro",],niche.rt[4],ratetype="all",dis,shape,main="tro-tro")
	
##plot目的进化速率
	# options("install.lock"=FALSE)
	# install.packages("castor")
	fam=tapply(biodis$genus,biodis$Family_E,sample,1)
	tre=drop.tip(tre0, tip=subset(tre0$tip.label,!tre0$tip.label%in%fam))
	tre$tip.label=names(fam[match(tre$tip.label,fam)])		
	order.lst=unique(biodis$Order)
	node=c()#每个calde的节点编号
	for (i in 1:length(order.lst)){
		sp=unique(biodis[biodis$Order==order.lst[i],"Family_E"])
		node.t=castor::get_mrca_of_set(tre,sp)
		node=c(node,node.t)
	}
	names(node)=order.lst	
	niche.rt=paste("rate",c("bio01","bio05","bio06","bio12","bio16","bio17"),"mean",sep=".")	
	rate=log(tapply(abs(biodis[,niche.rt[1]]),biodis$Order,mean))
	col=data.frame(rate=rate[order(rate)],col=colorRampPalette(c("palegreen3","sandybrown"))(length(rate)))
	p=ggtree(tre,layout="circular")
	p2=collapse(p, node[1], 'max',fill =col[match(names(node)[1],rownames(col)),"col"])
	for (i in 2:length(node)) p2=p2%>% collapse(node[i], 'max',fill =col[match(names(node)[i],rownames(col)),"col"])
	
##每个地理单元生态位宽度、位置和谱系信号显著高于整体的物种，其生态位进化速率与其他类群的差异（所有显著高于整体的物种的生态位进化速率与整体进化速率的最高置信区间的差异 的平均值 ）
	load("Rate.Rdata")
	vars=paste("conf",c("rate.MAT.all","rate.MAP.all","rate.MAT.inc","rate.MAP.inc","rate.MAT.dec","rate.MAP.dec","width.MAT","width.MAP","pos.MAT","pos.MAP","physig.MAT","physig.MAP"),sep=".")
	adc=unique(Rate$Adcode99)	
	q.plot2=function(Rate,adc,cordacd,data.test,hQ=TRUE,var.t,drivers.t,var.tp=c("mean","cor")){
		if (hQ==TRUE) rate.h=Rate[Rate[,var.t]>0,] else rate.h=Rate[Rate[,var.t]<0,]
		dat=na.omit(rate.h[,c("Adcode99","Species_E1",var.t,drivers.t)])
		if (var.tp=="mean"){
			re.mean=tapply(dat[,drivers.t],dat$Adcode99,mean)
			re=data.frame(Adcode99=adc,mean=re.mean[match(adc,names(re.mean))])
		} else {
			mf.cor=function(i,dat) {dat.t=dat[dat$Adcode99==i,];return(cor(dat.t[,c(var.t,drivers.t)])[1,2])}
			re.cor=sapply(adc,mf.cor,dat)	
			re.cor[is.na(re.cor)]=0;names(re.cor)=adc
			re=data.frame(Adcode99=adc,cor=re.cor[match(adc,names(re.cor))])
		}			
		cordacd2=cbind(cordacd,re[match(cordacd$ADCODE99,re$Adcode99),])
		pC.t4=cbind(coord,cordacd2[match(rownames(coord),rownames(cordacd2)),])	
		data.test[]=pC.t4[,var.tp]
		test_spdf <- as(data.test, "SpatialPixelsDataFrame")
		test_df <- as.data.frame(test_spdf)
		colnames(test_df) <- c("value", "x", "y")	
		p=ggplot()+geom_raster(data = test_df, aes(x = x, y = y,fill = value)) +
		scale_fill_gradient2(midpoint = 0,low="red",high="blue")+
			labs(title = paste(var.t,drivers.t,var.tp,sep="_"))+
			theme_bw()+	
			theme(axis.title = element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
			panel.grid=element_blank(),panel.background = element_rect (fill = 'lightgray'),legend.position=c(0.1,0.3))
		return(p)
	}
	plotn=vector("list",6)
	for (i in c(1,3,5)) plotn[[i]]=q.plot2(Rate,adc,cordacd,data.test,hQ=TRUE,var.t=vars[i+6],drivers.t=vars[1],var.tp="mean")	
	for (i in c(2,4,6)) plotn[[i]]=q.plot2(Rate,adc,cordacd,data.test,hQ=TRUE,var.t=vars[i+6],drivers.t=vars[2],var.tp="mean")
	
	windows(80,50)
	ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],plotn[[5]],plotn[[6]],
			nrow=3,ncol = 2,widths=c(2,2,2),heights=c(2,2),common.legend=TRUE,legend = "right",hjust=0,align="hv")
		
## 整合科水平的信息 ---
	##年龄、分布范围、科内物种数、科内属的个数、净分化速率、草木本比例、热带类群比例、生态位大小对进化速率的散点图	##年龄、分布范围、科内物种数、科内属的个数、净分化速率、草木本比例、热带类群比例、进化速率、生态位大小、生态位的谱系信号对进化速率谱系信号的散点图
	biodis=as.data.frame(fread("Splev.new3.biomat.csv"))[,-1]	
	rownames(biodis)=biodis$genus	
	mf=function(x) length(x[x%in%"W"])/(length(x[x%in%"W"])+length(x[x%in%"H"]))
	woody.pro=tapply(biodis$lifefrom,biodis$Family_E,mf)
	
	mf.tro=function(x) length(x[x%in%"tro"])/(length(x[x%in%"tro"])+length(x[x%in%"tmp"]))
	tro.pro=tapply(biodis$to.tropic,biodis$Family_E,mf.tro)
	
	#family level physig
	physig=na.omit(read.csv("physig.all.csv")[,-1])
	rownames(physig)=physig$Family	
	vlist0=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")
	vlist=paste("rate",vlist0,sep=".")	
	vlst=c(vlist0,vlist)
	psig.lambda=do.call(cbind,lapply(1:length(vlst),function(i){
		corg.t=physig[,c(paste("Lambda",vlst[i],sep="_"),paste("Lambda",".p_",vlst[i],sep=""))]
		#corg.t[corg.t[,2]>0.05|corg.t[,1]>=1.2,1]=0
		corg=corg.t[,1]
		return(corg)
	}))
	colnames(psig.lambda)=paste("lmd",vlst,sep=".")	
	
	psig.K=do.call(cbind,lapply(1:length(vlst),function(i){
		corg.t=physig[,c(paste("K",vlst[i],sep="_"),paste("K",".p_",vlst[i],sep=""))]
		#corg.t[corg.t[,2]>0.05|corg.t[,1]>=1.2,1]=0
		corg=corg.t[,1]
		return(corg)
	}))
	colnames(psig.K)=paste("K",vlst,sep=".")	
	psigall=cbind(psig.lambda,psig.K)
	
	#family level nich & nich evl rate
	nich.fam=do.call(cbind,lapply(1:length(vlst),function(i){
		tapply(abs(biodis[,vlst[i]]),biodis$Family_E,median)
	}))
	colnames(nich.fam)=vlst
	
	#family age	
	famlylst=as.character(unique(biodis$Family_E))
	sp.age=unlist(lapply(1:length(famlylst),function(m){	
		sp=biodis[biodis$Family_E%in%famlylst[m],"genus"]
		tre2=drop.tip(tre, tip=subset(tre$tip.label,!tre$tip.label%in%sp))
		if (Ntip(tre2)<=2) {
			age=mean(biodis[biodis$Family_E%in%famlylst[m],"age"])
		} else {
			## calculate root length(crown) of the tree
			root.dist <- function(tree, method = c("n.node","branch.length"), type=c("tip","node","both"))
		{
		if (length(method) > 1) method <- method[1]
		if (length(type) > 1) type <- type[1]
		
		if (class(tree) != "phylo") stop("The tree is not a phylo tree")
		
		## find the root
		root.label <- unique(tree$edge[,1][!tree$edge[,1] %in% tree$edge[,2]])
		
		tree.edge <- rbind(tree$edge,c(0,root.label))
		ii <- order(tree.edge[,2])
		tree.edge <- tree.edge[ii,]
		
		## the results of root distance
		N.tip <- Ntip(tree)
		N.node <- tree$Nnode
		N.edge <- Nedge(tree)
		if (type == "tip") {
			N.ii <- 1:N.tip
			rd <- numeric(N.tip)
			names(rd) <- tree$tip.label
			}
		else if (type =="node") {
			N.ii <- N.tip + (1:N.node)
			rd <- numeric(N.node)
			names(rd) <- tree$node.label
			}
		else if (type == "both") {
			N.ii <- 1:(N.tip + N.node)
			rd <- numeric(N.tip + N.node)
			names(rd) <- c(tree$tip.label, tree$node.label)
			}
		
		## count the number of nodes from root to tips
		if (method == "n.node")
			{
			for (i in 1:length(N.ii))
				{
				tip <- descendant <- N.ii[i]
				for (j in 1:N.node)
					{
					if (tree.edge[descendant,1] == 0)
						{
						rd[i] <- 0
						break
						}
					if (tree.edge[descendant,1] == root.label) 
						{
						rd[i] <- j
						break
						}
					ancestor <- tree.edge[descendant,1]
					descendant <- tree.edge[ancestor,2]
					}
				}
			}
		else if (method == "branch.length")
			{
			if (is.null(tree$edge.length)) edge.length <- rep(1, times = N.tip + N.node)
			if (!is.null(tree$edge.length))	edge.length <- c(tree$edge.length, 0)[ii]
			
			for (i in 1:length(N.ii))
				{
				tip <- descendant <- N.ii[i]
				for (j in 1:N.node)
					{
					rd[i] <- rd[i] + edge.length[descendant]
					if (tree.edge[descendant,1] == root.label | tree.edge[descendant,1] == 0) 
						{
						#rd[i] <- rd[i] + edge.length[root.label]
						break
						}
					ancestor <- tree.edge[descendant,1]
					descendant <- tree.edge[ancestor,2]
					}
				}		
			}
		return(rd)
		}
		age=max(root.dist(tre2, method = "branch.length", type="both"))		
		}		
		return(age)
	}))
	names(sp.age)=famlylst
	
	#family rangesize
	dis.fam=unique(as.data.frame(fread("SpLevDis2.csv"))[,c("Adcode99","Species_E1","Shape_Area")])	
	rg.adc=tapply(dis.fam$Adcode99,dis.fam$Family_E,length)
	rg.area=tapply(dis.fam$Shape_Area,dis.fam$Family_E,sum)
	
	#family number of sp and genus 
	sp.num=tapply(biodis$genus,biodis$Family_E,length)
	biodis2=unique(biodis[,c("Family_E","Genus_E")])
	genus.num=tapply(biodis2$Genus_E,biodis2$Family_E,length)
		
	#净分化速率 蒋可
	rate.bar=read.csv("div_rate_fam_JangKe/Barbara diver rate.csv")
	rate.mag=read.csv("div_rate_fam_JangKe/magallon diver rate.csv")
	rate=cbind(rate.mag[,c("family","r0","r0.9")],rate.bar[match(rate.mag$family,rate.bar$family),c("r0","r0.9")])
	colnames(rate)=c("family","r0_mag","r0.9_mag","r0_bar","r0.9_bar")
	
	biostat.fam=data.frame(n.sp=sp.num,n.genus=genus.num,
		rg.adc=rg.adc[match(rownames(sp.num),names(rg.adc))],rg.area=rg.area[match(rownames(sp.num),names(rg.adc))],
		age=sp.age[match(rownames(sp.num),names(sp.age))],
		rate[match(rownames(sp.num),rate$family),-1],
		woody.pro=woody.pro[match(rownames(sp.num),names(woody.pro))],
		nich.fam,psigall[match(rownames(sp.num),rownames(psigall)),],
		tro.pro=tro.pro[match(rownames(sp.num),names(tro.pro))])		
	save(biostat.fam,file="biostat.fam.Rdata")
	 
	load("biostat.fam.Rdata")
	## compare fam age from different trees ---
	##magallon 2020 nee
	mag=read.csv("maglon2020_Ages.csv")
	fam.mag=tapply(mag$Crown_BEAST,mag$Family,max)
	##Li Dezhu 2019 Nature Plants	
	library(ape)
	li.tre=read.nexus("2881_dating_Species.tre")
	famlylst=as.character(unique(read.csv("lidezhu.csv")[,"Family"]))
	fam.li=unlist(lapply(1:length(famlylst),function(m){	
		tre2=drop.tip(tre, tip=subset(tre$tip.label,!stringr::str_detect(tre$tip.label,famlylst[m])))
		if (is.null(tre2)) {
			age=NA
		} else {
			## calculate root length(crown) of the tree
			root.dist <- function(tree, method = c("n.node","branch.length"), type=c("tip","node","both"))
		{
		if (length(method) > 1) method <- method[1]
		if (length(type) > 1) type <- type[1]
		
		if (class(tree) != "phylo") stop("The tree is not a phylo tree")
		
		## find the root
		root.label <- unique(tree$edge[,1][!tree$edge[,1] %in% tree$edge[,2]])
		
		tree.edge <- rbind(tree$edge,c(0,root.label))
		ii <- order(tree.edge[,2])
		tree.edge <- tree.edge[ii,]
		
		## the results of root distance
		N.tip <- Ntip(tree)
		N.node <- tree$Nnode
		N.edge <- Nedge(tree)
		if (type == "tip") {
			N.ii <- 1:N.tip
			rd <- numeric(N.tip)
			names(rd) <- tree$tip.label
			}
		else if (type =="node") {
			N.ii <- N.tip + (1:N.node)
			rd <- numeric(N.node)
			names(rd) <- tree$node.label
			}
		else if (type == "both") {
			N.ii <- 1:(N.tip + N.node)
			rd <- numeric(N.tip + N.node)
			names(rd) <- c(tree$tip.label, tree$node.label)
			}
		
		## count the number of nodes from root to tips
		if (method == "n.node")
			{
			for (i in 1:length(N.ii))
				{
				tip <- descendant <- N.ii[i]
				for (j in 1:N.node)
					{
					if (tree.edge[descendant,1] == 0)
						{
						rd[i] <- 0
						break
						}
					if (tree.edge[descendant,1] == root.label) 
						{
						rd[i] <- j
						break
						}
					ancestor <- tree.edge[descendant,1]
					descendant <- tree.edge[ancestor,2]
					}
				}
			}
		else if (method == "branch.length")
			{
			if (is.null(tree$edge.length)) edge.length <- rep(1, times = N.tip + N.node)
			if (!is.null(tree$edge.length))	edge.length <- c(tree$edge.length, 0)[ii]
			
			for (i in 1:length(N.ii))
				{
				tip <- descendant <- N.ii[i]
				for (j in 1:N.node)
					{
					rd[i] <- rd[i] + edge.length[descendant]
					if (tree.edge[descendant,1] == root.label | tree.edge[descendant,1] == 0) 
						{
						#rd[i] <- rd[i] + edge.length[root.label]
						break
						}
					ancestor <- tree.edge[descendant,1]
					descendant <- tree.edge[ancestor,2]
					}
				}		
			}
		return(rd)
		}
		age=max(root.dist(tre2, method = "branch.length", type="both"))		
		}		
		return(age)
	}))
	names(fam.li)=famlylst
	biostat.fam2=cbind(age.mag=fam.mag[rownames(biostat.fam)],age.li=fam.li[rownames(biostat.fam)],biostat.fam)
	
	par(mfrow = c(1,2),mar=c(0.5,4,0.5,1),oma=c(2,2,2,2))
	plot(na.omit(biostat.fam2[,c("age","age.mag")]))
	abline(coef = c(0,1))	
	plot(na.omit(biostat.fam2[,c("age","age.li")]))
	abline(coef = c(0,1))	
	cor(na.omit(biostat.fam2[,c("age","age.mag","age.li")]))
	age.sum=biostat.fam2[,c("age","age.mag","age.li")]
	write.csv(age.sum,"age.sum.csv")
		
	load("biostat.fam.Rdata")	
##plot fam level data	
	#1. scanter plot	
	getscater=function(xlog.list,xnolog.list,y.list,biostat.fam,log.x=TRUE,log.y=TRUE){
		ss.glm <- function(r.glm)
					{
					r.ss <- summary(r.glm)
					rsq <- 100*(r.ss$null.deviance-r.ss$deviance)/r.ss$null.deviance
					adj.rsq <- 100*(1-(r.ss$deviance/r.ss$df.residual)/(r.ss$null.deviance/r.ss$df.null))
					f.stat <- ((r.ss$null.deviance-r.ss$deviance)/(r.ss$df.null-
							r.ss$df.residual))/(r.ss$deviance/r.ss$df.residual)
					p <- pf(f.stat, r.ss$df.null-r.ss$df.residual, r.ss$df.residual, lower.tail=FALSE)

					return(c(r2=rsq,adj.r2=adj.rsq,p=p))
					}
		if (log.x==TRUE) {
			temp.x<-cbind(log(biostat.fam[,xlog.list]),biostat.fam[,xnolog.list])
			xlab <- c(paste("ln",xlog.list,sep=" "),xnolog.list)
		} else{
			temp.x<-biostat.fam[,c(xlog.list,xnolog.list)]
			xlab <- c(xlog.list,xnolog.list)
		}
		if (log.y==TRUE) {
			temp.y<-log(biostat.fam[,y.list])
			ylab <- paste("ln",c("bio01","bio05","bio06","bio12","bio16","bio17"),sep=" ")
		} else{
			temp.y<-biostat.fam[,y.list]
			ylab <- c("bio01","bio05","bio06","bio12","bio16","bio17")
		}
		temp=cbind(temp.x,temp.y)		
		x.list=c(xlog.list,xnolog.list)
		width <- 70; height <- 40
		windows(width=width, height=height)
		n.col <- length(x.list); n.row <- length(y.list)
		mylayout <- layout(matrix(1:((2+n.col)*(2+n.row)), ncol=2+n.col, byrow=T), 
				width=c(0.7*width/(2+n.col), rep(width/(2+n.col),times = n.col), 0.3*width/(2+n.col)),
				height=c(0.3*height/(2+n.row), rep(height/(2+n.row),times=n.row), 0.7*height/(2+n.row)))
		layout.show((2+n.col)*(2+n.row))
		
		glmr=c()
		for (i in 1:length(y.list)) {		
			glmr.t2=c()
			if (i == 1) {
				for (j in 1:(2+n.col)) {par(mar=c(0,0,0,0)); plot.new()}
				}	
			for (xi in 1:length(x.list)) {
				temp2=na.omit(temp[,c(x.list[xi],y.list[i])])				
				x <- temp2[,1]
				y <- temp2[,2]			
				if (xi == 1) {par(mar=c(0,0,0,0)); plot.new()}
			  par(mar=c(0.5,0.5,0,0), cex=1, cex.axis=0.8, cex.lab=1.2, mgp=c(2.4,0.3,0), tck=-0.04)
			  cex.axis1=1.4;cex1=1.3;tck1=-0.02;cex.lab1=1.8;cex.legend = 1.4;cex.txt=1.4
			   plot(y~x, data=temp, pch=19, xlab="", ylab="", cex=0.25, col='grey', cex.axis=cex.axis1, axes=F); box()
				 m <- glm(y~x,data=temp)
				 glmr.t=data.frame(x=x.list[xi],y=y.list[i],adj.r2=ss.glm(m)[2],p=ss.glm(m)[3])
				 a <- seq(min(x), max(x), 0.1)
			  if (i!=4) if (ss.glm(m)[3] < 0.05)  lines(x=a, y=predict(m, newdata=list(x=a)), col="black", lty=1, lwd=0.8)
			 if (xi == 1) {
					axis(side=2, label=T)
					mtext(side=2, text=ylab[i], line=1.6, cex=1, las=0)
					}
				if (i == length(y.list)) {
					axis(side=1, lwd = 0.5)
					mtext(side=1, text=xlab[xi], line=1.6, cex=1, las=0)
					}
			if (xi == 1 & i == 0.5+0.5*length(y.list)) {
						 mtext(side=2, text="proportion%", line=5.6, cex=1.2)
					 }
				if (xi == length(x.list)) {par(mar=c(0,0,0,0)); plot.new()}
			glmr.t2=rbind(glmr.t2,glmr.t)	
			}
		glmr=rbind(glmr,glmr.t2)		
		}
		return(glmr)
	}
	
	getscater2=function(x.list,y.list,biostat.fam,log.x=FALSE,log.y=TRUE){
		require(ggpubr)
		ss.glm <- function(r.glm)
					{
					r.ss <- summary(r.glm)
					rsq <- 100*(r.ss$null.deviance-r.ss$deviance)/r.ss$null.deviance
					adj.rsq <- 100*(1-(r.ss$deviance/r.ss$df.residual)/(r.ss$null.deviance/r.ss$df.null))
					f.stat <- ((r.ss$null.deviance-r.ss$deviance)/(r.ss$df.null-
							r.ss$df.residual))/(r.ss$deviance/r.ss$df.residual)
					p <- pf(f.stat, r.ss$df.null-r.ss$df.residual, r.ss$df.residual, lower.tail=FALSE)

					return(c(r2=rsq,adj.r2=adj.rsq,p=p))
					}
		if (log.x==TRUE) {temp.x<-log(biostat.fam[,x.list])} else{
			temp.x<-biostat.fam[,x.list]}
		if (log.y==TRUE) {temp.y<-log(biostat.fam[,y.list])} else{
			temp.y<-biostat.fam[,y.list]}
		temp=cbind(temp.x,temp.y)
		leg=c("a","b","c","d","e","f")
		plotn=vector("list",6)
		glmr=c()
		for (i in 1:6){
		plt=na.omit(temp[,c(x.list[i],y.list[i])])
		colnames(plt)=c("xvar","yvar")
		p=ggplot(data = plt,mapping=aes(x=xvar,y=yvar))+geom_point(size=2,colour="black",fill="#00BA38",shape=21)+		
		theme(panel.grid.major =element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black"),
			#axis.title.x = element_blank(),
			#axis.title.y = element_blank(),
			axis.text.x  = element_text(size=12),
			axis.text.y  = element_text(size=12))+			
			annotate("text",x=min(plt$xvar),y=max(plt$yvar),label=leg[i],size=8)
		if (log.x==TRUE) {p=p+xlab(paste(x.list[i],"(ln)",sep=""))} else{
			p=p+xlab(x.list[i])}
		if (log.y==TRUE) {p=p+ylab(paste(y.list[i],"(ln)",sep=""))} else{
			p=p+ylab(y.list[i])}				
		d=glm(yvar~xvar,data = plt)	
		if(ss.glm(d)[3]<0.05) p=p+geom_smooth(method = "glm",colour="red",fill="#00BD5F",size=1)
		if(ss.glm(d)[3]<0.1&ss.glm(d)[3]>0.05) p=p+geom_smooth(method = "glm",colour="red",,fill="#00BD5F",linetype="dashed")#marginal sig		
		plotn[[i]]=p
		glmr.t=data.frame(x=x.list[i],y=y.list[i],adj.r2=ss.glm(d)[2],p=ss.glm(d)[3])
		glmr=rbind(glmr,glmr.t)
		}
		windows(width=60, height=40)
		figure=ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],plotn[[5]],plotn[[6]],
		nrow=2,ncol = 3,widths=c(2,2,2),heights=c(2,2))	
		return(list(glmr,figure))
	}
		
	#去掉1Ma以下的科
	biostat.fam=biostat.fam[biostat.fam$age>=1,]
	niche=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")
	niche.lmd=paste("lmd",niche,sep=".")
	niche.K=paste("K",niche,sep=".")
	
	niche.rt=paste("rate",niche,sep=".")
	niche.rt.lmd=paste("lmd",niche.rt,sep=".")
	niche.rt.K=paste("K",niche.rt,sep=".")
	
	vlist.x.nlog=c("age","r0_mag","r0.9_mag","r0_bar","r0.9_bar","woody.pro","tro.pro")
	vlist.x.log=c("n.sp","n.genus","rg.adc","rg.area")
	
	rt.xlst=getscater(vlist.x.log,vlist.x.nlog,niche.rt,biostat.fam,log.x=TRUE,log.y=TRUE)
	K.xlst=getscater(vlist.x.log,vlist.x.nlog,niche.K,biostat.fam,log.x=TRUE,log.y=FALSE)	
	getscater2(niche.K,niche.rt,biostat.fam,log.x=FALSE,log.y=TRUE)	
	
	##2.SEM
	load("biostat.fam.Rdata")
	biostat.fam=biostat.fam[biostat.fam$age>=1,]
	library(lavaan)
	library(semPlot)
	library(magrittr)
	library(dplyr)
	niche=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")
	niche.lmd=paste("lmd",niche,sep=".")
	niche.K=paste("K",niche,sep=".")
	
	niche.rt=paste("rate",niche,sep=".")
	niche.rt.lmd=paste("lmd",niche.rt,sep=".")
	niche.rt.K=paste("K",niche.rt,sep=".")
	
	par(mfrow = c(2,3),mar=c(0,0,0.4,0),oma=c(0,0,0,0))
	AIC=c()	
	for(i in 1:length(niche.rt)){
		sem.data=na.omit(biostat.fam[,c(niche.rt[i],niche.K[i],"rg.area","woody.pro","tro.pro","age","r0.9_bar")])
		sem.data[,c(niche.rt[i],"rg.area")]=log(sem.data[,c(niche.rt[i],"rg.area")])
		colnames(sem.data)=c("NichRate","PhyloSig","Range","WoodyPro","TroPro","FamAge","DivRate")
		myModel1="DivRate~Range+WoodyPro+TroPro+FamAge+PhyloSig+NichRate
			NichRate~Range+WoodyPro+TroPro+FamAge+PhyloSig
			PhyloSig~Range+WoodyPro+TroPro+FamAge
			Range~~WoodyPro+TroPro+FamAge"			
		fit1 <- sem(myModel1, data = scale(sem.data))	

		myModel2="DivRate~Range+WoodyPro+TroPro+FamAge+PhyloSig
			NichRate~Range+WoodyPro+TroPro+FamAge+PhyloSig+DivRate
			PhyloSig~Range+WoodyPro+TroPro+FamAge
			Range~~WoodyPro+TroPro+FamAge"			
		fit2 <- sem(myModel2, data = scale(sem.data))
		#AIC值越小的模型表明越有可能准确地预测新数据
		if (AIC(fit1)<AIC(fit2)) {fit=fit1;best="niche2div"} else {fit=fit2;best="div2niche"}
		AIC.t=data.frame(SEM1=AIC(fit1),SEM2=AIC(fit2),Best=best)
		AIC=rbind(AIC,AIC.t)
		#display only significant path [R: lavaan, semPlot]
		#ref:https://stackoverflow.com/questions/51270032/how-can-i-display-only-significant-path-lines-on-a-path-diagram-r-lavaan-sem
		pvalue_cutoff <- 0.05
		obj <- semPlot:::semPlotModel(fit)
		# save a copy of the original, so we can compare it later and be sure we removed only what we intended to remove
		# original_Pars <- obj@Pars
		check_Pars <- obj@Pars %>% filter(!(edge %in% c("int","<->") | lhs == rhs)) # this is the list of paramater to sift thru
		keep_Pars <- obj@Pars %>% filter(edge %in% c("int","<->") | lhs == rhs) # this is the list of paramater to keep asis
		test_against <- standardizedSolution(fit) %>% filter(pvalue < pvalue_cutoff, rhs != lhs)
		test_against_rev <- test_against %>% rename(rhs2 = lhs,lhs = rhs) %>% rename(rhs = rhs2)
		checked_Pars <- check_Pars %>% semi_join(test_against, by = c("lhs", "rhs")) %>% bind_rows(check_Pars %>% semi_join(test_against_rev, by = c("lhs", "rhs")))
		obj@Pars <- keep_Pars %>% bind_rows(checked_Pars)
		#let's verify by looking at the list of the edges we removed from the object
		# anti_join(original_Pars,obj@Pars)
		semPaths(obj,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex =2,label.cex =3.5,borders=FALSE,nCharNodes=0,fade = F)		
		title(main=niche[i])
	}	
	
	##3. map phylosig in the map
	library(data.table)
	dis.fam=unique(as.data.frame(fread("SpLevDis2.csv"))[,c("Adcode99","Family_E")])	
	load("biostat.fam.Rdata")
	biostat.fam=biostat.fam[biostat.fam$age>=1,]	
	niche=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")
	niche.lmd=paste("lmd",niche,sep=".")
	niche.K=paste("K",niche,sep=".")
	
	library(sp)
	library(maptools)
	shape<-readShapeSpatial("county_islandrm.shp")			
	shape2=readShapeSpatial("coast_world/coast_world.shp")
	library(foreign)
	a=read.dbf("coast_world/Output.dbf")
	shape3=subset(shape2,shape2@data$File_ID%in%a$File_ID)	
	getColor = function(mapdata, provname, provcol, othercol){
			f = function(x, y) ifelse(x %in% y, which(y == x), 0);
			colIndex = sapply(mapdata@data$ADCODE99, f, provname);
			fg = c(othercol, provcol)[colIndex + 1];
			return(fg);
		}
	
	phylosig.p=function(biostat.fam,niche.sig,dis.fam,shape,weighted=TRUE,log_trans=FALSE){
		bios.fam=na.omit(biostat.fam[,c("n.sp",niche.sig)])
		if (log_trans==TRUE) bios.fam[,niche.sig]=log(bios.fam[,niche.sig])
		phylosig.dis=na.omit(cbind(dis.fam,bios.fam[match(dis.fam$Family_E,rownames(bios.fam)),]))		
		adc=as.character(unique(phylosig.dis$Adcode99))
		if (weighted==TRUE){
			phylosig=unlist(lapply(1:length(adc),function(j){
				tmp=phylosig.dis[phylosig.dis$Adcode99==adc[j],]
				sig=sum(tmp[,4]*tmp[,3])/sum(tmp[,3])
				return(sig)
			}))	
			names(phylosig)=adc
			#cap=paste(substr(niche.sig, 1,nchar(niche.sig)-5),"(wt)",sep="")
			cap=paste(niche.sig,"(wt)",sep="")			
		}else{
			phylosig=tapply(phylosig.dis[,4],phylosig.dis[,1],mean)
			#cap=substr(niche.sig, 1,nchar(niche.sig)-5)
			cap=niche.sig
		}			
		shp.fam=shape
		shp.fam@data=cbind(shp.fam@data,phylosig=phylosig[match(shp.fam@data$ADCODE99,names(phylosig))])
		
		pop0=shp.fam@data[,c("ADCODE99","phylosig")]				
		pop=pop0[pop0[,2]!=0,2];provname=as.character(pop0[pop0[,2]!=0,1])
		col=data.frame(cwe=pop[order(pop)],col=colorRampPalette(c("palegreen3","sandybrown"))(length(pop)))
		provcol=as.character(col[match(pop,col$cwe),2])					
		
		width <- 60; height <- 30
		windows(width=width, height=height)	
		par(fig=c(0,1,0,1),new=FALSE)
		plot(shp.fam, col = getColor(shp.fam, provname, provcol, "gray"),border = getColor(shp.fam, provname, provcol, "gray"))
		plot(shape3,border="darkgray",add=T)
		#text(-170,80,"a")
		box()
		#legend
		par(fig=c(0.075,0.2,0.05,0.55),new=TRUE)				
		barplot(as.matrix(rep(1,length(pop))),col=as.character(col$col),horiz=F,axes=F,border = NA,main =cap)
		axis(2,c(seq(1,length(pop),round(length(pop)/4)),length(pop)),signif(sort(pop)[c(seq(1,length(pop),round(length(pop)/4)),length(pop))],2))
		box()	
	}	
	
	phylosig.p(biostat.fam,niche.K[1],dis.fam,shape,weighted=FALSE)
	phylosig.p(biostat.fam,niche.K[2],dis.fam,shape,weighted=FALSE)
	phylosig.p(biostat.fam,niche.K[3],dis.fam,shape,weighted=FALSE)
	phylosig.p(biostat.fam,niche.K[4],dis.fam,shape,weighted=FALSE)
	phylosig.p(biostat.fam,niche.K[5],dis.fam,shape,weighted=FALSE)
	phylosig.p(biostat.fam,niche.K[6],dis.fam,shape,weighted=FALSE)
	
	phylosig.p(biostat.fam,niche.K[1],dis.fam,shape,weighted=TRUE)
	phylosig.p(biostat.fam,niche.K[2],dis.fam,shape,weighted=TRUE)
	phylosig.p(biostat.fam,niche.K[3],dis.fam,shape,weighted=TRUE)
	phylosig.p(biostat.fam,niche.K[4],dis.fam,shape,weighted=TRUE)
	phylosig.p(biostat.fam,niche.K[5],dis.fam,shape,weighted=TRUE)
	phylosig.p(biostat.fam,niche.K[6],dis.fam,shape,weighted=TRUE)
	
	#
	phylosig.p(biostat.fam,"rate.bio12.mean",dis.fam,shape,weighted=FALSE,log_trans=TRUE)
	
	phylosig.p(biostat.fam,"r0_mag",dis.fam,shape,weighted=TRUE)
	phylosig.p(biostat.fam,"r0.9_mag",dis.fam,shape,weighted=TRUE)
	phylosig.p(biostat.fam,"r0_bar",dis.fam,shape,weighted=TRUE)
	phylosig.p(biostat.fam,"r0.9_bar",dis.fam,shape,weighted=TRUE)
	
## 每个地理单元内生态位进化速率前25%和后25%中木本植物的比例以及主要是哪个clade贡献
	adc=unique(dis2$Adcode99)
	q.cal=function(i,dis2,niche.rt,type="all",n=0){
		dis3=subset(dis2[,c(niche.rt,"Species_E1","clades","lifefrom")],dis2$Adcode99==i)
		mf=function(dis.q,i,n,type){
			if (dim(dis.q)[1]==0) {
				re.t= data.frame(Adcode99=i,niche.quantile=n,type=type,clade=NA,woody.pro=NA)
			}else{
				clade.q=tapply(dis.q$Species_E1,dis.q$clades,length)
				woody.pro=dim(subset(dis.q,dis.q$lifefrom=="W"))[1]/dim(na.omit(dis.q))[1]*100
				re.t=data.frame(Adcode99=i,niche.quantile=n,type=type,clade=names(clade.q[which.max(clade.q)]),woody.pro=woody.pro)
			}			
			return(re.t)
		}		
		if (length(niche.rt)>1) {
			dis.t=subset(dis3,!(is.na(dis3[,niche.rt[1]]))&!is.na(dis3[,niche.rt[2]]))			
			dis.q1=subset(dis.t,dis.t[,1]!=dis.t[,2])
			dis.q2=subset(dis.q1,dis.q1[,"from.tropic"]%in%"tro")#热带起源的温带物种
			dis.q3=subset(dis.q1,dis.q1[,"from.tropic"]%in%"tmp")#温带起源的热带物种
			re=rbind(data.frame(niche="trans",mf(dis.q1,i,n,type)),
				data.frame(niche="tro.ori",mf(dis.q2,i,n,type)),data.frame(niche="tmp.ori",mf(dis.q3,i,n,type)))		
		} else {
			if (type%in%"all") dis.t=subset(dis3,!is.na(dis3[,niche.rt]))
			if (type%in%"inc") dis.t=subset(dis3,dis3[,niche.rt]>0)
			if (type%in%"dec") dis.t=subset(dis3,dis3[,niche.rt]<0)			
			if (n==0) dis.q=subset(dis.t,dis.t[,1]%in%c("W","tro")) else {
				dis.t[,1]=abs(dis.t[,1])
				if (n<0.5) dis.q=subset(dis.t,dis.t[,1]<=quantile(dis.t[,1],n)) 
				if (n>0.5) dis.q=subset(dis.t,dis.t[,1]>=quantile(dis.t[,1],n))  
				if (n==0.5) dis.q=subset(dis.t,dis.t[,1]>=quantile(dis.t[,1],0.475)&dis.t[,1]<=quantile(dis.t[,1],0.625))
			}			
			re=data.frame(niche=niche.rt,mf(dis.q,i,n,type))
		}		
		return(re)
	}
	
	lowQ=0.01;highQ=0.99 #n值越大，受clade内部物种数影响越大，是否想办法消除物种数的影响？比如重采样？但物种数越大，确实对该地区进化速率的权重越大
	#每个地理单元内，生态位进化速率最快/慢/中值的物种木本植物的占比以及主要来自哪个类群
	MAT=list(
		MAT.lowQ.all=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio01.mean",type="all",n=lowQ)),
		MAT.highQ.all=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio01.mean",type="all",n=highQ)),
		MAT.lowQ.inc=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio01.mean",type="inc",n=lowQ)),
		MAT.highQ.inc=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio01.mean",type="inc",n=highQ)),
		MAT.lowQ.dec=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio01.mean",type="dec",n=lowQ)),
		MAT.highQ.dec=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio01.mean",type="dec",n=highQ))		
	)	
	MAP=list(
		MAP.lowQ.all=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio12.mean",type="all",n=lowQ)),
		MAP.highQ.all=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio12.mean",type="all",n=highQ)),
		MAP.lowQ.inc=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio12.mean",type="inc",n=lowQ)),
		MAP.highQ.inc=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio12.mean",type="inc",n=highQ)),
		MAP.lowQ.dec=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio12.mean",type="dec",n=lowQ)),
		MAP.highQ.dec=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio12.mean",type="dec",n=highQ))		
	)	
	rate.mid=list(#每个单元格内，降水生态位进化速率在中值（0.475-0.625）的物种中，占比最大的类群
		MAT.midQ.all=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio01.mean",type="all",n=0.5)),
		MAT.midQ.inc=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio01.mean",type="inc",n=0.5)),
		MAT.midQ.dec=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio01.mean",type="dec",n=0.5)),
		MAT.midQ.all=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio01.mean",type="all",n=0.5)),
		MAT.midQ.inc=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio01.mean",type="inc",n=0.5)),
		MAT.midQ.dec=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="rate.bio01.mean",type="dec",n=0.5))
	)		
	niche.wid=list(#每个地理单元内，生态位最宽/窄/中值的物种木本植物的占比以及主要来自哪个类群
		SNBT.lowQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="SNBT",n=lowQ)),
		SNBP.lowQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="SNBP",n=lowQ)),
		SNBT.highQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="SNBT",n=highQ)),
		SNBP.highQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="SNBP",n=highQ)),
		SNBT.midQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="SNBT",n=0.5)),		
		SNBP.midQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="SNBP",n=0.5))		
	)	
	niche.pos=list(#每个地理单元内，生态位最高/低/中值的物种木本植物的占比以及主要来自哪个类群
		MAT.niche.lowQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="bio01.mean",n=lowQ)),
		MAP.niche.lowQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="bio12.mean",n=lowQ)),
		MAT.niche.highQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="bio01.mean",n=highQ)),
		MAP.niche.highQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="bio12.mean",n=highQ)),
		MAT.niche.midQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="bio01.mean",n=0.5)),		
		MAP.niche.midQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="bio12.mean",n=0.5))		
	)	
	physig=list(#每个地理单元内，生态位最保守/不保守/中值的物种木本植物的占比以及主要来自哪个类群
		MAT.physig.lowQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="K.bio01.mean",n=lowQ)),
		MAP.physig.lowQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="K.bio12.mean",n=lowQ)),
		MAT.physig.highQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="K.bio01.mean",n=highQ)),
		MAP.physig.highQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="K.bio12.mean",n=highQ)),
		MAT.physig.midQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="K.bio01.mean",n=0.5)),		
		MAP.physig.midQ=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="K.bio12.mean",n=0.5))		
	)
	#每个地理单元内，祖先与后代分布区不同的物种木本植物的占比以及主要来自哪个类群
	tro.trans=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt=c("from.tropic","to.tropic")))
	#每个地理单元内，木本植物主要来自哪个类群
	woody=do.call(rbind,lapply(adc,q.cal,dis2,niche.rt="lifefrom"))
	trans=list(
		dis.shift=subset(tro.trans,tro.trans$niche%in%"trans"),
		woody=woody,	
		tro2tmp=subset(tro.trans,tro.trans$niche%in%"tro.ori"),
		tmp2tro=subset(tro.trans,tro.trans$niche%in%"tmp.ori")			
	)	
	
	drivers=list(MAT,MAP,rate.mid,niche.wid,niche.pos,physig,trans)
	save(drivers, file="drivers0.05.Rdata")
	save(drivers, file="drivers0.01.Rdata")#lowQ=0.01;highQ=0.99
	#每个地理单元内，生态位进化速率最快的类群	
	tag=drivers[[7]]#trans	
	plotn=vector("list",4)	
	for (i in 1:4) plotn[[i]]=q.plot(tag[[i]],data.test,clade.n,cordacd,names(tag[i]),"clade")		
	windows(80,35)	
	ggarrange(ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			nrow=2,ncol = 2,widths=c(2,2,2),heights=c(2,2),common.legend=FALSE,hjust=0,align="hv"),leg,
		nrow=1,ncol = 2,widths=c(9,0.8),heights=2)
	
## 计算迁移过程中的地理隔离程度 ----	
	past.clim0=read.csv("past.clim.csv") #每个geounit对应的生物地理区
	past.clim0[,grepl("prep",colnames(past.clim0))]=365*past.clim0[,grepl("prep",colnames(past.clim0))]
	clim=read.csv("vars.csv")[,c("ADCODE99","Lon","Lat","bio01","bio05","bio06","bio12","bio16","bio17")]#Worldclim from luoao
	past.clim=data.frame(clim[match(past.clim0$adcode,clim$ADCODE99),c("bio01","bio12")],
		past.clim0[,which(!colnames(past.clim0)%in%c("temp_0","prep_0"))])
	colnames(past.clim)[1:2]=c("temp_0","prep_0")
	library(data.table)
	biomat=as.data.frame(fread("Splev.new2.biomat.csv"))[,-1]
	biomat2=as.matrix(biomat[,-1])
	rownames(biomat2)=biomat[,1]		
	dis.clim=as.data.frame(fread("past.env.clim.csv"))# the Euclidean climate distance of pair-wised geounits since 150Ma with 1 Ma interval,温度标准化，降水*365后标准化
	dis.geo=as.data.frame(fread("pastdis2.csv"))[,-c(1,3,4)] # the cosdistance of pair-wised geounits since 150Ma with 1 Ma interval; cost value land:ocean=100:10
	#subset(clim,clim[,"bio01"]<=niche["bio01.90"]&clim[,"bio01"]>=niche["bio01.10"])
	disbar=function(niche){
		age=round(niche["age"])
		niche.dis.ans=subset(past.clim[,c("adcode",paste("temp",age,sep="_"),paste("prep",age,sep="_"))],						
						past.clim[,paste("temp",age,sep="_")]<=niche["ans.bio01.90"]&
						past.clim[,paste("temp",age,sep="_")]>=niche["ans.bio01.10"]&
						past.clim[,paste("prep",age,sep="_")]<=niche["ans.bio12.90"]&
						past.clim[,paste("prep",age,sep="_")]>=niche["ans.bio12.10"]
						)		
		niche.dis.rel=subset(past.clim[,c("adcode","temp_0","prep_0")],						
						past.clim[,"temp_0"]<=niche["bio01.90"]&
						past.clim[,"temp_0"]>=niche["bio01.10"]&
						past.clim[,"prep_0"]<=niche["bio12.90"]&
						past.clim[,"prep_0"]>=niche["bio12.10"]
						)						
		adcode.rel=niche.dis.rel$adcode
		adcode.ans=niche.dis.ans$adcode
		
		range.rel=length(adcode.rel)
		range.ans=length(adcode.ans)
		
		range.chg=ifelse(range.ans!=0&range.rel!=0,range.rel/range.ans,NA)#现在分布区占古分布区的比例		
		range.chg.rate=ifelse(range.ans!=0&range.rel!=0,abs(range.rel-range.ans)/niche["age"],NA)#潜在分布区的变化速率
		
		if (range.rel!=0&range.ans!=0){			
			#现在和古分布区重叠部分占现在和古分布区较小的那个的面积小于10%，认为现在和古分布区不重叠			
			range.overlap=length(intersect(adcode.ans,adcode.rel))/min(range.ans,range.rel) 
			f <- function(x, y) paste(x,y,sep="_")
			adcd=unique(c(as.vector(outer(adcode.ans, adcode.rel, f)),as.vector(outer(adcode.rel,adcode.ans, f))))				
			dis.geo.t=as.vector(as.matrix(dis.geo[dis.geo$id%in%adcd,paste(0:age,sep="")]))				
			costdis=quantile(dis.geo.t,0.05)
			dis.clim.t=na.omit(as.vector(as.matrix(dis.clim[dis.clim$ID%in%adcd,paste("clim_",0:age,sep="")])))				
			costdis.clim=quantile(dis.clim.t,0.05)						
		} else {costdis=NA;costdis.clim=NA;range.overlap=NA}
		
		temp=cbind(costdis,costdis.clim,range.rel,range.ans,range.overlap,range.chg,range.chg.rate)	
		# re=rbind(re,temp)}
		return(temp)
	}
	
	require(parallel)
	no_cores <- detectCores() - 1
	mycl <- makePSOCKcluster(no_cores); 
	clusterExport(cl = mycl, varlist = c("dis.clim","dis.geo","past.clim"))		
	biodis=t(parApply(mycl,biomat2,1,disbar))
	colnames(biodis)=c("dis.geo","dis.clim","range.rel","range.ans","range.overlap","range.chg","range.chg.rate")
	biodis2=cbind(biomat,biodis[match(biomat$genus,rownames(biodis)),])	
	write.csv(biodis2,"Splev.new2.biodis.csv")#tree use 100 radom trees of polytomy resolved tree from Luoao in 20201010.	
	stopCluster(mycl)	
	
##祖先分布区重建，参考Islands contribute disproportionately high amounts of evolutionary diversity in passerine birds ----	
	##assign states to tips	
	phname=read.csv("cluster.average.merge.csv")	
	library(data.table)	
	dis=as.data.frame(fread("Spdat_an_isrm_splev.csv"))
	rownames(dis)=dis[,1]
	dis=dis[,-1]
	dis.bio=cbind(phname=phname[,"phname"],dis[match(phname$adcode_ave2,rownames(dis)),])	
	require(parallel)
	no_cores <- detectCores() - 1
	mycl <- makePSOCKcluster(no_cores); 	
	f1=function(x,phname){
	 tmp1=tapply(x,phname,function(x){length(x[x>0])})#每个物种在各个生物地理区的geounits个数
	 tmp2=sort(tmp1/sum(tmp1),decreasing=T)#每个物种在各个生物地理区占该物种总分布区的比例
	 name=ifelse(tmp2[1]>=0.8,names(tmp2[1]),NA)
	 return(name)
	}
	sp2bio=na.omit(parApply(cl=mycl,dis.bio[,-1],2,f1,dis.bio[,1]))#214447 out of 235249,91.15%	
	
	##参考Islands contribute disproportionately...获得扩散routine	
	node2bio0=lapply(1:100,function(i,sp2bio,mycl){
		require(ape)
		require(castor)
		require(parallel)
		tre_ori=get(load(paste("trees/",i,".Rdata",sep="")))			
		d=subset(tre_ori$tip.label,!tre_ori$tip.label%in%names(sp2bio))
		tre <- drop.tip(tre_ori, tip=d)		
		status=sp2bio[tre$tip.label]		
		stus0=map_to_state_space(status)
		stus=stus0$mapped_states
		names(stus)=names(status)		
		ans=asr_max_parsimony(tre,stus)	
		node2bio=ans$ancestral_likelihoods
		colnames(node2bio)=names(stus0$name2index)
		re=unique(tre$edge[,1])
		re=re[order(re,decreasing=F)]
		rownames(node2bio)=re
		pos=which(tre$edge[,2]<=Ntip(tre))
		node2bio.t=node2bio[match(unique(tre$edge[pos,1]),rownames(node2bio)),]
		f2=function(x){
			x=sort(x,decreasing=T)
			ans.stat=ifelse(max(x)>=0.95,names(x[1]),NA)
		}
		node2bio2=parApply(cl=mycl,node2bio.t,1,f2)		
		node2bio3=node2bio2[match(tre$edge[pos,1],names(node2bio2))]
		names(node2bio3)=tre$tip.label[tre$edge[pos,2]]
		return(node2bio3)
	},sp2bio,mycl)
	stopCluster(mycl)	
	
	re=do.call(cbind,node2bio0)	
	mf=function(x){
		v=x[!is.na(x)]
		uniqv=unique(v)
		uniqv[which.max(tabulate(match(v,uniqv)))]	
	}
	node2bio=apply(re,1,mf)		
	biodis=as.data.frame(fread("Splev.new2.biodis.csv"))[,-1]
	biostat=data.frame(biodis,from=node2bio[match(biodis$genus,names(node2bio))],to=sp2bio[match(biodis$genus,names(sp2bio))])	
	write.csv(biostat,"Splev.new.biostat.csv")#最终分析用表 #tree use 100 radom trees of polytomy resolved tree from Luoao in 20201010.

## 进化速率与地理和气候距离的关系	
	#SEM
	library(data.table)
	biodis=as.data.frame(fread("Splev.new2.biodis.csv"))[,-1]
	library(lavaan)
	library(semPlot)
	vlist0=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")
	vlist=paste("rate",vlist0,sep=".")
	sem.data=abs(na.omit(biodis[,c("age",vlist,"dis.geo","dis.clim","range.overlap")]))	
	varl=c("T","Tmax","Tmin","P","Pmax","Pmin")
	colnames(sem.data)=c("age",varl,"geo","clim","range.overlap")
	sem.data2=sem.data[sem.data$range.overlap<=0.1,]
	par(mfrow = c(2,3),mar=c(0.5,0.5,0.5,0.5),oma=c(2,2,2,2))	
	for(i in 1:length(varl)){
		myModel<- paste(varl[i],"~clim+geo+age
		clim~~geo+age",sep="")
		fit <- sem(myModel, data = log(sem.data2))
		semPaths(fit,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex =3,label.cex =3,borders=FALSE,nCharNodes=0)#标准化的载荷std	
	}	
	#cor plot
	library(ggcorrplot)
	cord=log(sem.data2)
	corr <- round(cor(cord), 2)
	p.mat <- cor_pmat(cord)#sig
	ggcorrplot(corr, hc.order = FALSE, type = "lower", p.mat = p.mat,outline.color = "gray",insig = "blank",lab = TRUE)
	
	#scatter plot rate~geo+clim
	library(data.table)
	biodis=as.data.frame(fread("Splev.new2.biodis.csv"))[,-1]
	library(ggplot2);library(scales)
	vlist0=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")
	vlist=paste("rate",vlist0,sep=".")
	biop=abs(na.omit(biodis[,c(vlist,"dis.geo","dis.clim","range.overlap")]))	
	p <-ggplot(data=biop, mapping =aes(x=dis.geo,y=dis.clim,colour=rate.bio01.mean)) +
	scale_y_continuous(trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+
	scale_x_continuous(trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+		
		geom_point(size=0.6,alpha=0.5) +
		scale_colour_gradient(low = 'blue', high = 'red',trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+
		stat_density2d(aes(alpha=..density..), geom="tile",contour=FALSE,show.legend=FALSE)+
		theme(axis.title.x = element_blank(),
			axis.title.y = element_text(size=15),
			axis.text.x  = element_text(size=15),
			axis.text.y  = element_text(size=15),
			panel.grid.major =element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black"),
			legend.position=c(0.8,0.2))
## 地理和气候距离对进化速率的相对重要性随物种年龄和rangesize的变化	 ---	
	library(data.table)
	biostat=as.data.frame(fread("Splev.new.biostat.csv"))[,-1]
	vlist0=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")
	vlist=paste("rate",vlist0,sep=".")
##hp
	ager=do.call(rbind,lapply(1:length(vlist),function(j,biostat){	
		a=na.omit(biostat[,c("genus",vlist[j], "dis.geo","dis.clim","from","to","age")])		
		b=transform(a,age.cat = cut(a$age,breaks =50,include.lowest=TRUE))		
		b[,vlist[j]]=abs(b[,vlist[j]])		
		age=tapply(b$genus,b$age.cat,length)
		age=age[!is.na(age)&age>=3]
		#偏回归
		hpr=c()
		for (i in 1:length(age)){
			biostat.t=b[b$age.cat%in%names(age)[i],]
			re=hier.part::hier.part(log(biostat.t[,vlist[j]]),log(biostat.t[,c("dis.geo","dis.clim")]),gof = "Rsqu",barplot = FALSE)			
			temp=data.frame(evl=vlist[j],age=mean(as.numeric(unlist(strsplit(names(age)[i],",|[[]|[]]|[(]"))[-1])),
				geo=re$IJ$I[1]*100,clim=re$IJ$I[2]*100,joint=re$IJ$J[1]*100,clim.perc=re$I.perc[2,])
			hpr=rbind(hpr,temp)			
		}		
		return(hpr)	
	},biostat))
	
	ranger=do.call(rbind,lapply(1:length(vlist),function(j,biostat){	
		a=na.omit(biostat[,c("genus",vlist[j],"dis.geo","dis.clim","from","to","rangesize")])		
		b=transform(a,range.cat = cut(a$rangesize,breaks = 50,include.lowest=TRUE))			
		b[,vlist[j]]=abs(b[,vlist[j]])		
		rangesize=tapply(b$genus,b$range.cat,length)
		rangesize=rangesize[!is.na(rangesize)&rangesize>=3]		
		#偏回归
		hpr=c()
		for (i in 1:length(rangesize)){
			biostat.t=b[b$range.cat%in%names(rangesize)[i],]
			re=hier.part::hier.part(log(biostat.t[,vlist[j]]),log(biostat.t[,c("dis.geo","dis.clim")]),gof = "Rsqu",barplot = FALSE)
			temp=data.frame(evl=vlist[j],rangesize=mean(as.numeric(unlist(strsplit(names(rangesize)[i],",|[[]|[]]|[(]"))[-1])),
				geo=re$IJ$I[1]*100,clim=re$IJ$I[2]*100,joint=re$IJ$J[1]*100,clim.perc=re$I.perc[2,])
			hpr=rbind(hpr,temp)			
		}		
		return(hpr)	
	},biostat))
	
	##作图	
	library(ggpubr)	
	##堆积柱状图	
	hp=function(i,hpr){
		temp=data.frame(hpr[,1:2],rsq=hpr[,i+2],condition=colnames(hpr)[i+2])
		return(temp)
	}
	hprp=do.call(rbind,lapply(1:3,hp,ager))
	ggplot(hprp, aes(x=age, y=rsq, fill=condition))+geom_bar(stat="identity")+facet_wrap( ~ evl,ncol = 1,scales='free_y')+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
	
	hprp.rg=do.call(rbind,lapply(1:3,hp,ranger))
	ggplot(hprp.rg, aes(x=rangesize, y=rsq, fill=condition))+geom_bar(stat="identity")+facet_wrap( ~ evl,ncol = 1)+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
	
	##clim.perc
	ggplot(ager, aes(x=age, y=clim.perc))+geom_point()+
	geom_smooth(method = "gam",colour="red",fill="lightgray",size=1)+
	geom_hline(yintercept = 50,colour = "darkgray", linetype = "twodash", size = 1)+
	facet_wrap( ~ evl,ncol =3,scales='fixed')+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
	
	ggplot(ranger, aes(x=rangesize, y=clim.perc))+geom_point()+
	geom_smooth(method = "gam",colour="red",fill="lightgray",size=1)+
	geom_hline(yintercept = 50,colour = "darkgray", linetype = "twodash", size = 1)+
	facet_wrap( ~ evl,ncol =3,scales='fixed')+
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
	
## 各个生物地理区进化速率和扩散限制的分析	 ---		
	library(scales)
	library(data.table)
	biostat=as.data.frame(fread("Splev.new.biostat.csv"))[,-1]
	vlist0=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")
	vlist=paste("rate",vlist0,sep=".")
	cont=unique(na.omit(biostat$from))
	contype=do.call(rbind,lapply(1:length(cont),function(i,cont){
		do.call(rbind,lapply(1:length(cont),function(j,i,cont){
			data.frame(from=cont[i],to=cont[j],contype=paste(cont[i],cont[j],sep="/"))
		},i,cont))
	},cont))
	
	biostat.t=biostat#[biostat$range.overlap<=0.1,]
	bardat=do.call(rbind,lapply(1:length(vlist),function(i){
		a=na.omit(biostat.t[,c("genus",vlist[i],"dis.geo","dis.clim","from","to")])				
		a[,vlist[i]]=abs(a[,vlist[i]])			
		b=cbind(a,type=paste(a$from,a$to,sep="/"))		
		d=b[b$type%in%contype$contype,]	
		N=tapply(d$genus,d$type,length)
		mean=tapply(d[,vlist[i]],d$type,mean)				
		return(data.frame(varnm=vlist[i],type=names(N),N,mean))
	})) 	
	biostat2=cbind(biostat,type=paste(biostat$from,biostat$to,sep="/"))	
	nm=unique(biostat2[,c("from","to","type")])
	bardat2=cbind(bardat,nm[match(bardat$type,nm$type),1:2])	
	bardat2$from=factor(bardat2$from,levels=cont,ordered=TRUE)
	bardat2$to=factor(bardat2$to,levels=cont,ordered=TRUE)
	
	##气泡图显示各生物地理区扩散频率和年龄	
	bionode3=unique(biostat2[(!is.na(biostat2$from))&(!is.na(biostat2$to))&(biostat2$from!=biostat2$to),])
	shf.age=na.omit(cbind(mean=tapply(bionode3$age,bionode3$type,mean),median=tapply(bionode3$age,bionode3$type,median),N=tapply(bionode3$genus,bionode3$type,length)))
	shf.age=cbind(shf.age,bionode3[match(rownames(shf.age),bionode3$type),c("from","to","type")])
	shf.age$from=factor(shf.age$from,levels=cont,ordered=TRUE)
	shf.age$to=factor(shf.age$to,levels=cont,ordered=TRUE)	
	ggplot(shf.age, aes(x=from,y=to,size=N,fill=mean)) + geom_point(shape=21,stroke =0.5)+
		scale_size(range = c(6, 20), name="Freq. of shifts")+
			geom_abline(intercept=0,slope=1,colour = "darkgray", linetype = "twodash", size = 1)+			
			scale_fill_gradient(low="green",high="red")+	
			scale_x_discrete(labels=c("Ntp","Afr","Ori","Aus","Hol","S-A","Pat","N-Z"))+
			scale_y_discrete(labels=c("Ntp","Afr","Ori","Aus","Hol","S-A","Pat","N-Z"))+
			labs(x='Source',y='Sink',fill="Mean node age (Ma)") +
			geom_text(aes(y=to,label=N,hjust=0.5), size=3.5,color="black",position = position_dodge(width=0.00),check_overlap = FALSE)+theme(#panel.grid.major.y=element_blank(),
			panel.grid.minor.y=element_blank(),
			axis.text.y = element_text(angle = 90, hjust = 0.5),    
			axis.line.y=element_line(linetype=1,color='black'),
			axis.line.x=element_line(linetype=1,color='black'),
			axis.ticks = element_line(linetype=2,color='black'),
			panel.grid=element_line(linetype=2,color='grey'),
			panel.background = element_blank(),
			legend.background = element_rect(fill = NA),
			legend.text=element_text(face="bold",size=10),
			legend.title=element_text(face="bold",size=12),
			axis.text=element_text(face="bold",size=11.5),
			axis.title=element_text(face="bold",size=11.5))
			
	###1.各生物地理区扩散的时间格局 -----		
	###某一时期的扩散事件占该时期所有node的比例	
	extractNodeAge <- function(tree, node=NULL) {
		root.dist <- function(tree){				
			## find the root
			root.label <- unique(tree$edge[,1][!tree$edge[,1] %in% tree$edge[,2]])				
			tree.edge <- rbind(tree$edge,c(0,root.label))
			ii <- order(tree.edge[,2])
			tree.edge <- tree.edge[ii,]				
			## the results of root distance
			N.tip <- Ntip(tree)
			N.node <- tree$Nnode
			N.edge <- Nedge(tree)
			N.ii <- 1:(N.tip + N.node)
			rd <- numeric(N.tip + N.node)
			names(rd) <- c(tree$tip.label, tree$node.label)				
			## count the number of nodes from root to tips
			if (is.null(tree$edge.length)) edge.length <- rep(1, times = N.tip + N.node)
			if (!is.null(tree$edge.length))	edge.length <- c(tree$edge.length, 0)[ii]					
				for (i in 1:length(N.ii)){
					tip <- descendant <- N.ii[i]
					for (j in 1:N.node){
						rd <- rd + edge.length[descendant]
						if (tree.edge[descendant,1] == root.label | tree.edge[descendant,1] == 0){
							#rd[i] <- rd[i] + edge.length[root.label]
							break}
						ancestor <- tree.edge[descendant,1]
						descendant <- tree.edge[ancestor,2]
						}
					}					
			return(rd)
			}
		rd <- root.dist(tree)
		depth <- max(rd, na.rm=TRUE)
		## find the root
	    root.label <- unique(tree$edge[,1][!tree$edge[,1] %in% tree$edge[,2]])
		tree.edge <- rbind(tree$edge,c(0,root.label))
		ii <- order(tree.edge[,2])
		tree.edge <- tree.edge[ii,]
		edge.length <- c(tree$edge.length, 0)[ii]
		age.crown <- as.numeric(c(depth - rd)[Ntip(tree) + c(1:Nnode(tree))])
		age.stem <- as.numeric(c(depth - rd + edge.length)[Ntip(tree) + c(1:Nnode(tree))])
		if (is.null(node)) {
			result <- data.frame(Node=1:Nnode(tree), crown=age.crown, stem=age.stem)					   
			} else {
				if (is.character(node)) node.pos <- which(tree$node.label %in% node)
					else if (is.numeric(node)) node.pos <- which(tree.edge[,2] %in% c(node + Ntip(tree)))
				result <- data.frame(Node=c(1:Nnode(tree))[tree.edge[node.pos,2]-Ntip(tree)],
				crown=age.crown[tree.edge[node.pos,2]-Ntip(tree)],
				stem=age.stem[tree.edge[node.pos,2]-Ntip(tree)])					   
				}
		return(result)
		}

	library(ape)
	tre_ori=get(load(paste("trees/",1,".Rdata",sep="")))			
	d=subset(tre_ori$tip.label,!tre_ori$tip.label%in%biostat$genus)
	tre <- drop.tip(tre_ori, tip=d)	
	time=extractNodeAge(tre)
	time$Node=time$Node+Ntip(tre)	
	save(time, file="Splev_nodeAge.Rdata")
	time2=get(load("Splev_nodeAge.Rdata"))
	pos=which(tre$edge[,2]<=Ntip(tre))#提取edge中每个物种对应的编号所在的行号	
	crown=time[match(tre$edge[pos,1],time$Node),]
	rownames(crown)=tre$tip.label[tre$edge[pos,2]]
	contype2=c()	
	for(i in 1:length(cont)){
		for(j in 1: length(cont)){
		cont.t=data.frame(type=paste(cont[i],cont[j],sep="-"),
			contype[(contype[,"from"]%in%cont[i]&contype[,"to"]%in%cont[j])|(contype[,"to"]%in%cont[i]&contype[,"from"]%in%cont[j]),])
		contype2=rbind(contype2,cont.t)
		}
	}	
	thes=10
	bionode=cbind(crown,age.cat=round(crown$crown/thes)*thes,biostat2[match(rownames(crown),biostat2$genus),c("from","to","type")])
	bionode2=unique(bionode[(!is.na(bionode$from))&(!is.na(bionode$to))&(bionode$from!=bionode$to),])
	time2=cbind(time,age.cat=round(time$crown/thes)*thes)	
	all.n.node=tapply(time2[,"Node"],time2[,"age.cat"],length)	
	type.nm=as.character(unique(bionode2$type))	
	shf.node=do.call(rbind,lapply(1:length(type.nm),function(i,bionode,type.nm){
		sub.node=bionode2[bionode2$type%in%type.nm[i],]
		shf.n.node=tapply(sub.node[,"Node"],sub.node[,"age.cat"],length)		
		n.shf=shf.n.node[match(names(all.n.node),names(shf.n.node))]
		names(n.shf)=names(all.n.node)
		n.shf[is.na(n.shf)]=0
		pro.shf=n.shf/all.n.node*100		
		re=data.frame(age.cat=as.numeric(names(all.n.node)),n.shf,pro.shf,from=as.character(unique(sub.node$from)),to=as.character(unique(sub.node$to)),type.nm=type.nm[i])
		return(re)
	},bionode2,type.nm))	
	shf.p=cbind(shf.node,type=contype2[match(shf.node$type.nm,contype2$contype),"type"])
	#作图
	tye=as.character(unique(shf.p$type))
	plotn=lapply(1:length(tye),function(i,shf.p){
	shf.p2=shf.p[shf.p[,"type"]%in%tye[i],]		
	p=ggplot(shf.p2, aes(x=age.cat, y=pro.shf,group =type.nm,color=type.nm))+ geom_line(size=1)+
		#scale_y_continuous(trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+	
		theme(axis.title.x = element_blank(),
					axis.title.y = element_blank(),
					axis.text.x  = element_text(size=10),
					axis.text.y  = element_text(size=10),
					panel.grid.major.x =element_line(linetype=2,color='grey'),
					panel.grid.major.y =element_blank(),
					panel.grid.minor = element_blank(),
					panel.background = element_blank(),
					axis.line = element_line(colour = "black"),
					legend.position="top",
					legend.title=element_blank())+
					scale_x_reverse()+
					guides(color=guide_legend(nrow=2))
 	return(p)
	},shf.p)
	#pp=c();for(i in 1:length(plotn)) pp=paste(pp,"plotn[[",i,"]],",sep="")	
	figure=ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],plotn[[5]],plotn[[6]],plotn[[7]],plotn[[8]],plotn[[9]],plotn[[10]],plotn[[11]],plotn[[12]],plotn[[13]],plotn[[14]],plotn[[15]],plotn[[16]],plotn[[17]],plotn[[18]],plotn[[19]],plotn[[20]],plotn[[21]],plotn[[22]],plotn[[23]],plotn[[24]],plotn[[25]],plotn[[26]],plotn[[27]],plotn[[28]],
	nrow=4,ncol =7,widths=rep(2,4),heights=rep(2,7),hjust=0,align="hv",common.legend=FALSE)
	annotate_figure(figure, bottom = text_grob("Time before present (Ma)", size=10),left = text_grob("Pro of shfs in given time(ln trans)", rot = 90,size=10))
	
	#all
	shf.n.node=tapply(bionode2[,"Node"],bionode2[,"age.cat"],length)		
	n.shf=shf.n.node[match(names(all.n.node),names(shf.n.node))]
	names(n.shf)=names(all.n.node)
	n.shf[is.na(n.shf)]=0
	pro.shf=n.shf/all.n.node*100		
	re=data.frame(age.cat=as.numeric(names(all.n.node)),n.shf,pro.shf,from=as.character(unique(sub.node$from)),to=as.character(unique(sub.node$to)),type.nm=type.nm[i])
	p1=ggplot(re, aes(x=age.cat, y=pro.shf))+ geom_line(size=1)+ 			
		theme(axis.text.x  = element_text(size=10),
					axis.text.y  = element_text(size=10),
					panel.grid.major.x =element_line(linetype=2,color='grey'),
					panel.grid.major.y =element_blank(),
					panel.grid.minor = element_blank(),
					panel.background = element_blank(),
					axis.line = element_line(colour = "black"),
					legend.position="top",
					legend.title=element_blank())+
					scale_x_reverse()+labs(x='Time before present (Ma)',y='Pro of shfs in given time')
		#画下面的五个地质历史时期的矩形
		plt=plt+annotate("rect", xmin=100, xmax=65.5, ymin=-4, ymax=0,fill= "#B6DCB6")+		 
		 annotate("rect", xmin=65.5, xmax=55.8,ymin=-4, ymax=0,fill= "#D2E9E1")+
		  annotate("rect", xmin=55.8, xmax=33.9, ymin=-4, ymax=0,fill= "#FBEDC9")+
		  annotate("rect", xmin=33.9, xmax=23.03,ymin=-4, ymax=0,fill="#F8DDA9")+
		  annotate("rect", xmin=23.03, xmax=15, ymin=-4, ymax=0,fill="#FCB6D0")+		 
		   #画下面的五个地质历史时期的矩形
		  annotate("text", x =65.5+(100-65.5)/2 , y = -2 ,label = "K2",size=3)+
		  geom_vline(xintercept=65.5,col="gray",size=1,linetype="longdash",alpha=0.5)+
		  annotate("text", x =55.8+(65.5-55.8)/2, y =-2, label = "E1",size=3) +
		  geom_vline(xintercept=55.8,col="gray",size=1,linetype="longdash",alpha=0.5)+
		  annotate("text", x =33.9+(55.8-33.9)/2, y =-2, label = "E2",size=3) +
		  geom_vline(xintercept=33.9,col="gray",size=1,linetype="longdash",alpha=0.5)+
		  annotate("text", x =23.03+(33.9-23.03)/2, y = -2, label = "E3",size=3) +
		  geom_vline(xintercept=23.03,col="gray",size=1,linetype="longdash",alpha=0.5)+
		  annotate("text", x =15+(23.03-15)/2, y =-2, label = "N1",size=3)
		  
	re2=rbind(data.frame(nodetype="Biome-shift",re[,c("age.cat","n.shf")]),
	data.frame(nodetype="Node Num.",age.cat=as.numeric(names(all.n.node)),n.shf=all.n.node))
	p2=ggplot(re2, aes(x=age.cat, y=n.shf,group=nodetype,color=nodetype))+ geom_line(size=1)+ 			
		scale_y_continuous(trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+
		theme(axis.text.x  = element_text(size=10),
					axis.text.y  = element_text(size=10),
					panel.grid.major.x =element_line(linetype=2,color='grey'),
					panel.grid.major.y =element_blank(),
					panel.grid.minor = element_blank(),
					panel.background = element_blank(),
					axis.line = element_line(colour = "black"),
					legend.position=c(0.2,0.9),
					legend.title=element_blank())+
					scale_x_reverse()+labs(x='Time before present (Ma)',y='Freq of shfs in given time')
	ggarrange(p1,p2,nrow=1,ncol =2,widths=c(2),heights=c(2,2),labels = c("A", "B"),common.legend=FALSE)	
	
	##发生迁移的类群的统计（Family）
	library(data.table)
	biostat=na.omit(as.data.frame(fread("Splev.new.biostat.csv"))[,c("genus","age","from","to")])		
	taxonmic=unique(as.data.frame(fread("SpLevDis2.csv"))[,c("Family_E","Genus_E","Species_E1")])	
	taxa=cbind(taxonmic[match(biostat$genus,taxonmic$Species_E1),],biostat)
	#"Neotropical" "African" "Oriental" "Australian" "Holarctic" "Saharo-Arabian" "Patagonia" "NewZealand"
	fam.sp=tapply(taxonmic$Species_E1,taxonmic$Family_E,length)
	genus.sp=tapply(taxonmic$Species_E1,taxonmic$Genus_E,length)		
	bioshift.taxa=function(taxa,fam.sp,genus.sp,biofrom,bioto){
		taxa.t=subset(taxa,taxa$from==biofrom&taxa$to==bioto)
		if(dim(taxa.t)[1]==0){
		fam.re="没有发生迁移的物种；"
		}else{
			fam.sp.n=as.matrix(sort(tapply(taxa.t$Species_E1,taxa.t$Family_E,length),decreasing=T))
			colnames(fam.sp.n)=dim(fam.sp.n)[1]	#colnames is the number of family contained bioshift species
			fam.sp.pro=fam.sp.n[,1]/fam.sp[rownames(fam.sp.n)]*100
			fam=cbind(fam.sp.n,pro=round(fam.sp.pro[rownames(fam.sp.n)],1))
			n=0
			while (sum(fam[1:n,1])<=sum(fam[,1])*0.5){
			n=n+1
			fam2=fam[1:n,]
			}
			if (n==0) fam2=fam
			fam.re=paste("发生迁移的物种来自",colnames(fam2)[1],"个科，但超50%的物种来自",dim(fam2)[1],"个科,分别为",sep="")
			for (i in 1:dim(fam2)[1]){
				if(i==1){
					fam.re=paste(fam.re,rownames(fam2)[i],"（",fam2[i,1],"种，占科内总物种数的",fam2[i,2],"%），",sep="")
				}else{
					fam.re=paste(fam.re,rownames(fam2)[i],"（",fam2[i,1],"种，占",fam2[i,2],"%），",sep="")
				}
			}
			# #genus
			# genus.sp.n=as.matrix(sort(tapply(taxa.t$Species_E1,taxa.t$Genus_E,length),decreasing=T))
			# colnames(genus.sp.n)=dim(genus.sp.n)[1]	#colnames is the number of family contained bioshift species
			# genus.sp.pro=genus.sp.n[,1]/genus.sp[rownames(genus.sp.n)]*100
			# genus=cbind(genus.sp.n,pro=round(genus.sp.pro[rownames(genus.sp.n)],1))
			# m=0
			# while (sum(genus[1:m,1])<=sum(genus[,1])*0.5){
			# m=m+1
			# genus2=genus[1:m,]
			# }
			# if (m==0) genus2=genus
			# fam.re=paste(fam.re,"发生迁移的物种来自",colnames(genus2)[1],"个属，但超50%的物种来自",dim(genus2)[1],"个属,分别为",sep="")
			# for (i in 1:dim(genus2)[1]){
				# if(i==1){
					# fam.re=paste(fam.re,rownames(genus2)[i],"（",genus2[i,1],"种，占属内总物种数的",genus2[i,2],"%），",sep="")
				# }else{
					# fam.re=paste(fam.re,rownames(genus2)[i],"（",genus2[i,1],"种，占",genus2[i,2],"%），",sep="")
				# }					
			# }
		}
		return(fam.re)
	}
	phname=as.character(unique(biostat$from))
	re=""
	for (i in 1:(length(phname)-1)){
		for (j in (i+1):length(phname)){
			print(paste(phname[i],phname[j],sep=" "))
			intro=paste("从",phname[i],"区到",phname[j],"区，",sep="")
			dat=bioshift.taxa(taxa,fam.sp,genus.sp,phname[i],phname[j])	
			dat2=bioshift.taxa(taxa,fam.sp,genus.sp,phname[j],phname[i])
			re=paste(re,intro,dat,"相反方向上，",dat2,sep="")
		}
	}	
	write.table(re,"a.txt")#with genus
	write.table(re,"b.txt")#without genus
	
	###2.比较各生物地理区进化速率 ---	
	#bubble plot
	p=function(bardat2,v){
		bardat.p=bardat2[bardat2$varnm%in%v,]
		p=ggplot(bardat.p, aes(x=from,y=to,size=N,fill=mean)) + geom_point(shape=21,stroke =0.5)+
		scale_size(range = c(6, 25), name="Freq. of shifts",trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+
			geom_abline(intercept=0,slope=1,colour = "darkgray", linetype = "twodash", size = 1)+			
			scale_fill_gradient(low="green",high="red",trans = log_trans(),breaks = trans_breaks("log", function(x) exp(x)),labels = trans_format("log", math_format(e^.x)))+	
			scale_x_discrete(labels=c("Ntp","Afr","Ori","Aus","Hol","S-A","Pat","N-Z"))+
			scale_y_discrete(labels=c("Ntp","Afr","Ori","Aus","Hol","S-A","Pat","N-Z"))+
			labs(x='Source',y='Sink',fill="Niche evl rate") +
			geom_text(aes(y=to,label=N,hjust=0.5), size=3.5,color="black",position = position_dodge(width=0.00),check_overlap = FALSE)+ 			
			facet_wrap( ~ varnm,ncol = 2,scale="free")+#theme(strip.text= element_blank())+#theme_minimal()+
			theme(#panel.grid.major.y=element_blank(),
			panel.grid.minor.y=element_blank(),
			axis.text.y = element_text(angle = 90, hjust = 0.5),    
			axis.line.y=element_line(linetype=1,color='black'),
			axis.line.x=element_line(linetype=1,color='black'),
			axis.ticks = element_line(linetype=2,color='black'),
			panel.grid=element_line(linetype=2,color='grey'),
			panel.background = element_blank(),
			legend.background = element_rect(fill = NA),
			legend.text=element_text(face="bold",size=10),
			legend.title=element_text(face="bold",size=12),
			axis.text=element_text(face="bold",size=11.5),
			axis.title=element_text(face="bold",size=11.5))
	}
	p1=p(bardat2,"rate.bio01.mean")	
	p2=p(bardat2,"rate.bio12.mean")	
	ggarrange(p1,p2,nrow=1,ncol =2,widths=c(2),heights=c(2,2),labels = c("A", "B"),common.legend=FALSE,legend = "right")	
	
	###3.地理和气候距离对进化速率的相对重要性	 ---	
	hpr=do.call(rbind,lapply(1:length(vlist),function(j,biostat,contype){
	#for(j in 1:length(vlist)){
		a=na.omit(biostat[,c("genus",vlist[j],"dis.geo","dis.clim","from","to")])		
		a[,vlist[j]]=abs(a[,vlist[j]])		
		b=cbind(a,type=paste(a$from,a$to,sep="/"))
		c=sort(tapply(b$genus,b$type,length),decreasing=T)
		d=c[names(c)%in%contype$contype]	
		
		hpr=c();glmr=c()
		for (i in 1:length(d)){
			biostat.t=b[b$type%in%names(d)[i],]
			#偏回归
			re=hier.part::hier.part(log(biostat.t[,vlist[j]]),log(biostat.t[,c("dis.geo","dis.clim")]),gof = "Rsqu",barplot = FALSE)			
			#glm
			ss.glm <- function(r.glm)
                {
                r.ss <- summary(r.glm)
                rsq <- 100*(r.ss$null.deviance-r.ss$deviance)/r.ss$null.deviance
                adj.rsq <- 100*(1-(r.ss$deviance/r.ss$df.residual)/(r.ss$null.deviance/r.ss$df.null))
                f.stat <- ((r.ss$null.deviance-r.ss$deviance)/(r.ss$df.null-
                        r.ss$df.residual))/(r.ss$deviance/r.ss$df.residual)
                p <- pf(f.stat, r.ss$df.null-r.ss$df.residual, r.ss$df.residual, lower.tail=FALSE)

                return(c(r2=rsq,adj.r2=adj.rsq,p=p))
                }
			re.geo=glm(log(biostat.t[,vlist[j]])~log(biostat.t[,"dis.geo"]),family = "gaussian")
			re.clim=glm(log(biostat.t[,vlist[j]])~log(biostat.t[,"dis.clim"]),family = "gaussian")			
			sg=ifelse(ss.glm(re.geo)[3]<0.05,abs(as.numeric(coef(re.geo)[2])*sd(log(biostat.t[,"dis.geo"]))/sd(log(biostat.t[,vlist[j]]))),0)		
			sc=ifelse(ss.glm(re.clim)[3]<0.05,abs(as.numeric(coef(re.clim)[2])*sd(log(biostat.t[,"dis.clim"]))/sd(log(biostat.t[,vlist[j]]))),0)			
			temp=data.frame(evl=vlist[j],type=names(d)[i],geo=re$IJ$I[1]*100,clim=re$IJ$I[2]*100,joint=re$IJ$J[1]*100,
			from=unique(biostat.t$from),to=unique(biostat.t$to),clim.perc=re$I.perc[2,],freq=d[i],geo.glmr=ss.glm(re.geo)[2],clim.glmr=ss.glm(re.clim)[2],sg=sg,sc=sc)
			hpr=rbind(hpr,temp)						
		}	
		hpr=hpr[order(as.character(hpr$type)),]	
		return(hpr)	
	},biostat,contype))		
			
	hprp=na.omit(cbind(hpr[,c(1:2,6:length(hpr))],rsq=apply(hpr[,3:5],1,sum),
		glm.direct=ifelse(hpr$clim.perc>50,ifelse(hpr$clim.glmr>=0,"+","-"),ifelse(hpr$geo.glmr>=0,"+","-"))))	
	hprp$from=factor(hprp$from,levels=cont,ordered=TRUE)
	hprp$to=factor(hprp$to,levels=cont,ordered=TRUE)	
	##作图
	pp=function(hprp,v){
		hprp2=hprp[hprp$evl%in%v,]
		p=ggplot(hprp2, aes(x=from,y=to,size=rsq,fill=clim.perc)) + geom_point(shape=21,stroke =0.5,color="black")+
			scale_size(range = c(6, 25), name="Rsq (%)")+
			geom_abline(intercept=0,slope=1,colour = "darkgray", linetype = "twodash", size = 1)+
			scale_fill_gradient2(midpoint = 50,breaks = c(0,50,100),labels = c("geo","0",'clim'),low="#c15dd5", mid ="gray",high="#92d050",limits=c(0,100))+	scale_x_discrete(labels=c("Ntp","Afr","Ori","Aus","Hol","S-A","Pat","N-Z"))+
			scale_y_discrete(labels=c("Ntp","Afr","Ori","Aus","Hol","S-A","Pat","N-Z"))+
			labs(x='Source',y='Sink',fill="clim dorminance") +
			geom_text(aes(y=to,label=round(rsq,1),hjust=0.5), size=3.5,color="black",position = position_dodge(width=0.00),check_overlap = FALSE)+
			facet_wrap( ~ evl,ncol = 2,scale="free")+#theme(strip.text= element_blank())+#theme_minimal()+
			theme(#panel.grid.major.y=element_blank(),
			panel.grid.minor.y=element_blank(),
			axis.text.y = element_text(angle = 90, hjust = 0.5),    
			axis.line.y=element_line(linetype=1,color='black'),
			axis.line.x=element_line(linetype=1,color='black'),
			axis.ticks = element_line(linetype=2,color='black'),
			panel.grid=element_line(linetype=2,color='grey'),
			panel.background = element_blank(),
			legend.background = element_rect(fill = NA),
			legend.text=element_text(face="bold",size=10),
			legend.title=element_text(face="bold",size=12),
			axis.text=element_text(face="bold",size=11.5),
			axis.title=element_text(face="bold",size=11.5))
	}
	p1=pp(hprp,"rate.bio01.mean")	
	p2=pp(hprp,"rate.bio12.mean")	
	ggarrange(p1,p2,nrow=1,ncol =2,widths=c(2),heights=c(2,2),labels = c("A", "B"),common.legend=TRUE,legend = "right")	