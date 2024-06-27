###################################################################################################
################################## FORMAL analysis ################################################
###################################################################################################

####################################################################
#############correct names in spdis and in phylogeny using WCVP ####
####################################################################
library(data.table);library(dplyr)
dis=fread("Spdis_A&G.txt") #distribution data from Zhiheng at May 1, 2024
# dis$dis=1
# dis$Species_E2 <- iconv(dis$Species_E2, from = "ISO-8859-1", to = "UTF-8")%>%do.call(rbind,lapply(stringr::str_squish))#correct strange characters in species names
# write.table(dis,"Spdis_A&G.txt")
geo=read.csv("Geo-isrm.csv")#remove the islands(area<2.5 wan km2) and Antarctic	                            
dis=dis %>% left_join(geo,by=c("Adcode99"="ADCODE99")) %>% filter(Division=="A"&(!is.na(Lon))) %>% select(Adcode99,dis,Species_E2) %>%distinct%>%na.omit
# species name standard
spnames <- rWCVPdata::wcvp_names
splist.mat = dis%>%select(Species_E2)%>%distinct %>% left_join(spnames,by=c("Species_E2"="taxon_name"),relationship ="many-to-many") %>% filter(taxon_status%in%c("Accepted","Synonym")) %>% 
	dplyr::select(Species_E2,powo_id,plant_name_id,accepted_plant_name_id,taxon_rank,taxon_status,primary_author)
sp.unmat=dis%>%select(Species_E2)%>%distinct %>% filter(!Species_E2%in%unique(splist.mat$Species_E2)) 

library(TNRS)	
res=TNRS(taxonomic_names = sp.unmat$Species_E2,sources = "wcvp")
res.unmat=res%>%filter(Overall_score<0.9|Taxonomic_status=="No opinion"|Name_matched_rank%in%"genus"|Genus_score<1)
res2=TNRS(taxonomic_names = res.unmat$Name_submitted)%>%filter((!Taxonomic_status%in%"No opinion")&(!Name_matched_rank%in%"genus")&Genus_score==1)
res.all=rbind(res[!res$Name_submitted%in%res.unmat$Name_submitted,],res2)
#write.csv(res.all,"res.all.csv")#mannually check
sp.mat=splist.mat[order(splist.mat$taxon_status),]%>% .[!duplicated(.[,'Species_E2']),]%>%dplyr::select('Species_E2','taxon_status','accepted_plant_name_id')%>%
	left_join(spnames[,c('plant_name_id','taxon_name')],by=c('accepted_plant_name_id'='plant_name_id'))
res.all=res.all[,c("Name_submitted",'Taxonomic_status',"Accepted_name_id","Accepted_name")]
colnames(sp.mat)=colnames(res.all)
spname.cor=rbind(sp.mat,res.all)
dis=dis%>%left_join(spname.cor,by=c("Species_E2"="Name_submitted"))
save(dis,file="spdis_A_namecorred20240506.rda")

#phylogeny name standard
require(ape)
tre0=read.tree(paste("FromAo_Smith_100treesV2/",1,".tre",sep=""))
sp=data.frame(tip=tre0$tip.label,Species_E2=gsub("_"," ",tre0$tip.label))			
splist.mat = sp %>% left_join(spnames,by=c("Species_E2"="taxon_name"),relationship ="many-to-many") %>% filter(taxon_status%in%c("Accepted","Synonym")) %>% 
	dplyr::select(Species_E2,powo_id,plant_name_id,accepted_plant_name_id,taxon_rank,taxon_status,primary_author)
sp.unmat=dis%>%select(Species_E2)%>%distinct %>% filter(!Species_E2%in%unique(splist.mat$Species_E2)) 

library(TNRS)	
res=TNRS(taxonomic_names = sp.unmat$Species_E2,sources = "wcvp")
res.unmat=res%>%filter(Overall_score<0.9|Taxonomic_status=="No opinion"|Name_matched_rank%in%"genus"|Genus_score<1)
res2=TNRS(taxonomic_names = res.unmat$Name_submitted)%>%filter((!Taxonomic_status%in%"No opinion")&(!Name_matched_rank%in%"genus")&Genus_score==1)
res.all=rbind(res[!res$Name_submitted%in%res.unmat$Name_submitted,],res2)
#write.csv(res.all,"res.all.csv")#mannually check
sp.mat=splist.mat[order(splist.mat$taxon_status),]%>% .[!duplicated(.[,'Species_E2']),]%>%dplyr::select('Species_E2','taxon_status','accepted_plant_name_id')%>%
	left_join(spnames[,c('plant_name_id','taxon_name')],by=c('accepted_plant_name_id'='plant_name_id'))
res.all=res.all[,c("Name_submitted",'Taxonomic_status',"Accepted_name_id","Accepted_name")]
colnames(sp.mat)=colnames(res.all)
spname.cor=rbind(sp.mat,res.all)
tipnames=sp%>%left_join(spname.cor,by=c("Species_E2"="Name_submitted"))%>%filter(!(is.na(Accepted_name)|Accepted_name==""))
save(tipnames,file="FromAo_Smith_100treesV2/tipnames.rda")	
	
# caculate climatic niche for each species	
load("spdis_A_namecorred20240506.rda")	
clim=read.csv("vars.csv")[,c("ADCODE99","bio01","bio05","bio06","bio12","bio16","bio17")]#Worldclim	
spdat=dis %>% filter(!(is.na(Accepted_name)|Accepted_name==""))%>%select(Adcode99,Accepted_name,dis)%>% distinct %>%
	tidyr::pivot_wider(names_from = Accepted_name, values_from = dis,values_fill = 0) %>% tibble::column_to_rownames("Adcode99")

library(parallel)
no_cores <- detectCores() - 1
mycl <- makePSOCKcluster(no_cores);
xx0<-parLapply(cl=mycl, X=1:ncol(spdat),get.niche,spdat,clim)
xx1<-do.call(rbind,xx0)
stopCluster(mycl)
write.csv(xx1,"Climate_niche.csv")

################################################################################
#############Estimate absolute rates of niche evolution for each species #######
################################################################################
library(data.table);library(dplyr)
	load("FromAo_Smith_100treesV2/tipnames.rda")
	clim.niche=as.data.frame(fread("Climate_niche.csv"))
	vlist=c(paste(rep(c("bio01","bio05","bio06","bio12","bio16","bio17"),each=3),c(".mean",".90",".10"),sep=""),"SNBT","sppWLNBT","SNBP","sppWLNBP")
	status=clim.niche[,vlist]
	rownames(status)=clim.niche$V1	
	tipnames=tipnames%>%filter(Accepted_name%in%rownames(status))
	
get.niche.rate=function(i,status,vlist,tipnames){
		require(ape)
		tre0=read.tree(paste("FromAo_Smith_100treesV2/",i,".tre",sep=""))		#100 radom trees of polytomy resolved tree from Smith et al.	
		tre=drop.tip(tre0, tip=subset(tre0$tip.label,!tre0$tip.label%in%tipnames$tip))	
		tre$tip.label=tipnames[match(tre$tip.label,tipnames$tip),"Accepted_name"]
		status=as.matrix(status[match(tre$tip.label,rownames(status)),])	
		ans=do.call(cbind,lapply(1:length(vlist),function(j){
			stat=status[,vlist[j]]			
			ans=castor::asr_squared_change_parsimony(tre,stat,weighted = TRUE,check_input = TRUE)$ance
			return(ans)
		}))			
		colnames(ans)=paste("ans",vlist,sep=".")
		re=unique(tre$edge[,1])
		re=re[order(re,decreasing=F)]
		rownames(ans)=re
		
		pos=which(tre$edge[,2]<=Ntip(tre))
		genus=data.frame(X1=tre$edge[pos,1],X2=tre$edge[pos,2],genus=tre$tip.label[tre$edge[pos,2]],age=tre$edge.length[pos])		
		biomat=cbind(genus,status[match(genus$genus,rownames(status)),],ans[match(genus$X1,rownames(ans)),])				
		biomat2=biomat[order(biomat$genus),]
		return(biomat2)	
	}
	niche.rate.cal=function(biomat0){
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
		biomat2=data.frame(biomat,rate)	
		return(biomat2)
	}
		
	library(parallel)
	no_cores <- detectCores() - 1
	mycl <- makePSOCKcluster(no_cores); 	
	biomat0=parLapply(cl=mycl, X=1:100,get.niche.rate,status,vlist,tipnames)
	stopCluster(mycl)
	biomat2=niche.rate.cal(biomat0)	
	write.csv(biomat2,"Splev.new2.biomat.V2.csv")
	
	##estimate niche evl rate use mocecular smith tree
	#Smith phylogenetic tree
	# ref:Smith, S. A., and J. W. Brown. 2018. Constructing a broadly inclusive seed plant phylogeny. American Journal of Botany 105(3): 1–13.
	# download from: https://github.com/FePhyFoFum/big_seed_plant_trees
	tre0=read.tree("GBOTB.tre")
	tre=drop.tip(tre0, tip=subset(tre0$tip.label,!tre0$tip.label%in%tipnames$tip))	
	tre$tip.label=tipnames[match(tre$tip.label,tipnames$tip),"Accepted_name"]
	status2=status[match(tre$tip.label,rownames(status)),]
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
	
################################################################################	
##evaluate different evlution model
################################################################################
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
	library(data.table);library(dplyr)
	load("/blue/matthewthomas1/yunpeng.liu/bioshift/tipnames.rda")
	clim.niche=as.data.frame(fread("/blue/matthewthomas1/yunpeng.liu/bioshift/Climate_niche.csv"))
	vlist=c(paste(rep(c("bio01","bio05","bio06","bio12","bio16","bio17"),each=3),c(".mean",".90",".10"),sep=""),"SNBT","sppWLNBT","SNBP","sppWLNBP")
	status=clim.niche[,vlist]
	rownames(status)=clim.niche$V1	
	tipnames=tipnames%>%filter(Accepted_name%in%rownames(status))
	#tre0=read.tree("/blue/matthewthomas1/yunpeng.liu/bioshift/1.tre")
	
	no_cores <- detectCores() - 1
	mycl <- makePSOCKcluster(no_cores); 	
	biomat0=parLapply(cl=mycl, X=1:100,function(i,status,vlist,tipnames){
		require(ape)
		require(geiger)
		require(phytools)
		tre0=read.tree(paste("FromAo_Smith_100treesV2/",i,".tre",sep=""))		#100 radom trees of polytomy resolved tree from Smith et al.		
		tre=drop.tip(tre0, tip=subset(tre0$tip.label,!tre0$tip.label%in%tipnames$tip))	
		tre$tip.label=tipnames[match(tre$tip.label,tipnames$tip),"Accepted_name"]
		status=as.matrix(status[match(tre$tip.label,rownames(status)),])
		ans=do.call(cbind,lapply(1:length(vlist),function(j){
			stat=status[,vlist[j]]			
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
		biomat=cbind(genus,status[match(genus$genus,rownames(status)),],ans[match(genus$X1,rownames(ans)),])				
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
	write.csv(biomat2,"Splev.new2.biomat.V2.LA.csv")	
	

## correlation among different niche indexes ----
	#cor mat
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
	#PIC
	library(dplyr);library(ape);library(ggcorrplot);library(data.table)
	dis=as.data.frame(fread("SpLevDis2.csv"))
	clim=read.csv("vars.csv")
	#bio04 = Temperature Seasonality (standard deviation ×100)
	#bio15 = Precipitation Seasonality (Coefficient of Variation)	
	lat=dis%>%left_join(clim[,c("ADCODE99","bio04","bio15")],by=c("Adcode99"="ADCODE99"))%>%
		group_by(Species_E1)%>%summarize(Lat.mean=mean(abs(Lat),na.rm=T),TSN=mean(bio04,na.rm=T),PSN=mean(bio15,na.rm=T))
	
	niche=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),"mean",sep=".")
	niche.rate=paste("rate",niche,sep=".")
	niche.wid=c("SNBT","SNBP")
	biodis2=biodis[,c("genus","age",niche,niche.rate,niche.wid)]%>%filter(age>1)%>%
		left_join(lat,by=c("genus"="Species_E1"))
	biodis2[,niche.rate]=abs(biodis2[,niche.rate])
	rownames(biodis2)=biodis2$genus
	tre0=read.tree("FromAo_Smith_100treesV2/1.tre")		
	tre=drop.tip(tre0, tip=subset(tre0$tip.label,!tre0$tip.label%in%biodis2$genus))		
	tre$node.label<-NULL
	#Compute the phylogenetically independent contrasts using the method described by Felsenstein (1985).
	get.pic=function(i,biodis2,tre){
		var.t=biodis2[,i];names(var.t)=rownames(biodis2)
		pic.X <- pic(var.t, tre)
		return(pic.X)
	}
	biodis.pic=do.call(cbind,lapply(c("Lat.mean","TSN","PSN",niche.wid,niche,niche.rate),get.pic,biodis2,tre))
	colnames(biodis.pic)=c("Latitude","TSN","PSN","SNBT","SNBP","MAT","MTWM","MTCM","MAP","MPWQ","MPDQ",
		"MAT.rate","MTWM.rate","MTCM.rate","MAP.rate","MPWQ.rate","MPDQ.rate")
	
	#cor plot
	corr <- round(cor(biodis.pic, 2))
	p.mat <- cor_pmat(biodis.pic)#sig
	corr.t1<-corr["MAT.rate",c("Latitude","TSN","SNBT","MAT")];p.mat.t1=p.mat["MAT.rate",c("Latitude","TSN","SNBT","MAT")]
	corr.t2<-corr["MTWM.rate",c("Latitude","TSN","SNBT","MTWM")];p.mat.t2=p.mat["MTWM.rate",c("Latitude","TSN","SNBT","MTWM")]
	corr.t3<-corr["MTCM.rate",c("Latitude","TSN","SNBT","MTCM")];p.mat.t3=p.mat["MTCM.rate",c("Latitude","TSN","SNBT","MTCM")]
	
	corr.p1<-corr["MAP.rate",c("Latitude","PSN","SNBP","MAP")];p.mat.p1=p.mat["MAP.rate",c("Latitude","PSN","SNBP","MAP")]
	corr.p2<-corr["MPWQ.rate",c("Latitude","PSN","SNBP","MPWQ")];p.mat.p2=p.mat["MPWQ.rate",c("Latitude","PSN","SNBP","MPWQ")]
	corr.p3<-corr["MPDQ.rate",c("Latitude","PSN","SNBP","MPDQ")];p.mat.p3=p.mat["MPDQ.rate",c("Latitude","PSN","SNBP","MPDQ")]
	
	names(corr.t1)=names(p.mat.t1)=c("Latitude","ClimVar","NichWid","NichMean");
	names(corr.t2)=names(p.mat.t2)=c("Latitude","ClimVar","NichWid","NichMean");
	names(corr.t3)=names(p.mat.t3)=c("Latitude","ClimVar","NichWid","NichMean");
	names(corr.p1)=names(p.mat.p1)=c("Latitude","ClimVar","NichWid","NichMean");
	names(corr.p2)=names(p.mat.p2)=c("Latitude","ClimVar","NichWid","NichMean");
	names(corr.p3)=names(p.mat.p3)=c("Latitude","ClimVar","NichWid","NichMean");
	
	corr2=rbind(corr.t1,corr.t2,corr.t3,corr.p1,corr.p2,corr.p3);
	p.mat2=rbind(p.mat.t1,p.mat.t2,p.mat.t3,p.mat.p1,p.mat.p2,p.mat.p3);
	rownames(corr2)=rownames(p.mat2)=c("MAT.rate","MTWM.rate","MTCM.rate","MAP.rate","MPWQ.rate","MPDQ.rate")
	ggcorrplot(corr2, hc.order = FALSE, type = "full", p.mat = p.mat2,outline.color = "gray",insig = "blank",lab = TRUE,colors = c("#6D9EC1", "white", "#E46726"))		
	
	
	##based on grid cell Mean
	library(dplyr)	
	dis=as.data.frame(fread("SpLevDis2.csv"))
	clim=read.csv("vars.csv")
	
	mat.inc=biodis%>%filter(rate.bio01.mean>0&age>1)%>%select(genus,rate.bio01.mean);colnames(mat.inc)[2]="MAT.inc"
	mat.dec=biodis%>%filter(rate.bio01.mean<0&age>1)%>%select(genus,rate.bio01.mean);colnames(mat.dec)[2]="MAT.dec"
	map.inc=biodis%>%filter(rate.bio12.mean>0&age>1)%>%select(genus,rate.bio12.mean);colnames(map.inc)[2]="MAP.inc"
	map.dec=biodis%>%filter(rate.bio12.mean<0&age>1)%>%select(genus,rate.bio12.mean);colnames(map.dec)[2]="MAP.dec"
	biodis2=list(biodis,mat.inc,mat.dec,map.inc,map.dec) %>% purrr::reduce(left_join,by='genus') 
	rate.in.grid=dis%>%left_join(clim[,c("ADCODE99","bio04","bio15")],by=c("Adcode99"="ADCODE99"))%>%
		left_join(biodis2,by=c("Species_E1"="genus"))%>%filter(age>1)%>%group_by(Adcode99)%>%
			summarize(Adcode99=unique(Adcode99),Latitude=mean(abs(Lat),na.rm=T),
			TSN=mean(bio04,na.rm=T),PSN=mean(bio15,na.rm=T),
			SNBT=mean(SNBT),SNBP=mean(SNBP),MAT=mean(bio01.mean),MTWM=mean(bio05.mean),MTCM=mean(bio06.mean),
			MAP=mean(bio12.mean),MPWQ=mean(bio16.mean),MPDQ=mean(bio17.mean),
			MAT.rate=mean(abs(rate.bio01.mean)),MTWM.rate=mean(abs(rate.bio05.mean)),MTCM.rate=mean(abs(rate.bio06.mean)),
			MAP.rate=mean(abs(rate.bio12.mean)),MPWQ.rate=mean(abs(rate.bio16.mean)),MPDQ.rate=mean(abs(rate.bio17.mean)),
			MAT.rate.inc=mean(abs(MAT.inc),na.rm=T),MAT.rate.dec=mean(abs(MAT.dec),na.rm=T),
			MAP.rate.inc=mean(abs(MAP.inc),na.rm=T),MAP.rate.dec=mean(abs(MAP.dec),na.rm=T)
		)
	##SEM
	library(lavaan)
	library(semPlot)
	library(dplyr)
	par(mfrow = c(2,2),mar=c(0,0,0.4,0),oma=c(0,0,0,0))	
	#AIC=c()
	rt.var=c("MAT.rate","MTWM.rate","MTCM.rate","MAP.rate","MPWQ.rate","MPDQ.rate")
	rt.var=c("MAT.rate.inc","MAT.rate.dec","MAP.rate.inc","MAP.rate.dec")
	for(i in 1:length(rt.var)){
		if (i<=2) sem.data=rate.in.grid[,c("Latitude","TSN","SNBT",rt.var[i])] else 
			sem.data=rate.in.grid[,c("Latitude","PSN","SNBP",rt.var[i])]
		colnames(sem.data)=c("Latitude","ClimVar","NichWid","NichRate")
		myModel1="NichRate ~ ClimVar + NichWid 
			NichWid ~ ClimVar
			"					
		fit1 <- sem(myModel1, data = scale(sem.data))
		fit=fit1
		standardizedSolution(fit) %>% dplyr::filter(!is.na(pvalue)) %>% arrange(desc(pvalue)) %>% mutate_if("is.numeric","round",3) %>% select(-ci.lower,-ci.upper,-z)
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
	
	
	# GLM and scanter plot
	getscater.grid=function(x.list,y.list,xlab,ylab,dat.grid,show.lab=FALSE){
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
		width <- 40; height <- 25
		windows(width=width, height=height)
		n.col <- length(x.list); n.row <- length(y.list)
		mylayout <- layout(matrix(1:((2+n.col)*(2+n.row)), ncol=2+n.col, byrow=T), 
			width=c(0.7*width/(2+n.col), rep(width/(2+n.col),times = n.col), 0.3*width/(2+n.col)),
			height=c(0.3*height/(2+n.row), rep(height/(2+n.row),times=n.row), 0.7*height/(2+n.row)))
		layout.show((2+n.col)*(2+n.row))			
		for (i in 1:length(y.list)) {		
			if (i == 1) {
				for (j in 1:(2+n.col)) {par(mar=c(0,0,0,0)); plot.new()}
				}	
				for (xi in 1:length(x.list)) {
					temp2=na.omit(temp[,c(x.list[xi],y.list[i])])				
					x <- temp2[,1]
					y <- temp2[,2]			
					if (xi == 1) {par(mar=c(0,0,0,0)); plot.new()}
					par(mar=c(0.5,0.5,0,0), cex=1.2, cex.axis=1.2, cex.lab=1.2, mgp=c(2.4,0.3,0), tck=-0.04)
					plot(y~x, data=temp2, pch=19, xlab="", ylab="", cex=0.5, col='grey',axes=F); box()
					m <- glm(y~x,data=temp2)
					P=ss.glm(m)[3]
					Px=ifelse(P<=0.001,"***",
							ifelse(P>0.001&P<=0.01,"**",
							ifelse(P>0.01&P<=0.05,"*",
							#ifelse(P>0.05&P<0.1,"",
							ifelse(P>0.05,"ns",P))))
					a <- seq(min(x), max(x), (max(x)-min(x))/10)
					lines(x=a, y=predict(m, newdata=list(x=a)), col="black", lty=1, lwd=0.8)					
					mtext(substitute(paste(R^2, " = ",N^Px), list(N=round(ss.glm(m)[2],1),Px=Px)),
						side=3,adj=0.1,line=-2,cex=1.5,col="black",font=3)		
					if (xi == 1) {
						axis(side=2, label=T)						
						if(length(y.list)==1)mtext(side=2, text=ylab, line=1.6, cex=1.2, las=0)
						if(length(y.list)>1)mtext(side=2, text=ylab[[i]], line=1.6, cex=1.2, las=0)
					}
					if (i == length(y.list)) {
						axis(side=1, lwd = 0.5)						
						if (show.lab) mtext(side=1, text=xlab[xi], line=1.6, cex=1.2, las=0)
					}
					if (xi == length(x.list)) {par(mar=c(0,0,0,0)); plot.new()}					
				}
			}			
		}
			
	xlab1=c("Latitude","Niche position")
	xlab2=c("Seasonality","Niche Width")
	getscater.grid(c("Latitude","MAT"),"MAT.rate",xlab1,substitute(paste('MAT Rate (°C ',Myr^-1,")")),rate.in.grid,show.lab=TRUE)
	getscater.grid(c("Latitude","MAP"),"MAP.rate",xlab1,substitute(paste('MAP Rate (mm ',Myr^-1,")")),rate.in.grid,show.lab=TRUE)	
	getscater.grid(c("TSN","SNBT"),"MAT.rate",xlab2,substitute(paste('MAT Rate (°C ',Myr^-1,")")),rate.in.grid,show.lab=TRUE)
	getscater.grid(c("PSN","SNBP"),"MAP.rate",xlab2,substitute(paste('MAP Rate (mm ',Myr^-1,")")),rate.in.grid,show.lab=TRUE)
		
	getscater.grid(c("Latitude","MTWM"),"MTWM.rate",xlab1,substitute(paste('MTWM Rate (°C ',Myr^-1,")")),rate.in.grid,show.lab=TRUE)
	getscater.grid(c("Latitude","MTCM"),"MTCM.rate",xlab1,substitute(paste('MTCM Rate (°C ',Myr^-1,")")),rate.in.grid,show.lab=TRUE)	
	getscater.grid(c("Latitude","MPWQ"),"MPWQ.rate",xlab1,substitute(paste('MPWQ Rate (mm ',Myr^-1,")")),rate.in.grid,show.lab=TRUE)
	getscater.grid(c("Latitude","MPDQ"),"MPDQ.rate",xlab1,substitute(paste('MPDQ Rate (mm ',Myr^-1,")")),rate.in.grid,show.lab=TRUE)
	
	getscater.grid(c("Latitude"),c("MAT.rate.inc","MAT.rate.dec","MAP.rate.inc","MAP.rate.dec"),xlab1,
		c(substitute(paste('δT+ Rate (°C ',Myr^-1,")")),substitute(paste('δT- Rate (°C ',Myr^-1,")")),
			substitute(paste('δP+ Rate (°C ',Myr^-1,")")),substitute(paste('δP- Rate (°C ',Myr^-1,")"))),
		rate.in.grid,show.lab=TRUE)
	
	getscater.grid(c("TSN","SNBT"),c("MAT.rate.inc","MAT.rate.dec"),xlab2,
		c(substitute(paste('δT+ Rate (°C ',Myr^-1,")")),substitute(paste('δT- Rate (°C ',Myr^-1,")"))),
		rate.in.grid,show.lab=TRUE)	
	getscater.grid(c("PSN","SNBP"),c("MAP.rate.inc","MAP.rate.dec"),xlab2,
		c(substitute(paste('δP+ Rate (mm ',Myr^-1,")")),substitute(paste('δP- Rate (mm ',Myr^-1,")"))),
		rate.in.grid,show.lab=TRUE)
	
	getscater.grid(c("TSN","SNBT"),c("MTWM.rate","MTCM.rate"),xlab2,
		c(substitute(paste('MTWM Rate (°C ',Myr^-1,")")),substitute(paste('MTCM Rate (°C ',Myr^-1,")"))),
		rate.in.grid,show.lab=TRUE)	
	getscater.grid(c("PSN","SNBP"),c("MPWQ.rate","MPDQ.rate"),xlab2,
		c(substitute(paste('MPWQ Rate (mm ',Myr^-1,")")),substitute(paste('MPDQ Rate (mm ',Myr^-1,")"))),
		rate.in.grid,show.lab=TRUE)	

	
## niche evolutionary patterns ----
	library(data.table)	;library(dplyr)
	#biodis=as.data.frame(fread("Splev.new3.biomat.csv"))[,-1]
	biodis=as.data.frame(fread("Splev.new2.biomat.V2.smith.csv"))[,-1] #estimate niche evl rate use mocecular smith tree
	biodis=as.data.frame(fread("Splev.new3.biomat.V2.csv"))[,-1]#estimate niche evl rate BM
	biodis=as.data.frame(fread("Splev.new3.biomat.V2.LA.csv"))[,-1]#estimate niche evl rate LA
	rownames(biodis)=biodis$genus			
	dis=get(load("spdis_A_namecorred20240506.rda"))
	niche=paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep="")
	niche.rt=paste("rate",niche,sep=".")
	library(sf)
	map1 <- read_sf('PhyloRealms.new/PhyloRealms.new.shp')
	
	#plot	
	rate.p=function(biodis,rate,ratetype,dis,map1,main){
		#rate=niche.rt[1];ratetype="inc";main="MAT"
		biodis=biodis[biodis$age>=1,]#182277 out of 231567	
		bios0=biodis[,rate]
		names(bios0)=rownames(biodis)
		if (rate%in%c("lifefrom","to.tropic","genus")&ratetype%in%"sprich"){				
			bios=bios0
			rate.dis=na.omit(cbind(dis,rate=bios[match(dis$Accepted_name,names(bios))]))
			if (rate%in%"genus") {nich.rate=tapply(rate.dis$rate,rate.dis$Adcode99,length)} else {
				pro.f=function(x) length(x[x%in%"W"|x%in%"tro"])/(length(x[x%in%"W"|x%in%"tro"])+length(x[x%in%"H"|x%in%"tmp"]))
				nich.rate=tapply(rate.dis$rate,rate.dis$Adcode99,pro.f)
			}									
		} else{			
			if(ratetype=="inc") bios=bios0[bios0>0]#current>ans
			if(ratetype=="dec")	bios=abs(bios0[bios0<0])#current<ans
			if(ratetype=="all")	bios=abs(bios0)				
			rate.dis=na.omit(cbind(dis,rate=bios[match(dis$Accepted_name,gsub("_"," ",names(bios)))]))		
			nich.rate=tapply(rate.dis$rate,rate.dis$Adcode99,mean)
		}	
		range01 <- function(x){(x-min(x))/(max(x)-min(x))}
		map2<-cbind(map1,nich.rate=nich.rate[match(map1$GeoID,names(nich.rate))]) 
			
		mdat.p=as.data.frame(map2)[,c("Lat","nich.rate")]
		mdat.p$nich.rate2=183+range01(mdat.p$nich.rate)*40
		n=ifelse(rate%in%c("rate.bio01.mean","rate.bio05.mean","rate.bio06.mean"),1,0)	
		p<-ggplot()+
		geom_rect(aes(ymin=min(map1$Lat)-10,ymax=max(map1$Lat)+10,
					  xmin=min(map1$Lon)-30,xmax=max(map1$Lon)+10),fill='#F2F7FB',col='gray75')+##E5F1FB
		geom_rect(aes(ymin=min(map1$Lat)-10,ymax=max(map1$Lat)+10,
					  xmin=max(map1$Lon)+10,xmax=max(map1$Lon)+55),fill='transparent',col='gray75')+
		geom_sf(data=map2,aes(x=Lon,y=Lat,fill=nich.rate),color="transparent")+ 
		scale_fill_gradient2(midpoint=mean(nich.rate),low="blue",mid="yellow",high="red",name=main)+
		geom_point(data=mdat.p,aes(y=Lat,x=nich.rate2,fill=nich.rate),shape=21,size=1.5,pch=1,show.legend = F)+
		# scale_y_continuous(expand = c(0,0),breaks = seq(-60,90,20),
						   # labels = paste0(seq(-60,90,20),"°"),position = "right")+
		scale_x_continuous(expand = c(0,0),breaks = 183+c(14,28),
						   labels = c(round((max(nich.rate)-min(nich.rate))/3,n),round(2*(max(nich.rate)-min(nich.rate))/3,n)),
						   position = 'bottom')+
		theme(#aspect.ratio=(diff(range(test_df$y))+5)/(diff(range(test_df$x))+10*max(nich.rate)+2),#+75*max(srlist$P)
			  #text = element_text(size =10),
			  panel.background=element_rect(fill='transparent'),
			  legend.background=element_rect(fill='transparent'),
			  panel.border = element_blank(),
			  axis.title = element_blank(),
			  axis.text = element_text(face="bold",size =15),
			  axis.text.x = element_text(angle=30),
			  axis.ticks = element_blank(),
			  panel.grid.minor = element_blank(),
			  #panel.grid.major =  element_line(colour = "darkgray", size = 0.5),
			  legend.text=element_text(face="bold",size=15),
			  legend.title=element_text(face="bold",size=18),
			  legend.position = c(0.09,0.38)) 	
		return(p)
	}
	
	#fig.1
	p.mat=rate.p(biodis,niche.rt[1],ratetype="all",dis,map1,main="MAT")	
	p.map=rate.p(biodis,niche.rt[4],ratetype="all",dis,map1,main="MAP")
	ggarrange(p.mat,p.map,nrow=2,ncol =1,labels=c("a","b"),font.label = list(size = 20))

	#fig.2
	#bio01
	p2=rate.p(biodis,niche.rt[1],ratetype="inc",dis,map1,main="δT+")
	p3=rate.p(biodis,niche.rt[1],ratetype="dec",dis,map1,main="δT-")
	#bio12
	p5=rate.p(biodis,niche.rt[4],ratetype="inc",dis,map1,main="δP+")
	p6=rate.p(biodis,niche.rt[4],ratetype="dec",dis,map1,main="δP-")
	ggarrange(p2,p5,p3,p6,nrow=2,ncol =2,labels="auto",font.label = list(size = 30))

	
	#bio05 Max Temperature of Warmest Month
	p1=rate.p(biodis,niche.rt[2],ratetype="all",dis,cordacd,map1,main="MTWM")
	#bio16 Precipitation of Wettest Quarter
	p2=rate.p(biodis,niche.rt[5],ratetype="all",dis,cordacd,map1,main="MPWQ")
	#bio06
	p3=rate.p(biodis,niche.rt[3],ratetype="all",dis,cordacd,map1,main="MTCM")
	#bio17
	p4=rate.p(biodis,niche.rt[6],ratetype="all",dis,cordacd,map1,main="MPDQ")
	ggarrange(p.mat,p.map,p1,p2,p3,p4,
		nrow=3,ncol =2,labels=c("a","b","c","d"),font.label = list(size = 20))

	
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
	p1=rate.p(biodis,"SNBT",ratetype="all",dis,cordacd,map1,main="Temperature\nwidth")
	p2=rate.p(biodis,"SNBP",ratetype="all",dis,cordacd,map1,main="Precipitation\nwidth")
	p3=rate.p(biodis,"bio01.mean",ratetype="all",dis,cordacd,map1,main="NichMean\n(Thermal)")
	p4=rate.p(biodis,"bio12.mean",ratetype="all",dis,cordacd,map1,main="NichMean\n(Hydrologic)")
	ggarrange(p1,p2,ncol=1,labels=c("a","b"),label.x=0.03,font.label = list(size = 20))
	
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
	#plot.clim(shape,shape3,clim.var="ART",main="ART","b")
	plot.clim(shape,shape3,clim.var="PSN",main="PSN","b")
	
###################################################################################################
################################## compare niche using data of this study and data from liu et al. ################################################
###################################################################################################
dis20km=get(load("dis.china.20km.Rdata"))
clim20km=read.csv("grid20km/climate.csv")[,c("ADCODE99","bio01","bio05","bio06","bio12","bio16","bio17")]
disGSU=get(load("dis.china.gsu.Rdata"))
climGSU=read.csv("vars.csv")[,c("ADCODE99","bio01","bio05","bio06","bio12","bio16","bio17")]#Worldclim from luoao
get.niche=function(i,dis,clim){
		x=dis[,i]
		adcd=rownames(dis)[which(as.logical(x))]	
		n<-length(adcd)
		data1=clim[clim$ADCODE99%in%adcd,]
		
		# climate mean
		AGG <- sapply(data1[,c("bio01","bio05","bio06","bio12","bio16","bio17")],function(x)
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
}	

nicheGSU=do.call(rbind,lapply(1:ncol(disGSU),get.niche,disGSU,climGSU))
niche20km=do.call(rbind,lapply(1:ncol(dis20km),get.niche,dis20km,clim20km))
nicheGSU$Acname=rownames(nicheGSU);niche20km$Acname=rownames(niche20km)
save(nicheGSU,file="niche.china.gsu.Rdata")
save(niche20km,file="niche.china.20km.Rdata")#4381 spp.

# compare niche between GSU and 20km
	library(data.table)
	library(dplyr)	
	vlist=paste(rep(c("bio01","bio05","bio06","bio12","bio16","bio17"),each=4),c("mean","median","90","10"),sep=".")
	var.lab=data.frame(vlist=vlist,lab=paste(rep(c("MAT","MTWM","MTCM","MAP","MPWQ","MPDQ"),each=4),c("mean","median","90 quant.","10 quant."),sep="."))
	datall=c()
	for (i in vlist){
		dat=nicheGSU[,c("Acname",i)]%>%left_join(niche20km[,c("Acname",i)],by="Acname")%>%distinct()
		dat$var.type=var.lab[var.lab$vlist==i,"lab"]
		colnames(dat)[2:3]=c("GSU","grid20km")
		datall=rbind(datall,na.omit(dat))
	}
## compare niche with Liu et al. NEE
	library(data.table);library(caper)
	clim.niche=as.data.frame(fread("Climate_niche.csv"))
	rownames(clim.niche)=clim.niche$V1
	#liu et al NEE
	niche.nee=read.csv("Dataset_S6.csv")
	niche.nee$bio1.mean.all=clim.niche[match(niche.nee$Species,rownames(clim.niche)),"bio01.mean"]
	#write.csv(niche.nee,"Dataset_S6.csv") #手工检验未匹配的物种,808 angiosperms and 709 matched
	tre0=read.tree("FromAo_Smith_100treesV2/1.tre")		
	tre0$node.label=NULL
	tre=keep.tip(tre0,tip=tre0$tip[tre0$tip%in%niche.nee$Species])
	
	# vlist=paste(rep(c("bio01","bio05","bio06","bio12","bio16","bio17"),each=4),c("mean","median","90","10"),sep=".")	
	# var.lab=data.frame(vlist=vlist,lab=paste(rep(c("MAT","MTWM","MTCM","MAP","MPWQ","MPDQ"),each=4),c("mean","median","90 quant.","10 quant."),sep="."))
	vlist=c("bio01.mean","bio12.mean")	
	var.lab=data.frame(vlist=vlist,lab=c("MAT","MAP"))
	datall=c();rsq=c()
	for (i in vlist){
		dat=cbind(var.type=var.lab[var.lab$vlist==i,"lab"],niche.nee[,c("Species","Clade",i)],clim.niche[match(niche.nee$Species,rownames(clim.niche)),i])
		colnames(dat)[4:5]=c("Nee","thispaper")
		dat=na.omit(dat)
		dat.phy <- comparative.data(tre, dat, Species, vcv=TRUE, vcv.dim=3)
		M <- pgls(Nee ~ thispaper, data=dat.phy, lambda='ML')
		#M=lm(Nee~thispaper,data=dat)
		dat.phy[[2]]$resid=resid(M)
		datall=rbind(datall,dat.phy[[2]])
		rsq=c(rsq,summary(M)$r.sq)
	}
	
	library(ggpubr)
	theme=theme(axis.text = element_text(size=12,color='black',angle=0,hjust=1),
		axis.title = element_text(size=15,color='black',angle=0),
		panel.grid.minor.y=element_blank(),
		axis.line.y=element_line(linetype=1,color='black'),
		axis.line.x=element_line(linetype=1,color='black'),
		axis.ticks = element_line(linetype=2,color='black'),
		panel.grid=element_line(linetype=2,color='grey'),
		panel.background = element_blank())
	pdat.mat=datall[datall$var.type%in%"MAT",]
	p1=ggplot(pdat.mat,aes(x=Nee, y=thispaper))+
	labs(x='MAT Niche in Liu et al. (°C) ',
		y='MAT Niche in this study (°C) ')+#,fill="Clade in \nLiu et al. (2022)") +	
	geom_abline(intercept=0,slope=1, linetype = "twodash", col="red",size=1)+
	geom_point(size=2,shape=21,fill="#00BD5F",color="black",alpha=0.5)+
	annotate("text", x=20 , y= -5,size=6,label=substitute(paste(R^2, " = ",N^"***"), list(N=round(rsq[1]*100,2))))+
	theme
	p2=ggplot(pdat.mat, aes(x=resid))+
		geom_histogram(color="darkblue", fill="lightblue")+
		geom_vline(aes(xintercept=0),color="blue", linetype="dashed", size=1)+
		labs(x="Residuals of MAT Niche",y="Frequency")+
		theme	
	pdat.map=datall[datall$var.type%in%"MAP",]
	p3=ggplot(pdat.map,aes(x=Nee, y=thispaper))+
	labs(x="MAP Niche in Liu et al. (mm)",y="MAP Niche in this study (mm)")+#,fill="Clade in \nLiu et al. (2022)") +	
	geom_abline(intercept=0,slope=1, linetype = "twodash", col="red",size=1)+
	geom_point(size=2,shape=21,fill="#00BD5F",color="black",alpha=0.5)+
	annotate("text", x=3500 , y= 200,size=6,label=substitute(paste(R^2, " = ",N^"***"), list(N=round(rsq[2]*100,2))))+
	theme
	p4=ggplot(pdat.mat, aes(x=resid))+
		geom_histogram(color="darkblue", fill="lightblue")+
		geom_vline(aes(xintercept=0),color="blue", linetype="dashed", size=1)+
		labs(x="Residuals of MAP Niche",y="Frequency")+
		theme
	ggarrange(p1,p2,p3,p4,
		nrow=2,ncol =2,widths=c(1,1),heights=c(1,1),labels="auto",font.label = list(size = 20))

# compare niche evl rate
	library(dplyr)
	niche20km=get(load("niche.china.20km.Rdata"))
	nicheGSU=get(load("niche.china.gsu.Rdata"))
	vlist=c(paste(c("bio01","bio05","bio06","bio12","bio16","bio17"),".mean",sep=""))
	status20km=niche20km[,vlist];rownames(status20km)=gsub(" ", "_",rownames(status20km))
	statusGSU=nicheGSU[,vlist];rownames(statusGSU)=gsub(" ", "_",rownames(statusGSU))
	
	niche.nee=read.csv("Dataset_S6.csv")	
	status.nee=niche.nee[,vlist];rownames(status.nee)=gsub(" ", "_",niche.nee$Species)
	
	get.niche.rate=function(i,status,vlist){
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
	}
	niche.rate.cal=function(biomat0){
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
		biomat2=data.frame(biomat,rate)	
		return(biomat2)
	}
	
	library(parallel)
	no_cores <- detectCores() - 1
	mycl <- makePSOCKcluster(no_cores); 	
	biomat0.gsu=parLapply(cl=mycl, X=1:100,get.niche.rate,statusGSU,vlist)
	biomat0.20km=parLapply(cl=mycl, X=1:100,get.niche.rate,status20km,vlist)
	biomat0.nee=parLapply(cl=mycl, X=1:100,get.niche.rate,status.nee,vlist)
	stopCluster(mycl)
	biomat.gsu=niche.rate.cal(biomat0.gsu)
	biomat.20km=niche.rate.cal(biomat0.20km)
	biomat.nee=niche.rate.cal(biomat0.nee)
	save(biomat.gsu,file="biomat.china.gsu.Rdata")
	save(biomat.20km,file="biomat.china.20km.Rdata")
	save(biomat.nee,file="biomat.nee.Rdata")
	
	#plot niche evl Rate
	library(rgdal);library(ggpubr);library(dplyr);library(sf)	
	biomat.gsu=get(load("biomat.china.gsu.Rdata"))
	biomat.20km=get(load("biomat.china.20km.Rdata"))
	dis20km=get(load("dis.china.20km.Rdata"))
	disGSU=get(load("dis.china.gsu.Rdata"))
	dis20km2=data.frame(Acname=rep(colnames(dis20km),each=nrow(dis20km)),
		Adcode99=rep(rownames(dis20km),ncol(dis20km)),dis=as.vector(dis20km)) %>%filter(dis>0) %>%select(-dis)
	disGSU2=data.frame(Acname=rep(colnames(disGSU),each=nrow(disGSU)),
		Adcode99=rep(rownames(disGSU),ncol(disGSU)),dis=as.vector(disGSU)) %>%filter(dis>0) %>%select(-dis)
	
	map.gsu <- readOGR(dsn = "county_islandrm.dbf",stringsAsFactors=FALSE);
	map1.gsu<-ggplot2::fortify(map.gsu)
	cordacd.gsu<-as.data.frame(map.gsu@data)%>%mutate(ADCODE99=as.numeric(ADCODE99));
	cordacd.gsu$id<-as.character(0:(length(unique(map1.gsu$id))-1))
	
	map20km <- st_read("grid20km/grid20kmChina.shp") %>%st_transform(st_crs(map.gsu)$proj4string)%>%as('Spatial')
	colnames(map20km@data)[c(2,4:5)]=c("ADCODE99","Lon","Lat")
	map1.20km<-ggplot2::fortify(map20km)
	cordacd20km<-as.data.frame(map20km@data)%>%mutate(ADCODE99=as.numeric(ADCODE99));
	cordacd20km$id<-unique(map1.20km$id)
	
	rate="rate.bio12.mean";biodis=biomat.gsu;dis=disGSU2	;main="MAT GSU";cordacd=cordacd.gsu;map1=map1.gsu
	rate="rate.bio01.mean";biodis=biomat.20km;dis=dis20km2	;main="MAT";cordacd=cordacd20km;map1=map1.20km
	rate.p=function(biodis,rate,dis,cordacd,map1,main,na.rm=T){
		#rate=niche.rt[1];ratetype="inc"
		#biodis=biodis[biodis$age>=1,]#3394 spp.
		bios0=biodis[,rate]
		names(bios0)=rownames(biodis)
				
		bios=abs(bios0)				
		rate.dis=cbind(dis,rate=bios[match(dis$Acname,gsub("_", " ",names(bios)))])		
		nich.rate=tapply(rate.dis$rate,rate.dis$Adcode99,mean,na.rm=T)
		
		range01 <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}			
		mdat=cbind(cordacd,nich.rate=nich.rate[match(cordacd$ADCODE99,names(nich.rate))])
		map2<-inner_join(map1,mdat)
		map2=map2[map2$lat>15,]
		if(na.rm==T) {
			map2=map2[!is.na(map2$nich.rate),]
			mdat=mdat[!is.na(mdat$nich.rate),]
			}
			
		mdat.p=mdat[,c("Lat","nich.rate")]
		mdat.p$nich.rate2=max(map2$long)+range01(mdat.p$nich.rate)*8
		n=ifelse(rate%in%c("rate.bio01.mean","rate.bio05.mean","rate.bio06.mean"),1,0)			
		p<-ggplot()+
		geom_rect(aes(ymin=min(map2$lat)-5,ymax=max(map2$lat)+5,
					  xmin=min(map2$long),xmax=max(map2$long)),fill='#F2F7FB',col='gray75')+##E5F1FB
		geom_rect(aes(ymin=min(map2$lat)-5,ymax=max(map2$lat)+5,
					  xmin=max(map2$long),xmax=max(map2$long)+9),fill='transparent',col='gray75')+
		geom_polygon(data=map2,aes(x=long,y=lat,group=group,fill=nich.rate),color="transparent")+ 
		scale_fill_gradient2(midpoint=10,low="blue",mid="yellow",high="red",name=main,na.value = "grey50")+
		geom_point(data=mdat.p,aes(y=Lat,x=nich.rate2,fill=nich.rate),shape=21,size=1.5,pch=1,show.legend = F)+
		scale_y_continuous(expand = c(0,0),breaks = seq(10,60,15),
						   labels = paste0(seq(10,60,15),"°"),position = "right")+
		scale_x_continuous(expand = c(0,0),breaks = max(map2$long)+c(3,6),
						   labels = c(round((max(nich.rate,na.rm=T)-min(nich.rate,na.rm=T))/3,n),round(2*(max(nich.rate,na.rm=T)-min(nich.rate,na.rm=T))/3,n)),
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
			  legend.text=element_text(face="bold",size=8),
			  legend.title=element_text(face="bold",size=10),
			  legend.position = c(0.77,0.28)) 	
		return(p)
	}
	
	p1=rate.p(biomat.gsu,"rate.bio01.mean",disGSU2,cordacd.gsu,map1.gsu,"MAT \nGSU")
	p2=rate.p(biomat.20km,"rate.bio01.mean",dis20km2,cordacd20km,map1.20km,"MAT \n20km",na.rm=F)
	p3=rate.p(biomat.gsu,"rate.bio12.mean",disGSU2,cordacd.gsu,map1.gsu,"MAP \nGSU")
	p4=rate.p(biomat.20km,"rate.bio12.mean",dis20km2,cordacd20km,map1.20km,"MAP \n20km",na.rm=F)	
	ggarrange(p1,p2,p3,p4,
		nrow=2,ncol =2,widths=c(1,1),heights=c(1,1),labels=c("a","b","c","d"),font.label = list(size = 20))

## compare relationship of lat with NEE
## Liu et al. evaluated the relationship for each species, while here we focus on whether the mean nich evl within each geo units
library(data.table);library(dplyr)
	dis.nee=read.csv("Dataset_S12.csv")
	biomat.nee=get(load("biomat.nee.Rdata"))
	biomatGSU=as.data.frame(fread("Splev.new3.biomat.V2.csv"))[,-1]%>%filter(genus%in%biomat.nee$genus)	
	disGSU=as.data.frame(fread("SpLevDis2.csv")[,c("Adcode99","Species_E1","Lon","Lat")])%>%filter(Species_E1%in%biomat.nee$genus)
	biomat.nee=biomat.nee%>%filter(genus%in%biomatGSU$genus)
	dis.nee=dis.nee%>%filter(Species%in%biomatGSU$genus)
	
	biodis.nee=dis.nee%>%left_join(biomat.nee,by=c("Species"="genus"))%>%filter(age>1&!is.na(rate.bio01.mean))%>%
		select(Species,Longitude,Latitude,starts_with("rate.")&ends_with(".mean"))
	vlist=c(paste("rate",c("bio01","bio05","bio06","bio12","bio16","bio17"),"mean",sep="."))
	biodis.nee[,vlist]=abs(biodis.nee[,vlist])
	#overlay with GSU
	library(sf);library(sp);library(maptools)
	geo=read.csv("Geo-isrm.csv")
	# grid50=readShapeSpatial("area/50km/grid_50km.shp")	
	# villages = st_as_sf(grid50)%>% st_set_crs("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=km")	
	# nee.points = st_as_sf(biodis.nee,coords=c("Longitude","Latitude"))%>% st_set_crs("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=km")
	# map.int = st_intersection(nee.points,villages)%>%as.data.frame()%>%
		# select(Species,ADCODE99,starts_with("rate.")&ends_with(".mean"))%>%distinct()
	# write.csv(map.int,"map.int.csv")
	
	map.int=read.csv("map.int.csv")[,-1]
	nich.pattern.nee = map.int%>%group_by(ADCODE99)%>%
		summarize(MAT=mean(rate.bio01.mean),MAP=mean(rate.bio12.mean),.groups="keep")%>%
		left_join(geo,by="ADCODE99")%>%na.omit()
	cor.mat.nee=round(cor(abs(nich.pattern.nee[,"Lat"]),nich.pattern.nee[,"MAT"]),2)
	cor.map.nee=round(-cor(abs(nich.pattern.nee[,"Lat"]),nich.pattern.nee[,"MAP"]),2)
	#niche based on GSU
	biodisGSU=disGSU%>%left_join(biomatGSU,by=c("Species_E1"="genus"))%>%filter(age>1&!is.na(rate.bio01.mean))%>%
		select(Species_E1,Adcode99,Lon,Lat,starts_with("rate.")&ends_with(".mean"))
	biodisGSU[,vlist]=abs(biodisGSU[,vlist])	
	nich.pattern.gsu=biodisGSU%>%group_by(Lon,Lat)%>%
		summarize(MAT=mean(rate.bio01.mean),MAP=mean(rate.bio12.mean),.groups="keep")
	cor.mat.gsu=round(cor(abs(nich.pattern.gsu[,"Lat"]),nich.pattern.gsu[,"MAT"]),2)
	cor.map.gsu=round(cor(abs(nich.pattern.gsu[,"Lat"]),nich.pattern.gsu[,"MAP"]),2)
	#plot the results
	library(ggpubr)
	library(ggtext)
	theme=theme(axis.text = element_text(size=12,color='black',angle=0,hjust=1),
		axis.title = element_text(size=13,color='black',angle=0),		
		panel.grid.minor.y=element_blank(),
		axis.line.y=element_line(linetype=1,color='black'),
		axis.line.x=element_line(linetype=1,color='black'),
		axis.ticks = element_line(linetype=2,color='black'),
		panel.grid=element_line(linetype=2,color='grey'),
		panel.background = element_blank())
	p1=ggplot(nich.pattern.nee,aes(x=abs(Lat), y=MAT,fill=MAT))+
	labs(x="",y=substitute(paste('MAT Rate in Liu et al. °C ',Myr^-1)))+#,fill="Clade in \nLiu et al. (2022)") +	
	geom_point(size=2,shape=21,color="black",alpha=0.5,show.legend=FALSE)+
	scale_fill_gradient(low="blue",high="red")+
	geom_richtext(x = 20,y = 2.7,fill = NA,label.colour=NA,
                  label = paste("<i><b>r</i></b> = ",cor.mat.nee),size=8)+theme
    p2=ggplot(nich.pattern.gsu,aes(x=abs(Lat), y=MAT,fill=MAT))+
	labs(x="",y=substitute(paste('MAT Rate in this study °C ',Myr^-1)))+#,fill="Clade in \nLiu et al. (2022)") +	
	geom_point(size=2,shape=21,color="black",alpha=0.5,show.legend=FALSE)+
	scale_fill_gradient(low="blue",high="red")+
	geom_richtext(x = 20,y = 3.3,fill = NA,label.colour=NA,
                  label = paste("<i><b>r</i></b> = ",cor.mat.gsu),size=8)+theme
	
	p3=ggplot(nich.pattern.nee,aes(x=abs(Lat), y=MAP,fill=MAP))+
	labs(x="Latitude",y=substitute(paste('MAP Rate in Liu et al. mm ',Myr^-1)))+#,fill="Clade in \nLiu et al. (2022)") +	
	geom_point(size=2,shape=21,color="black",alpha=0.5,show.legend=FALSE)+
	scale_fill_gradient(low="blue",high="red")+
	geom_richtext(x = 20,y = 240,fill = NA,label.colour=NA,
                  label = paste("<i><b>r</i></b> = ",cor.map.nee),size=8)+theme
	p4=ggplot(nich.pattern.gsu,aes(x=abs(Lat), y=MAP,fill=MAP))+
	labs(x="Latitude",y=substitute(paste('MAP Rate in this study mm ',Myr^-1)))+#,fill="Clade in \nLiu et al. (2022)") +	
	geom_point(size=2,shape=21,color="black",alpha=0.5,show.legend=FALSE)+
	scale_fill_gradient(low="blue",high="red")+
	geom_richtext(x = 20,y = 500,fill = NA,label.colour=NA,
                  label = paste("<i><b>r</i></b> = ",cor.map.gsu),size=8)+theme
	ggarrange(p1,p2,p3,p4,
		nrow=2,ncol =2,labels="auto",font.label = list(size = 18))

##compare based on PIC
	library(data.table);library(dplyr)
	dis.nee=read.csv("Dataset_S12.csv")
	biomat.nee=get(load("biomat.nee.Rdata"))
	biomatGSU=as.data.frame(fread("Splev.new3.biomat.V2.csv"))[,-1]%>%filter(genus%in%biomat.nee$genus)	
	disGSU=as.data.frame(fread("SpLevDis2.csv")[,c("Adcode99","Species_E1","Lon","Lat")])%>%filter(Species_E1%in%biomat.nee$genus)
	biomat.nee=biomat.nee%>%filter(genus%in%biomatGSU$genus)
	dis.nee=dis.nee%>%filter(Species%in%biomatGSU$genus)
	vlist=c(paste("rate",c("bio01","bio05","bio06","bio12","bio16","bio17"),"mean",sep="."))
		
	rate.lat.nee=dis.nee%>%group_by(Species)%>%summarize(Lat=mean(Latitude,na.rm=T))
	rate.lat.GSU=disGSU%>%group_by(Species_E1)%>%summarize(Lat=mean(Lat,na.rm=T))
	biodis.GSU=biomatGSU%>%left_join(rate.lat.GSU,by=c("genus"="Species_E1"))%>%filter(age>1)
	biodis.nee=biomat.nee%>%left_join(rate.lat.nee,by=c("genus"="Species"))%>%na.omit()%>%filter(age>1)
	#pic
	library(ape)
	tre=read.tree("FromAo_Smith_100treesV2/1.tre")		
	tre$node.label<-NULL
	#Compute the phylogenetically independent contrasts using the method described by Felsenstein (1985).
	get.pic=function(i,biodis2,tre){
		var.t=biodis2[,i];names(var.t)=biodis2$genus
		tre2=drop.tip(tre, tip=subset(tre$tip.label,!tre$tip.label%in%biodis2$genus))
		pic.X <- pic(var.t, tre2)
		return(pic.X)
	}
	pic.gsu=do.call(cbind,lapply(c("Lat","rate.bio01.mean","rate.bio12.mean"),get.pic,biodis.GSU,tre))%>%as.data.frame()
	pic.nee=do.call(cbind,lapply(c("Lat","rate.bio01.mean","rate.bio12.mean"),get.pic,biodis.nee,tre))%>%as.data.frame()
	colnames(pic.gsu)=colnames(pic.nee)=c("Latitude","MAT","MAP")
	pic.nee$Latitude=abs(pic.nee$Latitude)
	pic.gsu$Latitude=abs(pic.gsu$Latitude)
	round(cor(pic.nee)[1,2],2);cor.test(pic.nee[,1],pic.nee[,2])$p.value
	round(cor(pic.gsu)[1,2],2);cor.test(pic.gsu[,1],pic.gsu[,2])$p.value
	
	round(cor(pic.nee)[1,3],2);cor.test(pic.nee[,1],pic.nee[,3])$p.value
	round(cor(pic.gsu)[1,3],2);cor.test(pic.gsu[,1],pic.gsu[,3])$p.value
	# #overlay with GSU
	# biodis.nee=dis.nee%>%left_join(biomat.nee,by=c("Species"="genus"))%>%filter(age>1&!is.na(rate.bio01.mean))%>%
		# select(Species,Longitude,Latitude,starts_with("rate.")&ends_with(".mean"))
	# biodis.nee[,vlist]=abs(biodis.nee[,vlist])
	
	# library(sf)
	# geo=read.csv("Geo-isrm.csv")
	# grid50=readShapeSpatial("area/50km/grid_50km.shp")	
	# villages = st_as_sf(grid50)%>% st_set_crs("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=km")	
	# nee.points = st_as_sf(biodis.nee,coords=c("Longitude","Latitude"))%>% st_set_crs("+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=km")
	# map.int = st_intersection(nee.points,villages)%>%as.data.frame()%>%
		# select(Species,ADCODE99,starts_with("rate.")&ends_with(".mean"))%>%distinct()
	# map.int=read.csv("map.int.csv")[,-1]
	# nich.pattern.nee = map.int%>%group_by(ADCODE99)%>%
		# summarize(MAT=mean(rate.bio01.mean),MAP=mean(rate.bio12.mean),.groups="keep")%>%
		# left_join(geo,by="ADCODE99")%>%na.omit()
	# write.csv(map.int,"map.int.csv")
	# plot(nich.pattern.nee[,c("Lat","MAT")])
	# nich.pattern.nee.n=nich.pattern.nee%>%filter(Lat>0)
	# nich.pattern.nee.s=nich.pattern.nee%>%filter(Lat<0)
	# cor(nich.pattern.nee.n[,c("Lat","MAT")])
	# cor(nich.pattern.nee.s[,c("Lat","MAT")])
	# cor(abs(nich.pattern.nee[,"Lat"]),nich.pattern.nee[,"MAT"])
	
	# #niche based on GSU
	# biodisGSU=disGSU%>%left_join(biomatGSU,by=c("Species_E1"="genus"))%>%filter(age>1&!is.na(rate.bio01.mean))%>%
		# select(Species_E1,Adcode99,Lon,Lat,starts_with("rate.")&ends_with(".mean"))
	# biodisGSU[,vlist]=abs(biodisGSU[,vlist])	
	# nich.pattern.gsu=biodisGSU%>%group_by(Lon,Lat)%>%
		# summarize(MAT=mean(rate.bio01.mean),MAP=mean(rate.bio12.mean),.groups="keep")
	# plot(nich.pattern.gsu[,c("Lat","MAT")])
	# nich.pattern.gsu.n=nich.pattern.gsu%>%filter(Lat>0)
	# nich.pattern.gsu.s=nich.pattern.gsu%>%filter(Lat<0)
	# cor(nich.pattern.gsu.n[,c("Lat","MAT")])
	# cor(nich.pattern.gsu.s[,c("Lat","MAT")])
	# cor(abs(nich.pattern.gsu[,"Lat"]),nich.pattern.gsu[,"MAT"])
	
	#plot the results
	library(ggpubr)
	theme=theme(axis.text = element_text(size=12,color='black',angle=0,hjust=1),
		axis.title = element_text(size=13,color='black',angle=0),		
		panel.grid.minor.y=element_blank(),
		axis.line.y=element_line(linetype=1,color='black'),
		axis.line.x=element_line(linetype=1,color='black'),
		axis.ticks = element_line(linetype=2,color='black'),
		panel.grid=element_line(linetype=2,color='grey'),
		panel.background = element_blank())
	p1=ggplot(pic.nee,aes(x=Latitude, y=log(MAT),fill=log(MAT)))+
	labs(x="",y="PIC of Niche evolutionary rates \n(MAT) based on Liu et al. (2020)")+#,fill="Clade in \nLiu et al. (2022)") +	
	geom_point(size=2,shape=21,color="black",alpha=0.5,show.legend=FALSE)+
	scale_fill_gradient(low="blue",high="red")+
	geom_smooth(method = "glm",colour="red",fill="#619CFF",size=1)+
	annotate("text", x=20 , y= 2.7,size=6,label=paste("Pearson r = ",round(cor(pic.nee)[1,2],2),
		"\np = ",round(cor.test(pic.nee[,1],pic.nee[,2])$p.value,2)))+theme
	p2=ggplot(nich.pattern.gsu,aes(x=abs(Lat), y=MAT,fill=MAT))+	
	labs(x="Latitude",y=substitute(paste('MAT Niche in this study °C ',Myr^-1,sep=" ")))+#,fill="Clade in \nLiu et al. (2022)") +	
	geom_point(size=2,shape=21,color="black",alpha=0.5,show.legend=FALSE)+
	scale_fill_gradient(low="blue",high="red")+
	#geom_vline(xintercept=0, linetype = "twodash", col="red",size=1)+
	geom_smooth(method = "glm",colour="red",fill="#619CFF",size=1)+
	annotate("text", x=20 , y= 3,size=6,label="Pearson r = 0.55")+theme
	ggarrange(p1,p2,
		nrow=2,ncol =1,labels=c("a","b"),font.label = list(size = 18))
		
