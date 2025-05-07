



rm(list=ls())
setwd()
library(data.table)
set.seed(1)
data=fread("yak_10k.txt",sep="\t",head=T)

rep=30
#Obtain the missing file
for(j in 1:rep){
 print(paste(j," 0.05 & ",rep," times",sep=""))
 X=data
 n=nrow(X)
 m=ncol(X)
 dp=m*n
 uv=runif(dp)
 mr=0.05
 missing=uv<mr
 index.m=matrix(missing,n,m)
 X[index.m]=NA
 index0=X==0
 index1=X==1
 index2=X==2
 indexna=is.na(X)
 fwrite(X,file=paste("10k_0.05_NA",j,".txt",sep=""),quote=F,sep="\t",col.name=T,row.name=F)
}
#Obtain the missing cluster file
for(i in 1:rep){
 print(paste(i," & 0.05 ",rep," times",sep=""))
 X=data
 n=nrow(X)
 m=ncol(X)
 dp=m*n
 uv=runif(dp)
 mr=0.05
 missing=uv<mr
 index.m=matrix(missing,n,m)
 X[index.m]=NA
 index0=X==0
 index1=X==1
 index2=X==2
 indexna=is.na(X)
 X[index0]="0"
 X[index1]="1"
 X[index2]="2"
 X[indexna]="1"
 fwrite(X,file=paste("10k_0.05_cluster",i,".txt",sep=""),quote=F,sep="\t",col.name=T,row.name=F)
}

#Imputation by Beagle 5.4
n=30
all.m=NULL
all.c=NULL
system.time({
all.l=NULL
for(i in 1:n){
print(paste(" 0.05 & ",i," times",sep=""))#改缺失率
#i=1
 data1=read.table(paste("10k_0.05_NA",i,".txt",sep=""),
                       sep="\t",
                       header=T)#读取含有缺失值的10K文件,缺失率和文件路径
 data2=fread("yak_10k.txt",sep="\t")#读取10K原始数据
 raw_path=getwd()
 dir.create("mid")
 setwd("mid")
 index0=data1==0
 index1=data1==1
 index2=data1==2
 indexna=is.na(data1)
 X=data1
 X[index0]="AA"
 X[index1]="AG"
 X[index2]="GG"
 X[indexna]="NN"
 hmp=t(X)
 N=nrow(hmp)
 rs=seq(1,N)
 alleles=rep("A/G",N)
 chrom=rep(2,N)
 pos=seq(1,N)
 hmp5.11=matrix(NA,N,7)
 out.hmp=cbind(rs,alleles,chrom,pos,hmp5.11,hmp)
 colnames(out.hmp)[5:11]=c("strand","assembly","center","protLSID","assayLSID","panelLSID","QCcode")
 colnames(out.hmp)[12:ncol(out.hmp)]=paste0("V",1:(ncol(out.hmp)-11))
 write.table(out.hmp,"mid1.hmp",col.names=T,sep="\t",row.names=F,quote=F)
 rm(hmp,out.hmp,X)
 system(paste("perl",
  "Tassel/run_pipeline.pl",
  "-Xms512m -Xmx2g",
  paste("-fork1 -h",
   hapmap_file="mid1.hmp",
   "-export -exportType VCF -runfork1"),
   sep=" "))
 system("java -Xmx50g -jar beagle.22Jul22.46e.jar gt=mid1.vcf out=mid3")
 lj=gzfile("mid3.vcf.gz","rt")
 content=readLines(lj)
 close(lj)
 output_file="mid3.vcf"
 writeLines(content, output_file)
 vcf_path = file.path(getwd(),"mid3.vcf")
 snp_path = file.path(getwd(),"mid")
 command = paste("perl",
     "vcf2hpm2_wy_v1.pl", 
     shQuote(vcf_path), 
     shQuote(snp_path), 
     sep = " ")		 
 system(command)
 X4=data.table::fread("mid2.txt",sep=" ")
 X5=as.matrix(X4[,-c(1:12)])
 X6=t(X5)
 index.AA=X6=="AA"
 index.AG=X6=="AG"
 index.GA=X6=="GA"
 index.GG=X6=="GG"
 nr=nrow(X6)
 nc=ncol(X6)
 genotype.n=matrix(0,nr,nc)
 genotype.n[index.AA]=0
 genotype.n[index.AG]=1
 genotype.n[index.GA]=1
 genotype.n[index.GG]=2
 X.bgl=genotype.n
 colnames(X.bgl)=colnames(data2)
 setwd(raw_path)
 finished_file_path = paste0(raw_path,"/Beagle_wc")
 dir.create(finished_file_path,showWarnings = FALSE)
 setwd(finished_file_path)
 fwrite(X.bgl,file=paste("10k_0.05_Beagle_wc",i,".txt",sep=""),quote=F,sep="\t",col.name=T,row.name=F)
 #文件名称
 setwd(raw_path)
 index.NA=is.na(data1)
 real.genotype=as.data.frame(data2)[index.NA]
 imputed.genotype=as.data.frame(X.bgl)[index.NA]
 sum_data=sum(real.genotype==imputed.genotype)
 matching=sum_data/length(real.genotype)
 cor1=cor(real.genotype,imputed.genotype)
 all.m=append(all.m,matching)
 all.c=append(all.c,cor1)
 unlink("mid",recursive=T)
}
})
write.csv(all.m,"10k_0.05_matching_Beagle.csv",quote=FALSE)#改文件名称
write.csv(all.c,"10k_0.05_cor_Beagle.csv",quote=FALSE)

#Imputation by KBeagle
Impute.subgroup=function(XYcluster){
     print("subgroup will read file...")
     X=read.table(XYcluster,header = TRUE,na.strings="NA",sep="\t")
     #X=read.table("/home/student/workshop/GXY/YAK/10K/kbeagle/mid_1/cluster_1.txt",header = TRUE,na.strings="NA",sep="\t")
     print("subgroup has done in reading file.")
     X[X==0]="AA"
     X[X==1]="AG"
     X[X==2]="GG"
     X[is.na(X)]="NN"
	 print("subgroup is preparing to convert the file Format.")
     hmp=t(X)
     N=nrow(hmp)
     rs=seq(1,N)
     alleles=rep("A/G",N)
     chrom=rep(2,N)
     pos=seq(1,N)
     hmp5.11=matrix(NA,N,7)
     out.hmp=cbind(rs,alleles,chrom,pos,hmp5.11,hmp)
     colnames(out.hmp)[5:11]=c("strand","assembly","center","protLSID","assayLSID","panelLSID","QCcode")
     colnames(out.hmp)[12:ncol(out.hmp)]=paste0("V",1:(ncol(out.hmp)-11))
     write.table(out.hmp,"mid1.hmp",col.names=T,sep="\t",row.names=F,quote=F)
     rm(hmp,out.hmp,X)
	 print("subgroup has converted from num format to hapmap format.")
	 tassel_script = file.path(getwd(),"run_pipeline.pl")
     system(paste("perl",
         tassel_script,
         "-Xms512m -Xmx2g",
         "-fork1 -h",
         hapmap_file="mid1.hmp",
         "-export -exportType VCF -runfork1",
         sep=" "))
     print("subgroup has converted from hapmap format to VCF format.")
     system("java -Xmx50g -jar beagle.22Jul22.46e.jar gt=mid1.vcf out=mid3")
	 print("The file has been completed imputation.")
	 tryCatch({
         lj = gzfile("mid3.vcf.gz", "rt")
    }, error = function(e) {
         print("Error occurred while opening the file.")
         system("java -Xmx20g -jar beagle.22Jul22.46e.jar gt=mid1.vcf out=mid3")
         lj = gzfile("mid3.vcf.gz", "rt")
    })
     content = readLines(lj)  
     close(lj)
     output_file = "mid3.vcf"
     writeLines(content, output_file)
     vcf_path = file.path(getwd(),"mid3.vcf")
     snp_path = file.path(getwd(),"mid")
     command = paste("perl",
         "vcf2hpm2_wy_v1.pl", 
         shQuote(vcf_path), 
         shQuote(snp_path), 
         sep = " ")		 
     system(command)
	 print("subgroup has converted from VCF format to hapmap format.")
     X4=data.table::fread("mid2.txt",sep=" ")
     X5=as.matrix(X4[,-c(1:12)]) 
     X6=t(X5)
     index.AA=X6=="AA"
     index.AG=X6=="AG"
	 index.GA=X6=="GA"
     index.GG=X6=="GG"
     nr=nrow(X6)
     nc=ncol(X6)
     genotype.n=matrix(0,nr,nc)
     genotype.n[index.AA]=0
     genotype.n[index.AG]=1
	 genotype.n[index.GA]=1
     genotype.n[index.GG]=2
	 XYcluster=gsub(".txt",".kbeagle",XYcluster)
	 print(dim(genotype.n))
     write.table(genotype.n,file=XYcluster,quote=F,sep="\t",col.names=F,row.names=F)
     print("subgroup has finished.")
}
library(data.table)
library(factoextra)
library(cluster)
library(ggplot2)
library(parallel)
library(foreach)
library(doParallel)
data3=data.table::fread("yak_10k.txt",sep="\t",header=T)
nr=nrow(data3)
b1=data.frame(IID=1:nr)
n=30
set.seed(1)
all.m=NULL
all.c=NULL
system.time({
for(k in 1:n){
	 data=read.table(paste("10k_0.05_cluster",k,".txt",sep=""),
                          sep="\t",
                          header=T)
     data1=read.table(paste("10k_0.05_NA",k,".txt",sep=""),
                          sep="\t",
                          header=T)     
	 
     sil_width=numeric(10)
     for(kn in 2:10){
      km=kmeans(data,centers=kn,nstart=25)
      ss=silhouette(km$cluster,dist(data))
      sil_width[kn]=mean(ss[,3])
    }
     KN=which.max(sil_width) 

     data_IID=data.frame(IID=b1,data)
     data1_IID=data.frame(IID=b1,data1)
	 data_kmeans=kmeans(data,centers=KN,nstart=25)
     t=table(data_IID$IID,data_kmeans$cluster)
     colnames(t)=paste("type",1:KN,sep="")
	 t1=cbind(b1,apply(t,2,as.numeric))
	 
     for(i in 1:KN){
         d=t1[grep("1",t1[,i+1]),]
         type=merge(d,data1_IID,all=FALSE,by="IID")
         type_result=type[,-c(2:(KN+1))]
         X=type_result[,-1]
         fwrite(X,file=paste("cluster_",i,".txt",sep=""),quote=F,sep="\t",col.name=T,row.name=F)
    }
	 
     for(n in 1:KN){
	     dir_name = paste0("mid_",n)
		 if(!file.exists(dir_name)){
             dir.create(dir_name)
        }
		 raw_cluster_path = paste(getwd(),"/cluster_",n,".txt",sep="")
		 now_cluster_path = paste0(getwd(),"/mid_",n)
		 file.copy(from = raw_cluster_path,to = now_cluster_path)
		 raw_tassel_1 = paste0(getwd(),"/sTASSEL.jar")
		 file.copy(from = raw_tassel_1,to = now_cluster_path)
		 raw_tassel_2 = paste0(getwd(),"/run_pipeline.pl")
		 file.copy(from = raw_tassel_2,to = now_cluster_path)
		 raw_tassel_3 = paste0(getwd(),"/lib")
		 file.copy(from=raw_tassel_3,to=now_cluster_path,recursive = TRUE)
		 raw_perl = paste0(getwd(),"/vcf2hpm2_wy_v1.pl")
		 file.copy(from = raw_perl,to = now_cluster_path)
		 raw_beagle = paste0(getwd(),"/beagle.22Jul22.46e.jar")
		 file.copy(from = raw_beagle,to = now_cluster_path)
	}
	 
	 num_cores=detectCores()
     cl=makeCluster(num_cores)
	 registerDoParallel(cl)
     raw_path=getwd()
     XYparallel=foreach(j = 1:KN) %dopar% {
       parallel_path=paste0(raw_path,"/mid_",j)
       dir.create(parallel_path,showWarnings = FALSE)
       setwd(parallel_path)
       data_cluster=file.path(getwd(), paste0("cluster_", j, ".txt"))
       result=Impute.subgroup(XYcluster = data_cluster)
    }
	 stopCluster(cl)
	
     setwd(raw_path)
     for(n in 1:KN){
	   raw_cluster_path=paste0(getwd(),"/mid_",n,"/cluster_",n,".kbeagle")
	   now_cluster_path=getwd()
	   file.copy(from=raw_cluster_path,to = now_cluster_path)
    }
	 	
     X.bgl2 = NULL
     for(i in 1:KN){
         cluster_some=data.table::fread(paste("cluster_",i,".kbeagle",sep=""),sep="\t")#删除head=T
         d=t1[grep("1",t1[,i+1]),]
         type=merge(d,data1_IID,all=FALSE,by="IID")
         type_result=type[,-c(2:(KN+1))]
         X.bgl=cbind(type_result[,1],cluster_some)
         X.bgl2=rbind(X.bgl2,X.bgl)
    }

     colnames(X.bgl2)=as.character(c("IID",colnames(data1)))
     X2=data.frame(X.bgl2)
	 loc=match(b1$IID,X2$IID)
	 W=X2[loc,]
     W1=W[,-1]
	 finished_file_path = paste0(raw_path,"/KBeagle_wc")
     dir.create(finished_file_path,showWarnings = FALSE)
     setwd(finished_file_path)
	 fwrite(W1,file=paste("10k_0.05_KBeagle_wc",k,".txt",sep=""),quote=F,sep="\t",col.name=T,row.name=F)
	 setwd(raw_path)
	 for(n in 1:KN){
	     unlink(paste0(raw_path,"/mid_",n),recursive = TRUE)
		 file.remove(paste0(raw_path,"/cluster_",n,".txt"))
		 file.remove(paste0(raw_path,"/cluster_",n,".kbeagle"))
	}
	 index.NA=is.na(data1)
	 real.genotype=as.data.frame(data3)[index.NA]
     imputed.genotype=as.data.frame(W1)[index.NA]
	 cor1=cor(real.genotype,imputed.genotype)
     sum_data=sum(real.genotype==imputed.genotype)
     matching=sum_data/length(real.genotype)
     all.m=append(all.m,matching)
	 all.c=append(all.c,cor1)
}
})
write.csv(all.m,"10k_0.05_matching_KBeagle.csv",quote=FALSE)
write.csv(all.c,"10k_0.05_cor_KBeagle.csv",quote=FALSE)

#GP
G2P=function(data_na,data,NQTN,h2){
 colnames(data_na)=NULL
 index_NA=apply(data_na,2,sum) 
 NA_col=which(is.na(index_NA))
 variance=apply(data_na,2,var,na.rm=TRUE)
 variable_SNPs=which(variance != 0)           
 filtered_cols=intersect(NA_col,variable_SNPs)
 QTN.position=sample(filtered_cols, NQTN, replace = FALSE) 
 SNPQ=as.matrix(data[,QTN.position,with=FALSE])
 addeffect=rnorm(NQTN,0,1)
 effect=SNPQ%*%addeffect
 effectvar=var(effect)
 residualvar=(effectvar-h2*effectvar)/h2
 nr=nrow(data)
 residual=rnorm(nr,0,sqrt(residualvar))
 y=effect+residual
 return(list(y=y,QTN.position=QTN.position))
}
n=30
mean_all1=NULL
mean_all2=NULL
mean_all3=NULL
for(i in 1:n){
     #i=1
	 data1=fread("yak_10k.txt")
 	 data2=read.table(paste("10k_0.05_Beagle_wc",i,".txt",sep=""),
                       sep="\t",
					   header=T)
 	 colnames(data2)=gsub("^X(?=\\d)", "", colnames(data2), perl = TRUE)					   
 	 data3=read.table(paste("10k_0.05_KBeagle_wc",i,".txt",sep=""),
                       sep="\t",
					   header=T)
	 colnames(data3)=gsub("^X(?=\\d)", "", colnames(data3), perl = TRUE)
 	 data_na=read.table(paste("10k_0.05_NA",i,".txt",sep=""),
                       sep="\t",
					   header=T) #na文件
 	 colnames(data_na)=gsub("^X(?=\\d)", "", colnames(data_na), perl = TRUE)
     nr_ID=nrow(data_na)
	 ID_col=seq(1,nr_ID)
	 
	 split_names=strsplit(colnames(data1), "__")
 	 first_part=sapply(split_names, function(x) as.numeric(x[1]))
 	 second_part=sapply(split_names, function(x) as.numeric(x[2]))
 	 sorted_index=order(first_part, second_part)
 	 data1=data1[,..sorted_index]
 	 data2=data2[,sorted_index]
	 data3=data3[,sorted_index]
 	 data_na=data_na[,sorted_index]

	 NQTN=20
	 h2=0.75
	 data=data1
	 phe=G2P(data_na,data,NQTN,h2)
	 y=phe$y
	 myY=cbind(taxa=ID_col,trait=y)
	 colnames(myY)=c("taxa","trait")
	 myGD1=cbind(taxa=ID_col,data1)
	 myGD2=cbind(taxa=ID_col,data2)
	 myGD3=cbind(taxa=ID_col,data3)

	 #生成myGM
	 SNP_col=as.character(colnames(data_na))
	 raw_data=fread("80yaks.hmp.txt",header=F)
	 raw_data$V1=trimws(as.character(raw_data$V1))
   	 SNP_col=trimws(SNP_col)
 	 matching_rows=raw_data[raw_data$V1 %in% SNP_col,]
 	 SNP_information=matching_rows[,c(1, 3, 4)]
	 colnames(SNP_information)=c("SNP","Chromosome","Position")
	 myGM=as.data.frame(SNP_information)
 
	 cor_all1=NULL
     cor_all2=NULL
     cor_all3=NULL
	 nfold=5
	 nr2=nrow(myY)
	 sets=sample(cut(1:nr2,nfold,labels=FALSE),nr2)
	 for(j in 1:nfold){
		 testing=myY 
	     testing[sets==j,2]=NA 
		 V1=colnames(testing)[2] 
		 training=testing[!is.na(testing[,2]),]
		 training=data.frame(training)
		 testing_file=testing[is.na(testing[,2]),]
		 real_indice=match(testing_file[,"taxa"],myY[,"taxa"])
		 real_value=myY[real_indice,"trait"]
		   
		 myGAPIT1=GAPIT(Y=training,
             GD=myGD1,
             GM=myGM,
			 Random.model=FALSE,
             model=c("gblup"),
			 Geno.View.output=FALSE)
		 Pred1=myGAPIT1$Pred
		 pre_indice1=match(testing_file[,"taxa"],Pred1$Taxa)
		 pre_value1=Pred1$Prediction[pre_indice1]
		 cor1=cor(pre_value1,real_value)
         cor_all1=cbind(cor1,cor_all1)		 
         
		 myGAPIT2=GAPIT(Y=training,
             GD=myGD2,
             GM=myGM,
			 Random.model=FALSE,
             model=c("gblup"),
			 buspred=TRUE,
			 lmpred=TRUE,
			 Geno.View.output=FALSE)
		 Pred2=myGAPIT2$Pred
		 pre_indice2=match(testing_file[,"taxa"],Pred2$Taxa)
		 pre_value2=Pred2$Prediction[pre_indice1]
		 cor2=cor(pre_value2,real_value)
		 cor_all2=cbind(cor2,cor_all2)	

		 myGAPIT3=GAPIT(Y=training,
             GD=myGD3,
             GM=myGM,
			 Random.model=FALSE,
             model=c("gblup"),
			 buspred=TRUE,
			 lmpred=TRUE,
			 Geno.View.output=FALSE)
		 Pred3=myGAPIT3$Pred
		 pre_indice3=match(testing_file[,"taxa"],Pred3$Taxa)
		 pre_value3=Pred3$Prediction[pre_indice3]
		 cor3=cor(pre_value3,real_value)
		 cor_all3=cbind(cor3,cor_all3)		 
        }
	 mean_all1=cbind(cor_all1,mean_all1)
	 mean_all2=cbind(cor_all2,mean_all2)
	 mean_all3=cbind(cor_all3,mean_all3)
}
write.table(mean_all1,"10_acc_real.csv",sep=",",quote=F,row.names=F,col.names=F)
write.table(mean_all2,"10k_acc_beagle.csv",sep=",",quote=F,row.names=F,col.names=F)
write.table(mean_all3,"10k_acc_kbeagle.csv",sep=",",quote=F,row.names=F,col.names=F)
h2_all1=NULL
h2_all2=NULL
h2_all3=NULL
for(i in 1:n){
 data1=fread("yak_10k.txt")
 data2=read.table(paste("10k_0.05_Beagle_wc",i,".txt",sep=""),
                       sep="\t",
					   header=T)
 colnames(data2)=gsub("^X(?=\\d)", "", colnames(data2), perl = TRUE)					   
 data3=read.table(paste("10k_0.05_KBeagle_wc",i,".txt",sep=""),
                       sep="\t",
					   header=T)
 colnames(data3)=gsub("^X(?=\\d)", "", colnames(data3), perl = TRUE)
 data_na=read.table(paste("10k_0.05_NA",i,".txt",sep=""),
                       sep="\t",
					   header=T)
 colnames(data_na)=gsub("^X(?=\\d)", "", colnames(data_na), perl = TRUE)
 	 
 split_names=strsplit(colnames(data1), "__")
 first_part=sapply(split_names, function(x) as.numeric(x[1]))
 second_part=sapply(split_names, function(x) as.numeric(x[2]))
 sorted_index=order(first_part, second_part)
 data1=data1[,..sorted_index]
 data2=data2[,sorted_index]
 data3=data3[,sorted_index]
 data_na=data_na[,sorted_index]

 NQTN=20
 h2=0.75
 data=data1
 phe=G2P(data_na,data,NQTN,h2)
 y=phe$y
  
 Vp=var(y)
 G1=A.mat(data1) 
 gblup_result1=mixed.solve(y=y,K=G1)
 GEBV1=gblup_result1$u
 Vg1=var(GEBV1)
 h2_1=Vg1/Vp
 h2_all1=rbind(h2_all1,h2_1)

 G2=A.mat(data2) 
 gblup_result2=mixed.solve(y=y,K=G2)
 GEBV2=gblup_result2$u
 Vg2=var(GEBV2)
 h2_2=Vg2/Vp
 h2_all2=rbind(h2_all2,h2_2)

 G3=A.mat(data3) 
 gblup_result3=mixed.solve(y=y,K=G3)
 GEBV3=gblup_result3$u
 Vg3=var(GEBV3)
 h2_3=Vg3/Vp
 h2_all3=rbind(h2_all3,h2_3)
}
write.table(h2_all1,"10K_real_h2.csv",sep=",",quote=F,row.names=F,col.names=T)
write.table(h2_all2,"10K_beagle_h2.csv",sep=",",quote=F,row.names=F,col.names=T)
write.table(h2_all3,"10K_kbeagle_h2.csv",sep=",",quote=F,row.names=F,col.names=T)

#GWAS
library(ggplot2)
G2P=function(data_na,data,NQTN,h2){
 colnames(data_na)=NULL
 index_NA=apply(data_na,2,sum) 
 NA_col=which(is.na(index_NA))
 variance=apply(data_na,2,var,na.rm=TRUE)
 variable_SNPs=which(variance != 0)           
 filtered_cols=intersect(NA_col,variable_SNPs)
 QTN.position=sample(filtered_cols, NQTN, replace = FALSE) 
 SNPQ=as.matrix(data[,QTN.position,with=FALSE])
 addeffect=rnorm(NQTN,0,1)
 effect=SNPQ%*%addeffect
 effectvar=var(effect)
 residualvar=(effectvar-h2*effectvar)/h2
 nr=nrow(data)
 residual=rnorm(nr,0,sqrt(residualvar))
 y=effect+residual
 return(list(y=y,QTN.position=QTN.position))
}
n=30
FDR1_30=NULL
FDR2_30=NULL
FDR3_30=NULL
for(i in 1:n){
 set.seed(1)
 data1=fread("yak_10k.txt")
 data2=read.table(paste("10k_0.05_Beagle_wc",i,".txt",sep=""),
                       sep="\t",
					   header=T)
 colnames(data2)=gsub("^X(?=\\d)", "", colnames(data2), perl = TRUE)					   
 data3=read.table(paste("10k_0.05_KBeagle_wc",i,".txt",sep=""),
                       sep="\t",
					   header=T)
 colnames(data3)=gsub("^X(?=\\d)", "", colnames(data3), perl = TRUE)
 data_na=read.table(paste("10k_0.05_NA",i,".txt",sep=""),
                       sep="\t",
					   header=T)
 colnames(data_na)=gsub("^X(?=\\d)", "", colnames(data_na), perl = TRUE)
 
 split_names=strsplit(colnames(data1), "__")
 first_part=sapply(split_names, function(x) as.numeric(x[1]))
 second_part=sapply(split_names, function(x) as.numeric(x[2]))
 sorted_index=order(first_part, second_part)
 data1=data1[,..sorted_index]
 data2=data2[,sorted_index]
 data3=data3[,sorted_index]
 data_na=data_na[,sorted_index]

 QTN=NULL
 FDR1_all=NULL
 FDR2_all=NULL
 FDR3_all=NULL
 
 SNP_col=as.character(colnames(data_na))
 raw_data=fread("80yaks.hmp.txt",header=F)
 raw_data$V1=trimws(as.character(raw_data$V1))
 SNP_col=trimws(SNP_col)
 matching_rows=raw_data[raw_data$V1 %in% SNP_col, ]
 SNP_information=matching_rows[, c(1, 3, 4)]
 colnames(SNP_information)=c("SNP","Chromosome","Position")
 myGM=as.data.frame(SNP_information)

 for(j in 1:20){
  NQTN=20
  h2=0.75
  data=data1
  phe=G2P(data_na,data,NQTN,h2)
  y=phe$y
  myseqQTN=phe$QTN.position
  QTN.position=myseqQTN
  QTN=rbind(myseqQTN,QTN)
  taxa=seq(1,nrow(data1))
  myY=data.frame(taxa,trait=y)

  myGD1=data.frame(taxa,data1)
  myGAPIT1=GAPIT(Y=myY,GD=myGD1,GM=myGM,PCA.total=0,Geno.View.output=F,model="Blink",Multiple_analysis=F,QTN.position=myseqQTN)
  myGWAS1=myGAPIT1$GWAS[,c(1:4)]
  myStat1=GAPIT.FDR.TypeI(WS=1,GM=SNP_information,seqQTN=myseqQTN,GWAS=myGWAS1,maxOut=nrow(myGWAS1))
  FDR1=myStat1$FDR
  FDR1_all=cbind(FDR1_all,FDR1)
 
  myGD2=data.frame(taxa,data2)
  myGAPIT2=GAPIT(Y=myY,GD=myGD2,GM=myGM,PCA.total=0,Geno.View.output=F,model="Blink",Multiple_analysis=F,QTN.position=myseqQTN)
  myGWAS2=myGAPIT2$GWAS[,c(1:4)]
  myStat2=GAPIT.FDR.TypeI(WS=1,GM=SNP_information,seqQTN=myseqQTN,GWAS=myGWAS2,maxOut=nrow(myGWAS2))
  FDR2=myStat2$FDR
  FDR2_all=cbind(FDR2_all,FDR2)

  myGD3=data.frame(taxa,data3)
  myGAPIT3=GAPIT(Y=myY,GD=myGD3,GM=myGM,PCA.total=0,Geno.View.output=F,model="Blink",Multiple_analysis=F,QTN.position=myseqQTN)
  myGWAS3=myGAPIT3$GWAS[,c(1:4)]
  myStat3=GAPIT.FDR.TypeI(WS=1,GM=SNP_information,seqQTN=myseqQTN,GWAS=myGWAS3,maxOut=nrow(myGWAS3))
  FDR3=myStat3$FDR
  FDR3_all=cbind(FDR3_all,FDR3)
 }
 FDR1_30=cbind(FDR1_all,FDR1_30)
 FDR2_30=cbind(FDR2_all,FDR2_30)
 FDR3_30=cbind(FDR3_all,FDR3_30)
 }
power1=myStat1$Power
FDR1_30=apply(FDR1_30,1,mean)
plot_data1=data.frame(power=power1,FDR=FDR1_30,group="non-missing")
write.table(plot_data1,"10k_nonmissing_gwas.csv",sep=",",quote=F,row.names=F,col.names=T)
power2=myStat2$Power
FDR2_30=apply(FDR2_30,1,mean)
plot_data2=data.frame(power=power2,FDR=FDR2_30,group="Beagle 5.4")
write.table(plot_data2,"10k_beagle_gwas.csv",sep=",",quote=F,row.names=F,col.names=T)
power3=myStat3$Power
FDR3_30=apply(FDR3_30,1,mean)
plot_data3=data.frame(power=power3,FDR=FDR3_30,group="KBeagle")
write.table(plot_data3,"10k_kbeagle_gwas.csv",sep=",",quote=F,row.names=F,col.names=T)







