library(MASS)
#Get parameters
args <- commandArgs(TRUE)
fileAddr <- args[1]
output <- args[2]

#Start to fit model
file<-read.table(fileAddr)

#Determine sample size for parameter estimation
intercept<-c()
coefficient<-c()
theta <- c()
#print(dim(file))	
#subset<-file[file[,"reads_in_cluster"]<=quantile(file$reads_in_cluster,0.99999)]
read_len = file$V3 - file$V2
sub_data = data.frame(readCount=file$V5,readLen=read_len)

#print(dim(sub_data))	
#start to fit model, use try-catch
result<-tryCatch(
	m<-glm(readCount~1,data=sub_data,family=poisson,offset=log(readLen)),
	warning = function(war){
		#for "iternation limit and converge" will be the most likely warning
		return(NULL) 
		#print(paste("Signal is:",signal))
	} 
)#end of tryCatch
	
if(length(result)!=0){
	intercept <- summary(m)$coef[1]
	testLen = file$V3-file$V2
	#test_mu<-exp(final_intercept+final_coef*testLen)#jonghyun
	test_mu<-exp(intercept) #Andy
	pv<-ppois(file$V5,lambda=test_mu,lower.tail=F,log.p=T)
	#pbh<-log(p.adjust(exp(pv),method="BH"))

	}else{#use back-up method
		#print("Start back up method")
		subset<-file[file$V5>2,]
		#print(dim(subset))
		read_len <- subset$V3 - subset$V2
		#print(read_len)	
		sub_data = data.frame(readCount=subset$V5,readLen=read_len)
	
		#result<-tryCatch(
		m<-glm.nb(readCount~readLen,data=sub_data)#,
		#for "iternation limit and converge" will be the most likely warning
		#	return(NULL) 
		#	} )#end of tryCatch
		final_intercept<-summary(m)$coef[1,1]
		final_coef<-summary(m)$coef[2,1]
		final_theta<-m$theta

		#start to estimate p value
		testLen = file$V3-file$V2
		test_mu<-final_intercept+final_coef*testLen
		test_size<-1/final_theta #may need to mofidy accroding to JongHyun's formular
		pv<-pnbinom(file$V5,size=test_size,mu=test_mu,lower.tail=F,log.p=T)
		#pbh<-log(p.adjust(exp(pv),method="BH"))
#qv<-qvalue(pv)
} # end of back-up else

newframe=data.frame(chrom=file$V1,start=file$V2,stop=file$V3,id=file$V4,score=file$V5,strand=file$V6,pvalue=pv)#,p_BH=pbh)#,qvalue=qv)
write.table(newframe,output,sep="\t",quote=F,row.names=F)
