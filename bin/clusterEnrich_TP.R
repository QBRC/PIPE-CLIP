suppressMessages(library(VGAM))
suppressMessages(library(foreign))
suppressMessages(library(boot))

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

	
subset<-file[file$V5<=quantile(file$V5,0.9),]
read_len = subset$V3 - subset$V2
sub_data = data.frame(readCount=subset$V5,readLen=read_len)
#start to fit model, use try-catch
result<-tryCatch(
	m<-vglm(readCount~readLen,data=sub_data,family=pospoisson()),
	warning = function(war){
		#for "iternation limit and converge" will be the most likely warning
		return(NULL) 
		#print(paste("Signal is:",signal))
	} 
)#end of tryCatch
	
if(length(result)!=0){
	#m.lambda = exp(m@predictors)
	testLen = file$V3-file$V2
	test_mu<-exp(coefficients(m)[1]+coefficients(m)[2]*testLen)#jonghyun
	#test_mu<-coefficients(m)[1]+coefficients(m)[2]*testLen #Andy
	#pv<-ppois(file$V5,lambda=m.lambda,lower.tail=F)/(1-dpois(0,m.lambda))#jonghyun
	#pv<-ppois(file$V5,test_mu,lower.tail=F,log.p=T)#Andy
	
	pv<-ppois(file$V5,test_mu,lower.tail=F)/(1-dpois(0,test_mu))#jonghyun
	pv <- log(pv)	#jonghyun
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
