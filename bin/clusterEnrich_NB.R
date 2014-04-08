library(MASS)

#Get parameters
args <- commandArgs(TRUE)
fileAddr <- args[1]
output <- args[2]

#Start to fit model
file<-read.table(fileAddr)

#Determine sample size for parameter estimation
file_row_count<-dim(file)[1]
intercept<-c()
coefficient<-c()
theta <- c()

if(file_row_count > 100000)
{
	sample_count <- min(100000,file_row_count*0.25)
}else
{
	sample_count <- file_row_count
}
	
for(i in c(1:5)) #loop 5 times and get average of parameter
{
	subset<-file[sample(nrow(file),sample_count),]
	read_len = subset$V3 - subset$V2
	sub_data = data.frame(readCount=subset$V5,readLen=read_len)
	
	#start to fit model, use try-catch
	result<-tryCatch(
	m<-glm.nb(readCount~readLen,data=sub_data),
	warning = function(war){
		#for "iternation limit and converge" will be the most likely warning
		return(NULL) 
		#print(paste("Signal is:",signal))
	} 
)#end of tryCatch
	
	if(length(result)!=0){
	intercept <- c(intercept, summary(m)$coef[1,1])
	coefficient <- c(coefficient,summary(m)$coef[2,1])
	theta <- c(theta,m$theta)}
}# end of for
	if(length(intercept)>=3) # at least 3 times converge
	{
		final_intercept	= median(intercept)
		final_coef = median(coefficient)
		final_theta = median(theta)
	}else{#use back-up method
		#print("Start back up method")
		subset<-file[file$V5>2,]
		#print(dim(subset))
		read_len <- subset$V3 - subset$V2
		#print(read_len)	
		sub_data = data.frame(readCount=subset$V5,readLen=read_len)
	
		#result<-tryCatch(
		m<-glm.nb(readCount~readLen,data=sub_data)#,
		#warning = function(war){
		#for "iternation limit and converge" will be the most likely warning
		#	return(NULL) 
		#	} )#end of tryCatch
		final_intercept<-summary(m)$coef[1,1]
		final_coef<-summary(m)$coef[2,1]
		final_theta<-m$theta
	}#end of back-up else

#start to estimate p value
testLen = file$V3-file$V2
test_mu<-final_intercept+final_coef*testLen
test_size<-1/final_theta
pv<-pnbinom(file$V5,size=test_size,mu=test_mu,lower.tail=F,log.p=T)
#pbh<-log(p.adjust(exp(pv),method="BH"))
#qv<-qvalue(pv)
newframe=data.frame(chrom=file$V1,start=file$V2,stop=file$V3,id=file$V4,readsCount=file$V5,strand=file$V6,pvalue=pv)#,p_BH=pbh)#,qvalue=qv)
write.table(newframe,output,sep="\t",quote=F,row.names=F)

