#Require at least R3.0 and package VGAM and dependencies
#Programmer: jonghyun.yun@utsouthwestern.edu and beibei.chen@utsouthwestern.edu 
#Usage: Get all the reads count and length combinations of enriched clusters
#Input: BED file with score column as reads number
#Output: table with reads count, length, p value and fdr
#Last modified: 19 Dec 2013

#command line parameters:
#1:input file name. After running the code, two output files: input.ztnb and input.ztnblog will be created
#2:FDR cutoff
#3:Regression epsilon
#4:Regression step

#suppressMessages(require(MASS));
#suppressMessages(require(VGAM));

require(MASS)
require(VGAM)
args = commandArgs(TRUE);
data = read.table(args[1], sep = "\t");
len = data[,3] - data[,2];
read = data[,5];
cut = as.numeric(as.character(args[2]));
vglm_epsilon = as.numeric(as.character(args[3]));
vglm_step = as.numeric(as.character(args[4]));
# tabularize the data
wts = table(read,len);
temp.coln = as.integer(names(wts[1,]));
temp.rown = as.integer(names(wts[,1]));
drow = length(temp.rown);
dcol = length(temp.coln);

rown = matrix(rep(temp.rown,dcol),drow,dcol);
coln = t(matrix(rep(temp.coln,drow),dcol,drow));

wdx = which(wts>0);
tread = rown[wdx];
tlen = coln[wdx];
twts = wts[wdx];
#print("local polynomial regression")
# local polynomial regression: read ~ len
# smoothing parameter: 
alp = 0.95
# polynomial degree:
pd = 1
lregfit = loess(tread ~ tlen, weights = twts, span = alp, family ="symmetric", degree = pd);
#print("loess finished")
mu = lregfit$fitted;
mu[mu<=0] = min(exp(-4),mu[mu>0]); # bounded mean function
logmu = log(mu);


# compute p-values using N(0,1) for (tread[-cdx] - mu[-cdx]) / sqrt(mu[-cdx])
nb.pvalue=numeric(length(tread));
zscore = (tread-mu)/sqrt(mu); 
cdx = which(zscore < 4);
nb.pvalue[-cdx] = pnorm(zscore[-cdx],lower.tail=FALSE);

tread_fit = tread[cdx];
tlen_fit = tlen[cdx];
twts_fit = twts[cdx];

lregfit = loess(tread_fit ~ tlen_fit, weights = twts_fit, span = alp, family ="symmetric", degree = pd);
mu = lregfit$fitted;
mu[mu<=0] = min(exp(-4),mu[mu>0]); # bounded mean function
logmu = log(mu);

#print("Start vgam binomial");
#print(paste("Epsilon",vglm_epsilon));
#print(paste("Step",vglm_step));
# negative binomail regression with the known predictor log(mu)

intercept1 = c(-1)#seq(-1,0,1);
intercept2 = c(-1)#seq(-1,0,1);
biggest_likelihood = -99999999;
khat = 0;
muhat = 0;
converge_flag = FALSE;
errmsg = paste("Error/warning messages for run using Epsilon:",vglm_epsilon,"and step size:",vglm_step,".");
options(warn=1);

for (i in intercept1)
{
  for (j in intercept2)
  {
    #print(paste("Coef start",i,j));
    #print(paste("Recent likelihood",biggest_likelihood))
    #print(paste("NB parameters",khat,muhat))
    nb<-tryCatch({vglm(tread_fit ~ 1, posnegbinomial(), weight = twts_fit, maxit = 200, trace = FALSE, step = vglm_step,offset = logmu,epsilon=vglm_epsilon,silent=FALSE,coefstart=c(i,j))},warning = function(w){
    #print(paste("typeof warning",typeof(w["message"])))
    warnmsg = strsplit(toString(w["message"])," iteration")[[1]][1];
    if(grepl("convergence not obtained",warnmsg))
    {
      newmsg = paste("Using coefstart c(",i,",",j,"), convergence not obtained in 200 iterarions");
      return(newmsg);
    }
    }, error = function(e){
      newmsg = paste("Using coefstart c(",i,",",j,"), program failed to converge.");
      return(newmsg)
    }, finally = {
    })# end of try-catch
    if (length(nb)>0)
    {
      if(typeof(nb)=="S4")
      {
        #print("Converged, get regression class")
        if(logLik(nb)>=biggest_likelihood & Coef(nb)["size"]<=100)
        {
          biggest_likelihood = logLik(nb);
          khat <- Coef(nb)["size"];
          muhat <- Coef(nb)["munb"];
        }
       
      } else if (typeof(nb)=="character")
      {
        errmsg = paste(errmsg,nb,sep="\n")
      }
    }
  }# end of loop intercept2
 }# end of loop intercept1
#print(paste("Biggest loglikelihood",biggest_likelihood));
#print(paste("NB parameters",khat,muhat))
#print("Finished fitting");
#print(paste("Error message for this run",errmsg))
#Finished fitting, calculate p value
if((biggest_likelihood==-99999999) && (khat==0) && (muhat ==0))
{
  #No model converged,exit
  errmsg = paste("N:No model converged for this run.",errmsg,sep="\n");
} else {
  # model parameters
	errmsg = paste("Y:Coverged for this run.",errmsg,sep="\n");
  nbsize = khat * mu;
  nbmu = muhat * mu;
  nb.pvalue[cdx] = pnbinom(tread_fit - 1, size = nbsize, mu = nbmu, lower.tail = FALSE) / (1-dnbinom(0, size = nbsize, mu = nbmu))
  nb.fdr = p.adjust(nb.pvalue,method='BH')
  # output
  nbdx = nb.fdr <= cut;
  #sum(twts[nbdx]); # the number of clusters whose FDR <= cut
  # read counts and cluster length of clusters whose FDR <= cut
  nb.out = as.data.frame(matrix(c(tread[nbdx],tlen[nbdx],nb.pvalue[nbdx],nb.fdr[nbdx]),sum(nbdx),4))
  names(nb.out) = c("read","length","p","fdr")
  outname = paste(args[1],".ztnb",sep="")
  write.table(nb.out,outname,sep="\t",quote=F,row.names=F);
}#end of output
#Out out error log file
#close(logfile);
write(errmsg,paste(args[1],".ztnblog",sep=""));
