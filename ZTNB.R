#Require at least R3.0 and package VGAM and dependencies
#Programmer: jonghyun.yun@utsouthwestern.edu and beibei.chen@utsouthwestern.edu 
#Usage: Get all the reads count and length combinations of enriched clusters
#Input: BED file with score column as reads number
#Output: table with reads count, length, p value and fdr
#Last modified: 19 Dec 2013

suppressMessages(require(MASS));
suppressMessages(require(VGAM));

args = commandArgs(TRUE);
data = read.table(args[1], sep = "\t");
len = data[,3] - data[,2];
read = data[,5];
cut = as.numeric(as.character(args[2]))
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

# local polynomial regression: read ~ len
# smoothing parameter: 
alp = 0.95
# polynomial degree:
pd = 1
lregfit = loess(tread ~ tlen, weights = twts, span = alp, family ="symmetric", degree = pd);
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

# negative binomail regression with the known predictor log(mu)
nb = vglm(tread_fit ~ 1, posnegbinomial(), weight = twts_fit, maxit = 200, trace = FALSE, offset = logmu,epsilon=10^(-3))
# model parameters
(khat <- Coef(nb)["size"]);
(muhat <- Coef(nb)["munb"]);
nbsize = khat * mu;
nbmu = muhat * mu;

# pvalues
nb.pvalue[cdx] = pnbinom(tread_fit - 1, size = nbsize, mu = nbmu, lower.tail = FALSE) / (1-dnbinom(0, size = nbsize, mu = nbmu))
nb.fdr = p.adjust(nb.pvalue,method='BH')

# output
nbdx = nb.fdr <= cut;
#sum(twts[nbdx]) # the number of clusters whose FDR <= cut
# read counts and cluster length of clusters whose FDR <= cut
nb.out = as.data.frame(matrix(c(tread[nbdx],tlen[nbdx],nb.pvalue[nbdx],nb.fdr[nbdx]),sum(nbdx),4))
names(nb.out) = c("read","length","p","fdr")
outname = paste(args[1],".ztnb",sep="")
write.table(nb.out,outname,sep="\t",quote=F,row.names=F);

