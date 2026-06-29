args <- commandArgs(trailingOnly=TRUE)
m.exp <- as.numeric(args[1])
skew <- as.numeric(args[2])
kurt <- as.numeric(args[3])


source('utils.R')
ns <- 10:25
mu.star <- mu.null <- 0
rsas <- rsas.init(skew=skew,kurtosis=kurt)
imbalance <- 1.5
by.n <- sapply(ns, function(n) {
    cat('.')
    m <- round(n^m.exp)
    stats <- replicate(n^4, {
        test.stats <- sapply(
            list(logodds=function(...)dgp.logodds(...,tau2.null=.25),
                 smd=function(...)dgp.smd(...,tau2.null=1,yi.distr='exp')),
            function(dgp) {
                yv <- dgp(n=n,m=m,imbalance=imbalance,re.distr=rsas)
                y <- yv$y; v <- yv$v; se2 <- yv$se2
                ma.LR(y=y,sigma2=se2,mu.null=mu.null)
            })
        c(logodds=test.stats[,'logodds'],smd=test.stats[,'smd'])
    })
    ecdfs <- apply(stats,1,function(x)ecdf(pchisq(x,df=1,lower=FALSE)))
    sapply(ecdfs,function(ecdf)max(abs(ecdf(knots(ecdf)) - knots(ecdf))))
})


filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(m.exp,skew=skew,kurt=kurt,by.n=by.n,ns=ns,file=filename)

