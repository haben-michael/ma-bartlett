
## Intro section--Bartlett correction in the Gaussian model 
save.path <- ''
png.file <- ''
filelist <- dir(save.path)
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
by.sim <- lapply(filelist, function(file) {
    load(paste0(save.path,file))
    by.n
})
by.sim <- simplify2array(by.sim)
load(paste0(save.path,filelist[1]))
by.n <- apply(by.sim,1:2,median)
by.n <- by.n[-grep('simple',rownames(by.n)),] 
png(png.file)
matplot(ns,t(by.n),log='xy',type='l',col=1,lty=1:3,xlab='n (log scale)',ylab='K-S distance (log scale)')
rate.lms <- apply(by.n,1,function(ks)lm(log10(ks)~log10(ns)))
rates <- sapply(rate.lms,function(rate.lm)coef(rate.lm)[[2]])
for(i in 1:length(rate.lms))abline(rate.lms[[i]],lty=i)
legend('topright',lty=1:4,legend=paste0(names(rates),' (',round(rates,2),')'))
print(length(filelist))


## Simulation section--effect of primary study size
data.path <- ''
out.path <- ''
filelist <- dir(data.path)
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
print(length(filelist))
by.sim <- lapply(filelist, function(file) {
    load(paste0(data.path,file))
    list(by.n=by.n,ns=ns,m.exp=m.exp)
})
long <- lapply(by.sim, function(sim) {
    by.n <- sim$by.n
    ns <- sim$ns
    rates <- apply(by.n,1,function(ks)coef(lm(log(ks)~log(ns)))[2])
    settings <- c(m.exp=sim$m.exp)
    data.frame(estimator=names(rates),rate=unname(rates),t(settings))
})
long <- Reduce(rbind,long)
long <- aggregate(cbind(rate) ~ ., FUN=median, data=long)
for(effect.type in c('logodds','smd')) {
    long.effect <- long[grep(effect.type,long$estimator),]
    out.file <- paste0(out.path,'m-',effect.type,'.png')
    png(out.file)
    plot(0,type='n',xlim=range(long.effect$m.exp),ylim=range(long.effect$rate),xlab=expression(paste(log[n]*m)),ylab='rate')
    estimators <- unique(long.effect$estimator)
    for(i in 1:length(estimators))
        lines(rate ~ m.exp, data=subset(long.effect,subset=(estimator==estimators[i])),lty=i)
    abline(h=-2,lty=3)
    legend.labels <- estimators |> gsub(pattern='(.+)test.stat$',replace='LR') |> gsub(pattern='(.+)test.stat.bartlett$',replace='LR Bartlett')
    legend('topright',lty=1:2,legend=legend.labels)
    dev.off()
}





## Simulation section--effect of deviations from normality
data.path <- ''
out.path <- ''
filelist <- dir(data.path)
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
print(length(filelist))
by.sim <- lapply(filelist, function(file) {
    load(paste0(data.path,file))
    list(by.n=by.n,ns=ns,m.exp=m.exp,skew=skew,kurt=kurt)
})
long <- lapply(by.sim, function(sim) {
    by.n <- sim$by.n
    ns <- sim$ns
    rates <- apply(by.n,1,function(ks)coef(lm(log(ks)~log(ns)))[2])
    settings <- c(m.exp=sim$m.exp,skew=sim$skew,kurt=sim$kurt)
    data.frame(estimator=names(rates),rate=unname(rates),t(settings))
})
long <- Reduce(rbind,long)
long <- aggregate(cbind(rate) ~ ., FUN=median, data=long)
long <- long[grep('bartlett',long$estimator),]
for(effect.type in c('logodds','smd')) {
    long.effect <- long[grep(effect.type,long$estimator),]
    out.file <- paste0(out.path,'normality-',effect.type,'.png')
    png(out.file)
    long.effect$kurt <- factor(long.effect$kurt,ordered=TRUE)
    kurts <- levels(long.effect$kurt)
    cols <- rev(gray.colors(length(kurts)))
    plot(0,type='n',xlab='skew of random effect',ylab='rate',xlim=range(long.effect$skew),ylim=range(long$rate))
    for(i in 1:length(kurts))
        lines(rate ~ skew, data=subset(long.effect,kurt==kurts[[i]]),col=cols[i])
    abline(h=-2,lty=2)
    legend('topright',col=cols,lty=1,legend=paste0('kurtosis=',kurts))
    dev.off()
}
