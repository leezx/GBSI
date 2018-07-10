[Practical Guide to Principal Component Analysis (PCA) in R & Python](https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/)











# too difficult for me, hard to interpret!

distance-based methods ( Mahalanobis )

projection pursuit 



kurtosis weights 

This enables us to assign weights to each component according to how likely we

think it is to reveal the outliers.





## [mvoutlier](https://cran.r-project.org/web/packages/mvoutlier/index.html): Multivariate Outlier Detection Based on Robust Methods

```R
pcout(x, makeplot = FALSE, explvar = 0.99, crit.M1 = 1/3, crit.c1 = 2.5,
crit.M2 = 1/4, crit.c2 = 0.99, cs = 0.25, outbound = 0.25, ...)
```

crit.M1 a numeric value between 0 and 1 indicating the quantile to be used as lower boundary for location outlier detection (default to 1/3) 

crit.c1 a positive numeric value used for determining the upper boundary for location outlier detection (default to 2.5) 

crit.M2 a numeric value between 0 and 1 indicating the quantile to be used as lower boundary for scatter outlier detection (default to 1/4) 

crit.c2 a numeric value between 0 and 1 indicating the quantile to be used as upper boundary for scatter outlier detection (default to 0.99)

cs a numeric value indicating the scaling constant for combined location and scatter weights (default to 0.25)

```R
# geochemical data from northern Europe
data(bsstop)
x=bsstop[,5:14]
# identify multivariate outliers
x.out=pcout(x,makeplot=T)
# visualize multivariate outliers in the map
op <- par(mfrow=c(1,2))
data(bss.background)
pbb(asp=1)
points(bsstop$XCOO,bsstop$YCOO,pch=16,col=x.out$wfinal01+2)
title("Outlier detection based on pcout")
legend("topleft",legend=c("potential outliers","regular observations"),pch=16,col=c(2,3))
# compare with outlier detection based on MCD:
x.mcd <- robustbase::covMcd(x)
pbb(asp=1)
points(bsstop$XCOO,bsstop$YCOO,pch=16,col=x.mcd$mcd.wt+2)
title("Outlier detection based on MCD")
legend("topleft",legend=c("potential outliers","regular observations"),pch=16,col=c(2,3))
par(op)
```

Result 

**wfinal01** 0/1 vector with final weights for each observation; weight 0 indicates potential multivariate outliers. 

**wfinal** numeric vector with final weights for each observation; small values indicate potential multivariate outliers.

**wloc** numeric vector with weights for each observation; small values indicate potential location outliers. 

**wscat** numeric vector with weights for each observation; small values indicate potential scatter outliers. 

**x.dist1** numeric vector with distances for location outlier detection. 

**x.dist2** numeric vector with distances for scatter outlier detection. 

**M1** upper boundary for assigning weight 1 in location outlier detection. 

**const1** lower boundary for assigning weight 0 in location outlier detection. 

**M2** upper boundary for assigning weight 1 in scatter outlier detection. 

**const2** lower boundary for assigning weight 0 in scatter outlier detection.

```R
function (x, makeplot = FALSE, explvar = 0.99, crit.M1 = 1/3, 
    crit.c1 = 2.5, crit.M2 = 1/4, crit.c2 = 0.99, cs = 0.25, 
    outbound = 0.25, ...) 
{
    p = ncol(x)
    n = nrow(x)
    x.mad = apply(x, 2, mad) # Median absolute deviation
    if (any(x.mad == 0)) 
        stop("More than 50% equal values in one or more variables!")
    x.sc <- scale(x, apply(x, 2, median), x.mad)
    x.sc <- 
    x.svd <- svd(scale(x.sc, TRUE, FALSE))
    a <- x.svd$d^2/(n - 1)
    p1 <- (1:p)[(cumsum(a)/sum(a) > explvar)][1]
    x.pc <- x.sc %*% x.svd$v[, 1:p1]
    xpc.sc <- scale(x.pc, apply(x.pc, 2, median), apply(x.pc, 
        2, mad))
    wp <- abs(apply(xpc.sc^4, 2, mean) - 3)
    xpcw.sc <- xpc.sc %*% diag(wp/sum(wp))
    xpc.norm <- sqrt(apply(xpcw.sc^2, 1, sum))
    x.dist1 <- xpc.norm * sqrt(qchisq(0.5, p1))/median(xpc.norm)
    M1 <- quantile(x.dist1, crit.M1)
    const1 <- median(x.dist1) + crit.c1 * mad(x.dist1)
    w1 <- (1 - ((x.dist1 - M1)/(const1 - M1))^2)^2
    w1[x.dist1 < M1] <- 1
    w1[x.dist1 > const1] <- 0
    xpc.norm <- sqrt(apply(xpc.sc^2, 1, sum))
    x.dist2 <- xpc.norm * sqrt(qchisq(0.5, p1))/median(xpc.norm)
    M2 <- sqrt(qchisq(crit.M2, p1))
    const2 <- sqrt(qchisq(crit.c2, p1))
    w2 <- (1 - ((x.dist2 - M2)/(const2 - M2))^2)^2
    w2[x.dist2 < M2] <- 1
    w2[x.dist2 > const2] <- 0
    wfinal <- (w1 + cs) * (w2 + cs)/((1 + cs)^2)
    wfinal01 <- round(wfinal + 0.5 - outbound)
    if (makeplot) {
        op <- par(mfrow = c(3, 2), mar = c(4, 4, 2, 2))
        on.exit(par(op))
        plot(x.dist1, xlab = "Index", ylab = "Distance (location)", 
            ...)
        abline(h = const1)
        abline(h = M1, lty = 2)
        plot(w1, xlab = "Index", ylab = "Weight (location)", 
            ylim = c(0, 1), ...)
        abline(h = 0)
        abline(h = 1, lty = 2)
        plot(x.dist2, xlab = "Index", ylab = "Distance (scatter)", 
            ...)
        abline(h = const2)
        abline(h = M2, lty = 2)
        plot(w2, xlab = "Index", ylab = "Weight (scatter)", ylim = c(0, 
            1), ...)
        abline(h = 0)
        abline(h = 1, lty = 2)
        plot(wfinal, xlab = "Index", ylab = "Weight (combined)", 
            ylim = c(0, 1), ...)
        abline(h = cs)
        plot(wfinal01, xlab = "Index", ylab = "Final 0/1 weight", 
            ylim = c(0, 1), ...)
    }
    list(wfinal01 = wfinal01, wfinal = wfinal, wloc = w1, wscat = w2, 
        x.dist1 = x.dist1, x.dist2 = x.dist2, M1 = M1, const1 = const1, 
        M2 = M2, const2 = const2)
}
```



x.out=pcout(t(simulate_data()),makeplot=FALSE)



