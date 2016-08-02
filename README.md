sim.glmm
========

R function to simulate from a generalized linear mixed model (GLMM).

*****
**_Since 2016-08-02 sim.glmm has been available in an R package, [GLMMmisc](https://github.com/pcdjohnson/GLMMmisc "GLMMmisc R package"):_**

https://github.com/pcdjohnson/GLMMmisc

**_I'll leave the old version here, so the code below will still work, but I recommend using the up-to-date versions in the GLMMmisc package._**
*****


The script can be sourced directly from github using the RCurl package:

    library(RCurl)
    options(RCurlOptions=list(cainfo=system.file("CurlSSL","cacert.pem",package="RCurl")))
    eval(expr=parse(text=getURL("https://raw.githubusercontent.com/pcdjohnson/sim.glmm/master/sim.glmm.R")))

See header of sim.glmm.R for details of how to use the function. Several examples are given below the function (but not run on sourcing). 

An article and tutorial on power analysis using this function are available here:
Johnson, P.C.D., Barry, S.J.E., Ferguson, H.M. & Müller, P. (2015). Power analysis for generalized linear mixed models in ecology and evolution. Methods in Ecology and Evolution. http://dx.doi.org/10.1111/2041-210X.12306
