sim.glmm
========

R function to simulate from a generalized linear mixed model (GLMM).

The script can be sourced directly from github using the RCurl package:

    library(RCurl)
    options(RCurlOptions=list(cainfo=system.file("CurlSSL","cacert.pem",package="RCurl")))
    eval(expr=parse(text=getURL("https://raw.githubusercontent.com/pcdjohnson/sim.glmm/master/sim.glmm.R")))

See header of sim.glmm.R for details of how to use the function. Several examples are given below the function (but not run on sourcing). 

An article and tutorial on power analysis using simulation are available here:
Johnson, P.C.D., Barry, S.J.E., Ferguson, H.M. & Müller, P. (2015). Power analysis for generalized linear mixed models in ecology and evolution. Methods in Ecology and Evolution. http://dx.doi.org/10.1111/2041-210X.12306
