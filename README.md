sim.glmm
========

R function to simulate from a generalized linear mixed model (GLMM).

The script can be sourced directly from github using the RCurl package:

    library(RCurl)
    options(RCurlOptions=list(cainfo=system.file("CurlSSL","cacert.pem",package="RCurl")))
    eval(expr=parse(text=getURL("https://raw.githubusercontent.com/pcdjohnson/sim.glmm/master/sim.glmm.R")))

See header of sim.glmm.R for details of how to use the function. Several examples are given below the function (but not run on sourcing).
