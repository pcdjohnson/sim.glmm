  # sim.glmm: Simulate responses from a GLMM.
  # Paul Johnson, Institute of BAHCM, University of Glasgow
  # 8th August 2014
  #
  #
  # DETAILS:
  #
  # The vector of responses is randomly simulated from a GLMM
  # and added to an input data set.
  # The values of the fixed effects, the random effects variances and 
  # covariances, and the response distribution are inputs. 
  #
  #
  # ARGUMENTS:
  #
  # mer.fit: A fitted GLMM object of class merMod. This includes models fitted by lmer 
  #   and glmer in the lme4 package. If mer.fit is supplied, no other argument should 
  #   be used, as all of the required inputs will be extracted from mer.fit.  
  #   
  # design.data: A data frame containing all the data except the response. Its columns 
  #   should correspond to names(fixed.eff), excluding the intercept, and names(rand.V).
  #
  # fixed.eff: A list of fixed effects. One element of the list must be called "intercept"
  #   or "(Intercept)".   
  #   The names of the other elements should correspond to variables in design.data.
  #   For example, to specify a model of the form y ~ 10 + 2*(sex=="Female") + 0.5*age
  #   you would use fixed.eff = list(intercept=10, sex=c("Male"=0,"Female"=2), age=0.5).
  #   This should work as long as design.data has a factor called "sex" with 
  #   levels "Male" and "Female" and a numeric variable called "age".
  #
  # rand.V: Either: (a) A vector of the variances of the random effects, where the
  #   names correspond to grouping factors in design.data; (b) A list of variance-covariance 
  #   matrices of the random effects, where the names of the list correspond to grouping factors
  #   in design.data. Currently only simple random effect structures are allowed: either
  #   random intercepts, or random intercepts-and-slopes. There is no limit on the number of
  #   random effects, and either crossed or nested structures are allowed.
  #   Option (b) allows covariances between random effects to be specified, which is necessary
  #   for random slopes-and-intercepts models because slopes and intercepts are almost 
  #   always correlated. ***Note that the function currently doesn't allow random slopes on 
  #   variables that are factors. See the examples for a simple workaround.***
  #   Where rand.V=NULL the resulting response will be simulated without 
  #   random effects, i.e. from a GLM.
  #
  # distribution: The response distribution. Currently has to be one of 
  #   "gaussian", "poisson", "binomial" and "negbinomial". For all but poisson 
  #   some additional information must be supplied:
  #   gaussian: SD must be suppied.
  #   binomial: For a binomial response (x successes out of n trials), design.data must have a
  #     column named "n" to specify the number of trials. For binary data (0 or 1),
  #     design.data$n should be a column of 1s.
  #   negative binomial: theta, the dispersion parameter, must be supplied. 
  #     theta is equivalent to "size" as defined in ?rnbinom. A negative binomial
  #     distribution with mean mu and dispersion parameter theta has
  #     variance mu + mu^2/theta.
  #  
  # SD: The residual standard deviation where distribution="gaussian".
  #
  # theta: The dispersion parameter where distribution="negbinomial".
  #
  # drop.effects: Deprecated, and only included for backward compatibility, so should be ignored.
  #
  # VALUE:
  # The input data frame with the response vector added as a new variable called "response".
  
    sim.glmm<-
      function(mer.fit=NULL, design.data=NULL, fixed.eff=NULL, rand.V=NULL,
        distribution=c("gaussian","poisson","binomial","negbinomial"),
        SD=NULL, theta=NULL, drop.effects=NULL)
      {

        # load required packages lme4 and MASS

          require(lme4)
          require(MASS)

        # if mer.fit is supplied, extract all inputs from it

          if(!is.null(mer.fit))
          {
            design.data <- model.frame(mer.fit)
            fixed.eff <- fixef(mer.fit)
            rand.V <- VarCorr(mer.fit)
            distribution <- as.character(family(mer.fit))[1]
            if(distribution == "binomial") design.data$n <- get("n",envir=slot(slot(mer.fit,"resp"),".xData"))
            SD <- attr(rand.V,"sc")
          }


        # add column of 1s. this is the "data value" for the intercept

          design.data$"(Intercept)"<-1


        # add the fixed effects to the data frame, first changing "intercept" to "(Intercept)"

          if(is.null(mer.fit))
          {
            names(fixed.eff)[names(fixed.eff)=="intercept"] <- "(Intercept)"
            fix.names<-names(fixed.eff)
            design.data[,paste(fix.names,".fixed",sep="")]<- 
              sapply(fix.names,
                function(fe)
                {
                  if(is.factor(design.data[,fe])) fixed.eff[[fe]][as.character(design.data[,fe])]
                    else fixed.eff[[fe]]*design.data[,fe]
                })
          } else 
            {
              design.data$combinedeffects.fixed <- model.matrix(mer.fit) %*% fixed.eff
              fix.names <- "combinedeffects"
            }
          eff.names <- paste(fix.names,".fixed",sep="")


        # if random effects are specified, simulate them 

          if(!is.null(rand.V))
          {
          
            # if rand.V is a vector, convert it to a list

              if(!is.list(rand.V)) rand.V <- as.list(rand.V)
              rand.names<-names(rand.V)

            # simulate the random effects

              a.list<-
                lapply(rand.names,
                  function(rn)
                  { 
                  	design.data[,rn] <- factor(design.data[,rn])
                    nt <- nlevels(design.data[,rn])
                    vcv <- rand.V[[rn]]
                    if(is.null(dim(vcv))) vcv <- structure(vcv, dim=c(1,1), dimnames=list("(Intercept)","(Intercept)"))
                    rownames(vcv)[rownames(vcv)=="intercept"] <- "(Intercept)"
                    colnames(vcv)[colnames(vcv)=="intercept"] <- "(Intercept)"
                    a <- mvrnorm(nt,rep(0,nrow(vcv)),vcv)
                    rownames(a) <- levels(design.data[,rn])
                    ax <- 
                      if(is.null(mer.fit)) {
                        a[design.data[,rn],colnames(vcv)] * design.data[,colnames(vcv)] } else {
                          a[design.data[,rn],colnames(vcv)] * model.matrix(mer.fit)[,colnames(vcv)] }
                    if(is.null(dim(ax))) dim(ax) <- c(nrow(design.data),nrow(vcv))
                    apply(ax,1,sum)
                  })
              names(a.list) <- rand.names


            # add the random effect residuals to the data frame
        
              design.data[,paste(rand.names,".random",sep="")]<- do.call("cbind",a.list[rand.names])
              eff.names <- c(eff.names,paste(rand.names,".random",sep=""))
          }


        # sum the fixed effects and random effect residuals to give the linear predictor
        # of the GLMM
        
          design.data$linear.predictor<-apply(design.data[,eff.names],1,sum)
          

        # simulate the response values from the linear predictor
  
          if(distribution[1]=="gaussian")
            design.data$response<-rnorm(n=1:nrow(design.data),mean=design.data$linear.predictor,sd=SD)
          if(distribution[1]=="poisson")
            design.data$response<-rpois(n=1:nrow(design.data),lambda=exp(design.data$linear.predictor))
          if(distribution[1]=="binomial")
            design.data$response<-rbinom(n=1:nrow(design.data),size=design.data$n,prob=plogis(design.data$linear.predictor))
          if(distribution[1]=="negbinomial")
            design.data$response<-rnbinom(n=1:nrow(design.data),size=theta,mu=exp(design.data$linear.predictor))
  

        # drop the intercept, fixed effects and random effect residuals from the data frame

          design.data<-design.data[,!(names(design.data) %in% c("(Intercept)","linear.predictor",eff.names))]
  

        # return the data frame including the simulated response vector 

          return(design.data)
      }



  # mor: Function to calculate median odds ratio (MOR) from a random effect variance in a logistic GLMM.
  # inv.mor: Inverse function to get variance from MOR.
  # For a poisson GLMM these functions will transform variances to median 
  # *rate* ratios (MRR), and vice versa.
  # see:
  #   Interpreting Parameters in the Logistic Regression Model with Random Effects
  #   Author(s): Klaus Larsen, J¯rgen Holm Petersen, Esben Budtz-J¯rgensen, Lars Endahl
  #   Biometrics, Vol. 56, No. 3 (Sep., 2000), pp. 909-914
  #
  #      mor=exp(sqrt(2*v)*qnorm(0.75))   (MOR function) 
  # =>   log(mor)=sqrt(2*v)*qnorm(0.75)
  # =>   (log(mor)/qnorm(0.75))^2=2*v
  # =>   v=((log(mor)/qnorm(0.75))^2)/2   (inverse MOR function)
  
    mor<-mrr<-function(v)exp(sqrt(2*v)*qnorm(0.75))
    inv.mor<-inv.mrr<-function(m)((log(m)/qnorm(0.75))^2)/2


if(F)
{

  # Poisson-lognormal example
  # simulate counts of tick burden on grouse chicks in a single year from a Poisson-lognormal GLMM,
  # loosely based on: 
  #   Elston et al. (2001). 
  #   Analysis of aggregation, a worked example: numbers of ticks on red grouse chicks. Parasitology, 122, 5639.
  #   http://abdn.ac.uk/lambin-group/Papers/Paper%2041%202001%20Elston%20Tick%20aggregation%20Parasitology.pdf
  # chicks are nested within broods, and broods within locations

    # simulate data set that defines sampling of chicks within broods within locations, 
    # assuming 3 chicks pr brood and 2 broods per location. N locations = 20.
      
      tickdata<-expand.grid(chick=1:3,brood=1:2,location=1:20)

    # make brood and chick ids unique (otherwise random effects will be simulated as crossed, not nested)

      tickdata$location<-factor(paste("loc",tickdata$location,sep=""))
      tickdata$brood<-factor(paste(tickdata$location,"brd",tickdata$brood,sep=""))
      tickdata$chick<-factor(paste(tickdata$brood,"chk",tickdata$chick,sep=""))
            
    # simulate tick counts with an average burden of 5 ticks per chick 
    # random effect variances are 2, 1 and 0.3 for location, brood and chick respectively
    
      tickdata<-
        sim.glmm(design.data = tickdata, 
          fixed.eff = list(intercept = log(5)),
          rand.V = c(location = 2, brood = 1, chick = 0.3),
          distribution = 'poisson')
      
    # plot counts and fit GLMM

      plot(response~jitter(as.numeric(location),factor=0.5),pch=21,bg=as.numeric(brood),data=tickdata)
      library(lme4)
      glmer(response~(1|location)+(1|brood)+(1|chick),family='poisson',data=tickdata)

    
  # lognormal-poisson example: trial of mosquito traps 
  
  # simulate mosquito abundance data from a field trial of 6 types of trap in 6 huts
  # huts A-F, weeks w1-w6 and experimental traps U-Z
  # six counts are taken in each hut on days 1-6 of each week
  # traps are rotated through huts weekly, 6 weeks every trap has been tested in every hut for 1 week (a Latin square design) 

    hut.data<-expand.grid(hut=LETTERS[1:6],week=paste("w",1:6,sep=""),obs=1:6)

  # rotate trap types through huts weekly 
  
    hut.data$trap<-factor(LETTERS[21:26][unlist(by(hut.data,hut.data$week,function(x) 1 + (0:5 + unique(as.numeric(x$week))) %% length(levels(x$week))))])

  # give each row a unique indentifier to allow lognormal overdispersion to be simulated      

    hut.data$row.id<-factor(paste("row",1:nrow(hut.data),sep=""))

  # simulate abundance data 

    hut.data<-
      sim.glmm(
        design.data=hut.data,
        fixed.eff=
          list(
            intercept=log(5),                     # average abundance is 5 mosquitoes per night in the control (reference) trap
            trap=
              log(                                # note that all the fixed effects are logged, because the parameters are estimated on the log scale
                c(U=1,                            # U is the reference category, so must have a relative rate of 1
                  V=3,                            # trap V catches 3 times as many mosquitoes as U, on average
                  W=1.5,X=1.5,Y=1.5,Z=1.5))),     # all the other traps catch 50% more than the control (U)
        rand.V=
          inv.mor(
            c(row.id=2,                           # the overdispersion median rate ratio (MRR) is 2
              hut=1.3,                            # there is variation in abundance between huts (MRR=1.3)
              week=1.5)),                         # there is also variation between weeks (MRR=1.5)
        distribution="poisson")                   # we are simulating a Poisson response
  
    # view and analyse hut data, testing for a difference between trap V and trap U

      par(mfrow=c(2,2))  
      hist(hut.data$response,xlab="Abundance")
      boxplot(response~trap,data=hut.data,ylab="Abundance",xlab="Trap")
      boxplot(response~hut,data=hut.data,ylab="Abundance",xlab="Hut")
      boxplot(response~week,data=hut.data,ylab="Abundance",xlab="Week")
    
      library(lme4)
      (mod.pois <-
        glmer(response~trap+(1|hut)+(1|week)+(1|row.id),family="poisson",
          data=hut.data))
      exp(fixef(mod.pois))

    # ... the "trapV" row of the "Pr(>|z|)" column of the fixed effects results table gives a p-value for
    # a test of the null hypothesis that U and V have the same abundance.
    # if you repeatedly run this simulation you should find that p < 0.05 close to 100% of the time,
    # that is, power is close to 100%. That could be considered wastefully excessive,  
    # and might motivate reducing the number of observations collected. 
    # however you should find that power is inadequate (~50%) for each of traps W-Z.
    

  # binomial example: simulate mortality data

  # now we are interesting in comparing mortality between the different traps, 
  # i.e. the number of n trapped mosquitoes that die.

  # we need a column called n to store the denominator (n mosquitoes cauhgt)

    hut.data$n<-hut.data$response

  # simulate the number that died

    hut.data<-
      sim.glmm(
        design.data=hut.data,
        fixed.eff=
          list(
            intercept=qlogis(0.7),                # mortality is 70% in the control (reference) trap
            trap=
              log(                                # note that all the fixed effects are logged, because the parameters are estimated on the log scale
                c(U=1,                            # U is the reference category, so must have an odds ratio of 1
                  V=2,                            # the odds of mortality is twice that in U
                  W=1.5,X=1.5,Y=1.5,Z=1.5))),     # in all the other traps the odds ratio is 1.5
        rand.V=
          inv.mor(
            c(row.id=2,                           # the overdispersion median odds ratio (MOR) is 2
              hut=1.3,                            # there is variation in mortality between huts (MOR=1.3)
              week=1.5)),                         # there is also variation between weeks (MOR=1.5)
        distribution="binomial")                  # we are simulating a binomial response
  
    # view and analyse hut data, testing for a difference between trap V and trap U

      par(mfrow=c(2,2))  
      hist(hut.data$response/hut.data$n,xlab="Mortality")
      boxplot(response/n~trap,data=hut.data,ylab="Mortality",xlab="Trap")
      boxplot(response/n~hut,data=hut.data,ylab="Mortality",xlab="Hut")
      boxplot(response/n~week,data=hut.data,ylab="Mortality",xlab="Week")
    
      library(lme4)
      (mod.bin<-glmer(cbind(response,n-response)~trap+(1|hut)+(1|week)+(1|row.id),
        family="binomial",data=hut.data))
      plogis(fixef(mod.bin)[1])   # estimated mortality in the control trap 
      exp(fixef(mod.bin)[-1])     # odds ratio estimates for the other traps


  # we could also simulate a gaussian response

    hut.data<-
      sim.glmm(
        design.data=hut.data,
        fixed.eff=
          list(
            intercept=10,                         # mean = 10 in the reference trap
            trap=
                c(U=0,                            # U is the reference category, so must have a regression coefficient of 0
                  V=1,W=1,X=1,Y=1,Z=1)),          # all the other traps raise the measure by 1 unit
        rand.V=c(hut=1,week=1),                   # there is variation between huts and between weeks (var=1)
        distribution="gaussian",                  # we are simulating a Gaussian response
        SD=2)                                     # the residual SD is 2

    # view and analyse hut data, testing for a difference between trap V and trap U

      par(mfrow=c(2,2))  
      hist(hut.data$response,xlab="Response")
      boxplot(response~trap,data=hut.data,ylab="Response",xlab="Trap")
      boxplot(response~hut,data=hut.data,ylab="Response",xlab="Hut")
      boxplot(response~week,data=hut.data,ylab="Response",xlab="Week")
    
      library(lme4)
      (mod.gaus<-lmer(response~trap+(1|hut)+(1|week),data=hut.data))
      fixef(mod.gaus)
    

  # returning to abundance, we can also simulate overdispersed counts from the negative binomial distribution

    hut.data<-
      sim.glmm(
        design.data=hut.data,
        fixed.eff=
          list(
            intercept=log(5),                     # average abundance is 5 mosquitoes per night in the control (reference) trap
            trap=
              log(                                # note that all the fixed effects are logged, because the parameters are estimated on the log scale
                c(U=1,                            # U is the reference category, so must have a relative rate of 1
                  V=3,                            # trap V catches 3 times as many mosquitoes as U, on average
                  W=1.5,X=1.5,Y=1.5,Z=1.5))),     # all the other traps catch 50% more than the control (U)
        rand.V=
          inv.mor(
            c(hut=1.3,                            # there is variation in abundance between huts (MRR=1.3)
              week=1.5)),                         # there is also variation between weeks (MRR=1.5)
        distribution="negbinomial",               # we are simulating a negative binomial response
        theta=0.5)                                # overdispersion is introduced via the theta parameter,  
                                                  # rather than as a random effect in the lognormal Poisson above


    # view and analyse hut data, testing for a difference between trap V and trap U

      # plot data

        par(mfrow=c(2,2))  
        hist(hut.data$response,xlab="Abundance") # similar amount of overdispersion to the lognormal Poisson example above
        boxplot(response~trap,data=hut.data,ylab="Abundance",xlab="Trap")
        boxplot(response~hut,data=hut.data,ylab="Abundance",xlab="Hut")
        boxplot(response~week,data=hut.data,ylab="Abundance",xlab="Week")

      # load glmmADMB package (http://glmmadmb.r-forge.r-project.org/) and prepare data

        library(glmmADMB)  # (note that this analysis failed when running glmmadmb on 32-bit R 2.15.3 for Mac. Worked fine on 64 bit).

      # fit negative binomial mixed model

        mod.nbin<-glmmadmb(response~trap+(1|hut)+(1|week),family="nbinom2",data=hut.data)
        summary(mod.nbin)
        mod.nbin$alpha # glmmadmb calls the overdispersion parameter "alpha" rather than "theta"
        exp(fixef(mod.nbin))

  
  # simulation from a random slopes model using the sleepstudy data 
  # from the lme4 package (see ?sleepstudy for details)

    # illustrate variation in slope between subjects

      library(lattice)
      library(lme4)
      xyplot(Reaction ~ Days | Subject, sleepstudy,
        panel=
          function(x,y){
            panel.xyplot(x,y)
            if(length(unique(x))>1) panel.abline(lm(y~x))
          })  

    # fit random slopes model

      fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

    # use the estimates from the fitted model to parameterise the simulation model
    # this can be done explicitly, by extracting the estimates and supplying them as arguments

      sim.glmm(design.data=sleepstudy, fixed.eff=fixef(fm1), rand.V=VarCorr(fm1), 
        distribution="gaussian", SD=attr(VarCorr(fm1),"sc"))

    # but if the model was fitted with lmer or glmer, the function can extract the estimates automatically

      sim.glmm(mer.fit=fm1)
      
    # check that the data is being simulated from the correct model by estimating the parameters 
    # from multiple simulated data sets and plotting the estimates with the input parameters 

      sim.res <- 
        sapply(1:100,function(i){
          print(i)
          sim.fm1  <- lmer(response ~ Days + (Days|Subject), sim.glmm(mer.fit=fm1))
          c(fixef(sim.fm1),unlist(VarCorr(sim.fm1)),SD=attr(VarCorr(sim.fm1),"sc"))
          })

      boxplot(t(sim.res), main="Boxplot of fixed and random effect\nestimates from 100 simulated data sets")
      points(c(fixef(fm1),unlist(VarCorr(fm1)),SD=attr(VarCorr(fm1),"sc")), pch="-", col="red", cex=4)
      legend("topright",legend="True values",pch="-",pt.cex=4, col="red")


  # the same example, but with a random slope on a factor fixed effect
  # currently sim.glmm doesn't handle random slopes for factors, so the 
  # following is a workaround
  
    # convert Days to a binary factor and fit model
  
      sleepstudy$fDays <- factor(sleepstudy$Days > 4.5, c(FALSE, TRUE), c("Lo", "Hi"))
      table(sleepstudy$fDays, sleepstudy$Subject)
      (fm1f <- lmer(Reaction ~ fDays + (fDays | Subject), sleepstudy))

    # use the estimates from the fitted model to parameterise the simulation model
      
      # this can be done directly from the merMod object: 

        sim.glmm(fm1f)
 
      # but if we wanted to change the parameters we would need to be able to specify the parameters individually 
      # which gives an error
       
        sim.glmm(design.data=sleepstudy, 
          fixed.eff=list(intercept=271.6, fDays=c(Lo=0, Hi=53.76)),
          rand.V=VarCorr(fm1f), 
          distribution="gaussian", SD=attr(VarCorr(fm1f),"sc"))

      # a simple workaround is to represent the factor as an indicator variable
      # (or variables if there are more than two levels):
      
        sleepstudy <- cbind(sleepstudy, model.matrix(~ fDays, data=sleepstudy))
      
      # the simulation code above should now work. i will apply this fix internally when i have time.
       

  # a poisson random slopes example
  # this example uses the Owls data which is in the glmmADMB package (see ?Owls for details)
  
    # load glmmADMB package (http://glmmadmb.r-forge.r-project.org/) and prepare data

      library(glmmADMB)
      Owls$obs <- factor(1:nrow(Owls))            # to fit observation-level random effect
      Owls$ArrivalTimeC <- Owls$ArrivalTime - 24   # centre the arrival times at midnight


    # illustrate variation in slope between nests

      xyplot(SiblingNegotiation ~ ArrivalTimeC | Nest, Owls,
        panel=
          function(x,y){
            panel.xyplot(x,y)
            if(length(unique(x))>1) panel.abline(lm(y~x))
          })  

    # fit random slopes model

      owlmod.rs <- 
        glmer(SiblingNegotiation ~ ArrivalTimeC + (ArrivalTimeC|Nest) + (1|obs),
          family="poisson", data=Owls)

    # fit simulate from fitted model and fit model on simulated data

      (sim.owlmod.rs <- 
        glmer(response ~ ArrivalTimeC + (ArrivalTimeC|Nest) + (1|obs),
          family="poisson", data=sim.glmm(owlmod.rs)))


    # check that the data is being simulated from the correct model by estimating the parameters 
    # from multiple simulated data sets and plotting the estimates with the input parameters 

      sim.res <- 
        sapply(1:20,function(i){
          print(i)
          sim.owlmod.rs  <- 
            glmer(response ~ ArrivalTimeC + (ArrivalTimeC|Nest) + (1|obs),
              family="poisson", data=sim.glmm(owlmod.rs))
          c(fixef(sim.owlmod.rs),unlist(VarCorr(sim.owlmod.rs)))
          })

      dev.off()
      boxplot(t(sim.res), main="Boxplot of fixed and random effect\nestimates from 20 simulated data sets")
      points(c(fixef(owlmod.rs),unlist(VarCorr(owlmod.rs))), pch="-", col="red", cex=4)
      legend("topright",legend="True values",pch="-",pt.cex=4, col="red")


}



