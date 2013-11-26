  # sim.glmm: simulate responses from a GLMM.
  # Paul Johnson, Institute of BAH&CM, University of Glasgow
  # 30th September 2013
  #
  # v4.2:  
  #   - changes from v1: now allows random effect variances to be supplied
  #     as a vector as well as a matrix
  #   - changes from v2: fixed bug
  #   - changes from v3: add negative binomial ditribution
  #   - changes from v4 and 4.1: 
  #       fixed bug in formation of var-covar matrices of dimension 1
  #       changed unix/mac end of line so that Windows notepad shows EOL
  #
  # DETAILS:
  # the vector of responses is randomly simulated from a GLMM
  # and added to an input data set.
  # the values of the fixed effects, the random effects variances and 
  # covariances, and the response distribution are inputs. 
  #
  # ARGUMENTS:
  # design.data: a data frame containing all the data except the response.
  # fixed.eff: a list of fixed effects. fixed.eff$intercept has to be supplied.
  #   the names of the other elements should correspond to variables in design.data.
  #   for example, to specify a model of the form y ~ 10 + 2*(sex=="Female") + 0.5*age
  #   you would use fixed.eff=list(intercept=10,sex=c("Male"=0,"Female"=2),age=0.5).
  #   this should work as long as design.data has a factor called sex with 
  #   levels "Male" and "Female" and a numeric variable called age.
  #   if there are no fixed effects (other than the intercept) then fixed.eff should have only
  #   one element, "intercept".
  # rand.V: either (a) a vector of the variances of the random effects, where the
  #   names correspond to factors in design.data; or (b) the variance-covariance 
  #   matrix of the random effects, where the column and 
  #   row names correspond to factors in design.data.
  #   option (b) allows covariances between random effects to be specified.
  #   if rand.V is null the resulting response will be simulated without 
  #   random effects, i.e. from a GLM.
  # distribution: specifies the response distribution. currently has to be 
  #   gaussian, poisson, binomial or negative binomial. for all but poisson 
  #   some additional information must be supplied:
  #   gaussian: SD must be suppied.
  #   binomial: for a binomial response (x successes out of n trials), design.data must have a
  #     column named n to specify the number of trials. for binary data (0 or 1)
  #     use design.data$n should be a column of 1s. i haven't dealt with the 
  #     (realistic) scenario where both x and n vary randomly, although
  #     i guess the n could be simulated first as poisson then the x as binomial.
  #     how to analyse this properly is also not clear to me. The approach I would
  #     take is to simulate variation in n but ignore it in the analysis. This
  #     seems to work pretty well.
  #   negative binomial: theta, the overdispersion parameter, must be supplied. 
  #     theta is equivalent to "size" as defined in ?rnbinom. a poisson distribution
  #     with mean mu has variance mu. a negative binomial distribution with mean mu 
  #     and overdispersion parameter theta has variance mu + mu^2/theta.  
  # SD: the residual standard deviation where the distribution is gaussian.
  # drop.effects. default is TRUE. if drop.effects is FALSE, the fixed effects and  
  #    the random effect residuals will be output as additional rows in the data frame.
  #
  # VALUE:
  # the input data frame with the response vector added as a new variable called "response".
  
    sim.glmm<-
      function(design.data,fixed.eff,rand.V=NULL,
        distribution=c("gaussian","poisson","binomial","negbinomial"),
    		SD=NULL,theta=NULL,drop.effects=TRUE)
      {

        # load required package MASS

          require(MASS)

        # add column of 1s. this is the "data value" for the intercept

          design.data$intercept<-1
          
        # if no random effects are specified, set the variances to zero

          if(is.null(rand.V))
          {
            rand.V<-matrix(0,dimnames=rep(list("dummy"),2))
            design.data$dummy<-1
          }
          
        # if rand.V is a vector, convert it to a variance covariance matrix

          if(is.null(dim(rand.V)))
          {
            rand.V<-structure(diag(rand.V,nrow=length(rand.V)),dimnames=list(names(rand.V),names(rand.V)))
          }

          rand.names<-colnames(rand.V)
          if(any(!sapply(design.data[,rand.names],is.factor))) stop("All random effect variables must be factors")

        # simulate the residuals of the random effects

          nt<-max(c(2,sapply(rand.names,function(r) nlevels(design.data[,r]))))
          b<-mvrnorm(nt,rep(0,nrow(rand.V)),rand.V)

        # add the random effect residuals to the data frame
        
          design.data[,paste(rand.names,".mean",sep="")]<-
            sapply(rand.names,function(r) b[as.numeric(design.data[,r]),r])

        # add the fixed effects to the data frame

          fix.names<-names(fixed.eff)
          design.data[,paste(fix.names,".eff",sep="")]<- 
            sapply(fix.names,
              function(fe)
              {
                if(is.factor(design.data[,fe])) fixed.eff[[fe]][as.character(design.data[,fe])]
                  else fixed.eff[[fe]]*design.data[,fe]
              })

        # sum the fixed effects and random effect residuals to give the linear predictor
        # of the GLMM
        
          eff.names<-c(paste(fix.names,".eff",sep=""),paste(rand.names,".mean",sep=""))
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
  
        # optionally drop the intercept, fixed effects and random effect residuals from the data frame

          if(drop.effects) design.data<-design.data[,!(names(design.data) %in% c("intercept","linear.predictor",eff.names))]
          design.data$dummy<-design.data$dummy.mean<-NULL
  
        # return the data frame including the simulated response vector 

          return(design.data)
      }



  # mor: function to calculate median odds ratio (MOR) from a variance component from a logistic GLMM.
  # inv.mor: inverse function to get variance from MOR.
  # for a poisson GLMM these functions will transform variances to median 
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
  #   Analysis of aggregation, a worked example: numbers of ticks on red grouse chicks. Parasitology, 122, 563–9.
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
      (mod.pois<-glmer(response~trap+(1|hut)+(1|week)+(1|row.id),family="poisson",data=hut.data))
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
      (mod.bin<-glmer(cbind(response,n-response)~trap+(1|hut)+(1|week)+(1|row.id),family="binomial",data=hut.data))
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
      (mod.gaus<-glmer(response~trap+(1|hut)+(1|week),family="gaussian",data=hut.data))
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

      par(mfrow=c(2,2))  
      hist(hut.data$response,xlab="Abundance") # similar amount of overdispersion to the lognormal Poisson example above
      boxplot(response~trap,data=hut.data,ylab="Abundance",xlab="Trap")
      boxplot(response~hut,data=hut.data,ylab="Abundance",xlab="Hut")
      boxplot(response~week,data=hut.data,ylab="Abundance",xlab="Week")

      library(glmmADMB)  # (note that this analysis failed when running glmmadmb on 32-bit R 2.15.3 for Mac. Worked fine on 64 bit).
      mod.nbin<-glmmadmb(response~trap+(1|hut)+(1|week),family="nbinom2",data=hut.data)
      summary(mod.nbin)
      mod.nbin$alpha # glmmadmb calls the overdispersion parameter "alpha" rather than "theta"
      exp(fixef(mod.nbin))

}
