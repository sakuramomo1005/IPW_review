########################################################################################################################################
#  * Mimic one simulation scenario in Hossain's paper
#    * scenario 1 in the paper
#  * k=50; M=50; number of clusters equals to 50 and number of participants
#    in each clusters eaquals to 50.
#  * cluster sizes are the same
#
# Started by Lanqiu Yao Sept 13th 2018
#
# Edited by Liz Turner Sept 15th 2018 to make the following corrections:
# 
# 1. Amend input sigma_delta (i.e. SD) to be sqrt(0.2) like Hossain (rather than the value of 0.2 in previous version since he used a variance, not SD, of 0.2)
# Note that this corresponds to the outcome RE variance i.e. will now decrease the outcome ICC to be like Hossain
# 2. In function "mypool" 
#     2A.  Amend "v <- mean(sd0,na.rm=TRUE)" to be "v <- mean(sd0^2,na.rm=TRUE)" 
#         since the average within-imputation variance should be the average of the variances 
#         (i.e. of squared standard error and not only of the standard errors) 
#     2B. Changed term v to be the more commonly used "W" since W refers to within variance
#     2C. Create a variable "num_actual <- num-na_times" so can use throughout
#     2D. Amend "v_hat <- v+(1+1/num)*B" to "v_hat <- W+(1+1/num_actual)*B"
#     2E. Amended the "(1/df_part)" part of "df_obs <- df_com*((df_com+1)/(df_com+3))*(1/df_part)" 
#         to (1+df_denom*(1/W))
#     2F. Amended "return(list(mean=m,std=v_hat))" after "if(print=='no')" to be "return(list(mean=m,std=v_hat,df_t=df_t,df_adj_t=df_adj_t))"
#     2G. Amended so that used sqrt_vhat for all inference since we need an estimate of the SE and not of SE^2 to be used
#     2H. Changed B <- sd(mean0,na.rm=TRUE) to B <- var(mean0,na.rm=TRUE)
#     2I. Need to use  v_hat_sqrt<-sqrt(v_hat) for inference, not v_hat
#     
# 3. Moved "Nimp <- 15" to be next to all the other fixed inputs
# 4. Saved output of df_adj_t_incorrect which is what we had previously programmed thinking to be the correct df_adj_t needed for inference from MMI
#
########################################################################################################################################

# Package Use: geepack

# Install packages 
install.packages("jomo")
install.packages("geepack")

# Load in library
library(jomo)
library(geepack)

# Functions
## 1. Expit function
expit <- function(x){y <- exp(x)/(1+exp(x));return(y)}

## 2. Pool function
# See that mean0 and sd0 are the vectors of point estiamtes and of SE of point estiamtes from the analysis of each of the "num" imputed data sets 
# i..e of the "num" complete data sets that were imputed
mypool <- function(mean0,sd0,num=5,print='no',J=50){
  
  na_times <- sum(is.na(mean0))
  num_actual <- num-na_times
  
  # see that we take in to account na_times in next commands for mean and sd since use "na.rm=TRUE"
  m <- mean(mean0,na.rm=TRUE) # this is the MI estimate of the beta parameter
  ##  v <- mean(sd0,na.rm=TRUE) ## THIS IS AMENDED IN NEXT LINE
  W <- mean(sd0^2,na.rm=TRUE) # estimate of average wihtin-imputation variance i.e. based on the SE^2 of the beta parameter from each fitted model
  ## B <- sd(mean0,na.rm=TRUE) # estimate of between-imputation variance i.e. empirical SD of the point estimates of beta parameter
  B <- var(mean0,na.rm=TRUE) # estimate of between-imputation variance i.e. empirical SD of the point estimates of beta parameter
  #v_hat <- v+(1+1/num)*B ## THIS IS AMENDED IN NEXT LINE
  v_hat <- W+(1+1/num_actual)*B # estimate of total variance i.e. will need to take the sqrt of this to use for inference in making confidence intervals etc.
  v_hat_sqrt<-sqrt(v_hat)
  
  # Testing based on naive normal distribution assumption
  #l <- m-1.96*v_hat
  #u <- m+1.96*v_hat
  l <- m-1.96*v_hat_sqrt
  u <- m+1.96*v_hat_sqrt
  
  
  # Testing based on standard results from MI literature
  # i.e. df of t distribution for testing based on standard results from MI literature
  df_denom <- (1+1/num_actual)*B
  df_part <- 1+W/df_denom
  df_t <- (num_actual-1)*df_part^2 # calculate the df of the t distribution we should compare this to
  #l_t <- m-qt(0.975,df_t)*v_hat
  #u_t <- m+qt(0.975,df_t)*v_hat
  l_t <- m-qt(0.975,df_t)*v_hat_sqrt
  u_t <- m+qt(0.975,df_t)*v_hat_sqrt
  
  # Testing based on results from MMI literature
  # df of t distribution for testing based on results in re adjustment for MMI feature
  # i.e. need to reduce the df based on a formula from Barnard and Rubin (1999)
  df_com <- 2*J - 2 #df for full data - calculated this way even have same # clusters in each arm
  # note that J is used to denote # clusters in each arm
  ## df_obs <- df_com*((df_com+1)/(df_com+3))*(1/df_part) ## AMENDED SINCE 1/df_part is wrong expression
  df_obs_incorrect <- df_com*((df_com+1)/(df_com+3))*(1/df_part) 
  parenthesis <- 1+df_denom*(1/W) ## see that this is 1+ [(Q+1)/Q]*(B/W) because df_denom is  (1+1/Q)*B or, equivalently, B*(Q+1)/Q
  df_obs <- df_com*((df_com+1)/(df_com+3))*(1/parenthesis) 
  df_adj_t_incorrect <- 1/(1/df_t + 1/df_obs_incorrect) 
  df_adj_t <- 1/(1/df_t + 1/df_obs) 
  
  # Now print all of them to screen
  print("Standard t df and Barnard/Rubin adjusted t df");print(c(df_t, df_adj_t_incorrect, df_adj_t))
  print("97.5% quantiles from standard t df and Barnard/Rubin adjusted t df");print(c(qt(0.975,df_t),qt(0.975,df_adj_t_incorrect), qt(0.975,df_adj_t)))
  #l_t <- m-qt(0.975,df_t)*v_hat
  #u_t <- m+qt(0.975,df_t)*v_hat
  #l_adj_t <- m-qt(0.975,df_adj_t)*v_hat
  #u_adj_t <- m+qt(0.975,df_adj_t)*v_hat
  
  l_t <- m-qt(0.975,df_t)*v_hat_sqrt
  u_t <- m+qt(0.975,df_t)*v_hat_sqrt
  l_adj_t <- m-qt(0.975,df_adj_t)*v_hat_sqrt
  u_adj_t <- m+qt(0.975,df_adj_t)*v_hat_sqrt


  if(print=='no'){
    #return(list(mean=m,std=v_hat,df_t=df_t,df_adj_t_incorrect=df_adj_t_incorrect,df_adj_t=df_adj_t)) ## AMENDED IN NEXT LINE
    return(list(mean=m,std=v_hat_sqrt,df_t=df_t,df_adj_t_incorrect=df_adj_t_incorrect,df_adj_t=df_adj_t))
  }
  if(print=='yes'){
    print('mean (95% CI)')
    print(paste(round(m,3)," (",round(l,3),',',round(u,3),')',sep=''))
    print(paste(round(m,3)," (",round(l_t,3),',',round(u_t,3),')',sep=''))
    print(paste(round(m,3)," (",round(l_adj_t,3),',',round(u_adj_t,3),')',sep=''))
    ## return(list(mean=m,std=v_hat,df_t=df_t,df_adj_t_incorrect=df_adj_t_incorrect,df_adj_t=df_adj_t)) ## AMENDED IN NEXT LINE
    return(list(mean=m,std=v_hat_sqrt,df_t=df_t,df_adj_t_incorrect=df_adj_t_incorrect,df_adj_t=df_adj_t))
  }
}

# 3. catch the warnings (since there are some non-convergence times)
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}


# The simulation
# handle the non-convergence with method 1 (the two methods mentioned in email)

# The parameters
b0 <- 0; b1 <- 1.36; mu_x <- 0
sigma_u <- sqrt(3.37)
sigma_a <- sqrt(0.18)
## sigma_delta <- 0.2
sigma_delta <- sqrt(0.2)
phi <- -1.34
si <- 1
J <- 50;L <- 50
Nimp <- 15 # set totall number of imputations for each incomplete data set


## Set the simulation objects to collect outputs

t_est <- c();t_std <- c()
means <- c();stds <- c()
df_t <- c();df_adj_t_incorrect <- c();df_adj_t <- c()
for(ts in 1:1){ # run 1000 times
  set.seed(ts)
  print(paste('ts',ts))
  # generate one dataset
  arm <- c();X <- c();Y <- c();cluster <- c();
  participant <- c();R <- c()
  for(i in c(0,1)){
    for(j in 1:J){
      a <- rnorm(1,mu_x,sigma_a)
      delta <- rnorm(1,0,sigma_delta)
      for(l in 1:L){
        u <- rnorm(1,0,sigma_u)
        x <- a+u
        pi_ijl <- expit(b0+b1*i+x+delta)
        r_ijl <- expit(phi+si*x)
        y <- rbinom(1,1,pi_ijl) # outcome for each participant
        r <- rbinom(1,1,r_ijl) # missing outcome indicatorfor each participant
        arm <- c(arm,i)
        cluster <- c(cluster,j)
        participant <- c(participant,l)
        X <- c(X,x)
        Y <- c(Y,y)
        R <- c(R,r)
      }
    }
  }
  
  # the data set
  dat <- data.frame(Y=Y,X=X,arm=arm,cluster=cluster,participant=participant,R=R)
  
  # calculate the true effect
  true <- myTryCatch(geeglm(formula=Y~arm,
                       id=cluster , data = dat,
                       family =  binomial("logit"),
                       corstr = "exchangeable"))
  if(is.null(true$warning) == 1 & is.null(true$error)==1){
    # save results
    std0 <- summary(true$value)$coefficients['arm','Std.err']
    beta0 <- summary(true$value)$coefficients['arm','Estimate']
    t_est <- c(t_est,beta0)
    t_std <- c(t_std,std0)
  }else{
    t_est <- c(t_est,NA)
    t_std <- c(t_std,NA)
  }
  
  # generate missing values
  dat$y_mis <- ifelse(dat$R==1,NA,dat$Y)
  
  # MMI
  data.miss <- dat
  y.cat <- data.frame(outcome=data.miss$y_mis)  # data frame for response variables with missing values
  y.numcat <- c(2)                                 # number of levels in outcome variable
  clus <- data.frame(clus=data.miss$cluster)          # data frame for clusters
  nobs <- dim(data.miss)[1]
  x <- data.frame(intercept=rep(1,nobs),covariate=data.miss$X,group=data.miss$arm)
  
  # run to generate Nimp full datasets
  # imp is full datasets 
  imp <- jomo1rancat(Y.cat=y.cat, Y.numcat=y.numcat, X=x,
                    clus=clus,nburn=100, nbetween=25, nimp=Nimp,output=0)
  est <- c();std <- c()
  # Analyze each of the full dataset with GEE
  for(i in 1:Nimp){
    print(c("Imputation",i))
    temp <- imp[imp$Imputation==i,]
    rownames(temp) <- NULL
    print("Dim of imputation outcome");print(length(temp$outcome))
    print("cross-tabulation of outcome and group");print(table(temp$outcome, temp$group))
    temp$outcome <- as.numeric(temp$outcome)-1
    print("cross-tabulation of outcome and group");print(table(temp$outcome, temp$group))
    ex <- myTryCatch(geeglm(formula=outcome~group,
                       id=clus , data = temp,
                       family =  binomial("logit"),
                       corstr = "exchangeable"))
    ## Save the results
    # deal with the non-convergence with method 1
    if(is.null(ex$warning)==1 & is.null(ex$error)==1){
      std1 <- summary(ex$value)$coefficients['group','Std.err']
      beta1 <- summary(ex$value)$coefficients['group','Estimate']
      est <- c(est,beta1)
      std <- c(std,std1)
    }else{
      est <- c(est,NA)
      std <- c(std,NA)
    }
  }
  # pool the results
  t <- mypool(est,std,num=Nimp, print="yes",J=J)
  means <- c(means,t$mean)
  stds <- c(stds,t$std)
  df_t <- c(df_t,t$df_t)
  df_adj_t_incorrect <- c(df_adj_t_incorrect,t$df_adj_t_incorrect)
  df_adj_t <- c(df_adj_t,t$df_adj_t)
}
# do some sanity checks on last stored simulated data set to see that approx 30% missing out of the 5000 participants in the study
sum((data.miss$y_mis == 1)|(data.miss$y_mis == 0),na.rm=TRUE)/5000


results <- list(t_est=t_est,t_std=t_std,means=means,stds=stds,df_t=df_t, df_adj_t_incorrect=df_adj_t_incorrect,df_adj_t=df_adj_t)

#setwd("~/Dropbox/Liz-Melanie/Education IPW Manuscript/Lanqiu simulations/6. Hossain replication analysis - summer 2018/3. Third version  Hossain replication Sept 14th 2018 - with geeM and geepack")
#save(results,file='mimic_geepack_m1_100REPLICATE_TEST.RData') 
save(results,file='mimic_geepack_m1_100REPLICATE_TEST_with_BandWcorrected.RData') 

#load('mimic_geepack_m1_100REPLICATE_TEST.RData')
#load('mimic_geepack_m1_100REPLICATE_TEST.RData') 

# Calculate the coverage
true_value <- mean(results$t_est,na.rm=TRUE)

par(mfrow=c(2,1))
hist(results$t_est)
sd(results$t_est)
hist(results$means)
sd(results$means)

# check how much non-convergence --> doens't look to be any 
sum(results$t_est >0 & results$t_est <100, na.rm=TRUE)

coverage <- sum((results$means-1.96*results$stds)<true_value & (results$means+1.96*results$stds)>true_value,na.rm=TRUE)
coverage_df_t <- sum((results$means-qt(0.975,results$df_t)*results$stds)<true_value & (results$means+qt(0.975,results$df_t)*results$stds)>true_value,na.rm=TRUE)
coverage_df_adj_t_incorrect <- sum((results$means-qt(0.975,results$df_adj_t_incorrect)*results$stds)<true_value & (results$means+qt(0.975,results$df_adj_t_incorrect)*results$stds)>true_value,na.rm=TRUE)
coverage_df_adj_t <- sum((results$means-qt(0.975,results$df_adj_t)*results$stds)<true_value & (results$means+qt(0.975,results$df_adj_t)*results$stds)>true_value,na.rm=TRUE)

cbind(results$means-qt(0.975,results$df_adj_t)*results$stds, true_value, results$means+qt(0.975,results$df_adj_t)*results$stds)

## corrected to take sqrt of v_hat
#coverage <- sum((results$means-1.96*sqrt(results$stds))<true_value & (results$means+1.96*sqrt(results$stds))>true_value,na.rm=TRUE)
#coverage_df_t <- sum((results$means-qt(0.975,results$df_t)*sqrt(results$stds))<true_value & (results$means+qt(0.975,results$df_t)*sqrt(results$stds))>true_value,na.rm=TRUE)
#coverage_df_adj_t_incorrect <- sum((results$means-qt(0.975,results$df_adj_t_incorrect)*sqrt(results$stds))<true_value & (results$means+qt(0.975,sqrt(results$df_adj_t_incorrect))*results$stds)>true_value,na.rm=TRUE)
#coverage_df_adj_t <- sum((results$means-qt(0.975,results$df_adj_t)*sqrt(results$stds))<true_value & (results$means+qt(0.975,results$df_adj_t)*sqrt(results$stds))>true_value,na.rm=TRUE)

# View all estimates of coverage of normal and three variations of t-distribution testing
coverage # 98.5%
coverage_df_t # 98.5%
coverage_df_adj_t_incorrect # 98.9%
coverage_df_adj_t 

# View all estimates of the 3 variations on degrees of freedom for t-distribution testing
mean(results$df_t) # 277.4635
mean(results$df_adj_t_incorrect)
mean(results$df_adj_t) # 20.41984

# over coverage
# there is no convergence problem. 

cbind(results$t_est, results$means, results$df_t, results$df_adj_t_incorrect, results$df_adj_t)














# The simulation
# handle the non-convergence with method 2
t_est <- c();t_std <- c()
means <- c();stds <- c()
df_t <- c();df_adj_t <- c()
for(ts in 1:2){ # run 1000 times
  set.seed(ts)
  print(paste('ts',ts))
  # generate one dataset
  arm <- c();X <- c();Y <- c();cluster <- c();participant <- c();R <- c()
  for(i in c(0,1)){
    for(j in 1:J){
      a <- rnorm(1,mu_x,sigma_a)
      delta <- rnorm(1,0,sigma_delta)
      for(l in 1:L){
        u <- rnorm(1,0,sigma_u)
        x <- a+u
        pi_ijl <- expit(b0+b1*i+x+delta)
        r_ijl <- expit(phi+si*x)
        y <- rbinom(1,1,pi_ijl)
        r <- rbinom(1,1,r_ijl)
        arm <- c(arm,i)
        cluster <- c(cluster,j)
        participant <- c(participant,l)
        X <- c(X,x)
        Y <- c(Y,y)
        R <- c(R,r)
      }
    }
  }
  
  # the data set
  dat <- data.frame(Y=Y,X=X,arm=arm,cluster=cluster,participant=participant,R=R)
  
  # calculate the true effect
  true <- myTryCatch(geeglm(formula=Y~arm,
                       id=cluster , data = dat,
                       family =  binomial("logit"),
                       corstr = "exchangeable"))
  if(is.null(true$warning) == 1 & is.null(true$error)==1){
    # save results
    std0 <- summary(true$value)$coefficients['arm','Std.err']
    beta0 <- summary(true$value)$coefficients['arm','Estimate']
    t_est <- c(t_est,beta0)
    t_std <- c(t_std,std0)
  }else{
    t_est <- c(t_est,NA)
    t_std <- c(t_std,NA)
  }
  
  # generate missing values
  dat$y_mis <- ifelse(dat$R==1,NA,dat$Y)
  
  # MMI
  data.miss <- dat
  y.cat <- data.frame(outcome=data.miss$y_mis)  # data frame for response variables with missing values
  y.numcat <- c(2)                                 # number of levels in outcome variable
  clus <- data.frame(clus=data.miss$cluster)          # data frame for clusters
  nobs <- dim(data.miss)[1]
  x <- data.frame(intercept=rep(1,nobs),covariate=data.miss$X,group=data.miss$arm)
  #Nimp=15
  Nimp <- 15
  # run to generate Nimp full datasets
  # imp is full datasets 
  imp <- jomo1rancat(Y.cat=y.cat, Y.numcat=y.numcat, X=x,
                    clus=clus,nburn=100, nbetween=25, nimp=Nimp,output=0)
  est <- c();std <- c()
  # Analyze each of the full dataset with GEE
  for(i in 1:Nimp){
    print(i)
    temp <- imp[imp$Imputation==i,]
    rownames(temp) <- NULL
    temp$outcome <- as.numeric(temp$outcome)-1
    ex <- myTryCatch(geeglm(formula=outcome~group,
                       id=clus , data = temp,
                       family =  binomial("logit"),
                       corstr = "exchangeable"))
    ## Save the results
    # deal with the non-convergence with method 1
    if(is.null(ex$warning)==1 & is.null(ex$error)==1){
      std1 <- summary(ex$value)$coefficients['group','Std.err']
      beta1 <- summary(ex$value)$coefficients['group','Estimate']
      est <- c(est,beta1)
      std <- c(std,std1)
    }else{
      est <- c(est,NA)
      std <- c(std,NA)
    }
  }
  # pool the results
  t <- mypool(est,std,num=Nimp, print="yes",J=J)
  if(sum(is.na(est))>0){
    means <- c(means,NA)
    stds <- c(stds,NA)
    df_t <- c(df_t,NA) 
    df_adj_t <- c(df_adj_t,NA) 
    print('na')
  }else{
    means <- c(means,t$mean)
    stds <- c(stds,t$std)
    df_t <- c(df_t,t$df_t) 
    df_adj_t <- c(df_adj_t,t$df_adj_t) 
    print(t)
  }
}

results <- list(t_est=t_est,t_std=t_std,means=means,stds=stds,df_t=df_t, df_adj_t=df_adj_t)

save(results,file='mimic_geepack_m2.RData') 

load('mimic_geepack_m2.RData')
# Calculate the coverage
true_value <- mean(results$t_est,na.rm=TRUE)
coverage <- sum((results$means-1.96*results$stds)<true_value & (results$means+1.96*results$stds)>true_value,na.rm=TRUE)
coverage_df_t <- sum((results$means-qt(0.975,results$df_t)*results$stds)<true_value & (results$means+qt(0.975,results$df_t)*results$stds)>true_value,na.rm=TRUE)
coverage_df_adj_t <- sum((results$means-qt(0.975,results$df_adj_t)*results$stds)<true_value & (results$means+qt(0.975,results$df_adj_t)*results$stds)>true_value,na.rm=TRUE)
coverage # 985
coverage_df_t # 985
coverage_df_adj_t # 989
coverage/sum(is.na(results$means)==0) # 98.5%
coverage_df_t/sum(is.na(results$means)==0) # 98.5%
coverage_df_adj_t /sum(is.na(results$means)==0)# 98.9%

mean(results$df_t,na.rm=TRUE) # 277.4635
mean(results$df_adj_t,na.rm=TRUE) # 20.41984




