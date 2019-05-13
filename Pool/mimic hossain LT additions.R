#######################################################################
#  * Mimic one simulation scenario in Hossain's paper
#    * scenario 1 in the paper
#  * k=50; M=50; number of clusters equals to 50 and number of participants
#    in each clusters eaquals to 50.
#  * cluster sizes are the same
#######################################################################

# Install packages - added by Liz
install.packages("jomo")
install.packages("geeM")
install.packages("geepack")

# Load in library
library(jomo)
library(geeM)
library(geepack)

# Functions
## 1. Expit function
expit=function(x){y=exp(x)/(1+exp(x));return(y)}

## 2. Pool function
mypool=function(mean0,sd0,num=5,print='no',J=50){
  m=mean(mean0,na.rm=TRUE)
  v=mean(sd0,na.rm=TRUE)
  B=sd(mean0,na.rm=TRUE)
  v_hat=v+(1+1/num)*B
  l=m-1.96*v_hat
  u=m+1.96*v_hat
  
  # Liz additions with df of t distribution for testing based on standard results from MI literature
  df_denom <- (1+1/num)*B
  df_part <- 1+v/df_denom
  df_t <- (num-1)*df_part^2 # LIZ EDIT to calculate the df of the t distribution we should compare this to
  l_t=m-qt(0.975,df_t)*v_hat
  u_t=m+qt(0.975,df_t)*v_hat
  
  
  # Liz additions with df of t distribution for testing based on results in Hossain re adjustment for MMI feature
  # i.e. need to reduce the df based on a formula from Barnard and Rubin (1999)
  df_com <- 2*J - 2 #df for full data - calculated this way even have same # clusters in each arm
  # note that J is used to denote # clusters in each arm
  df_obs <- df_com*((df_com+1)/(df_com+3))*(1/df_part)
  df_adj_t <- 1/(1/df_t + 1/df_obs) # see formula on pg. 7 of Hossain 2017
    
  print("Standard t df and Barnard/Rubin adjusted t df");print(c(df_t, df_adj_t))
  print("97.5% quantiles from standard t df and Barnard/Rubin adjusted t df");print(c(qt(0.975,df_t), qt(0.975,df_adj_t)))
  l_t=m-qt(0.975,df_t)*v_hat
  u_t=m+qt(0.975,df_t)*v_hat
  l_adj_t=m-qt(0.975,df_adj_t)*v_hat
  u_adj_t=m+qt(0.975,df_adj_t)*v_hat
  
  if(print=='no'){
    return(list(mean=m,std=v_hat))
  }
  if(print=='yes'){
    print('mean (95% CI)')
    print(paste(round(m,3)," (",round(l,3),',',round(u,3),')',sep=''))
    print(paste(round(m,3)," (",round(l_t,3),',',round(u_t,3),')',sep=''))
    print(paste(round(m,3)," (",round(l_adj_t,3),',',round(u_adj_t,3),')',sep=''))
    return(list(mean=m,std=v_hat,df_t=df_t,df_adj_t=df_adj_t))
  }
}

# The parameters
b0=0;b1=1.36
mu_x=0
sigma_u=sqrt(3.37)
sigma_a=sqrt(0.18)
sigma_delta=0.2
phi=-1.34
si=1
J=50;L=50

# The simulation
t_est=c();t_std=c()
means=c();stds=c()
df_t=c();df_adj_t=c()
#for(ts in 1:1000){ # run 1000 times
for(ts in 1:100){ # run 1000 times
  set.seed(ts)
  print(paste('ts',ts))
  # generate one dataset
  arm=c();X=c();Y=c();cluster=c();participant=c();R=c()
  for(i in c(0,1)){
    for(j in 1:J){
      a=rnorm(1,mu_x,sigma_a)
      delta=rnorm(1,0,sigma_delta)
      for(l in 1:L){
        u=rnorm(1,0,sigma_u)
        x=a+u
        pi_ijl=expit(b0+b1*i+x+delta)
        r_ijl=expit(phi+si*x)
        y=rbinom(1,1,pi_ijl)
        r=rbinom(1,1,r_ijl)
        arm=c(arm,i)
        cluster=c(cluster,j)
        participant=c(participant,l)
        X=c(X,x)
        Y=c(Y,y)
        R=c(R,r)
      }
    }
  }
  # the data set
  dat=data.frame(Y=Y,X=X,arm=arm,cluster=cluster,participant=participant,R=R)
  
  # calculate the true effect
  true=geem(formula=Y~arm,
            id=cluster , data = dat,
            family =  binomial("logit"),
            corstr = "exchangeable")
  # save results
  std0=summary(true)$se.robust[2]
  beta0=true$beta
  t0=beta0[2]
  t_est=c(t_est,t0)
  t_std=c(t_std,std0)
  
  # generate missing values
  dat$y_mis=ifelse(dat$R==1,NA,dat$Y)
  
  # MMI
  data.miss=dat
  y.cat= data.frame(outcome=data.miss$y_mis)  # data frame for response variables with missing values
  y.numcat=c(2)                                 # number of levels in outcome variable
  clus=data.frame(clus=data.miss$cluster)          # data frame for clusters
  nobs=dim(data.miss)[1]
  x= data.frame(intercept=rep(1,nobs),covariate=data.miss$X,group=data.miss$arm)
  #Nimp=15
  Nimp=15
  # run to generate Nimp full datasets
  # imp is full datasets 
  imp = jomo1rancat(Y.cat=y.cat, Y.numcat=y.numcat, X=x,
                    clus=clus,nburn=100, nbetween=25, nimp=Nimp,output=0)
  est=c();std=c()
  # Analyze each of the full dataset with GEE
  for(i in 1:Nimp){
    print(i)
    temp=imp[imp$Imputation==i,]
    rownames(temp)=NULL
    temp$outcome=as.numeric(temp$outcome)-1
    ex=geem(formula=outcome~group,
            id=clus , data = temp,
            family =  binomial("logit"),
            corstr = "exchangeable")
    ## Save the results
    std1=summary(ex)$se.robust[2]
    beta1=ex$beta
    t1=beta1[2]
    est=c(est,t1)
    std=c(std,std1)
  }
  # pool the results
  t=mypool(est,std,num=Nimp, print="yes",J=J)
  means=c(means,t$mean)
  stds=c(stds,t$std)
  df_t =c(df_t,t$df_t) #added by Liz
  df_adj_t =c(df_adj_t,t$df_adj_t) #added by Liz
  print(t)
}

results=list(t_est=t_est,t_std=t_std,means=means,stds=stds,df_t=df_t, df_adj_t=df_adj_t)



save(results,file='mimic.100replicates.09.04.2018.RData') 

load('mimic.100replicates.09.04.2018.RData')

# Calculate the coverage
true_value=mean(t_est)
coverage=sum((results$means-1.96*results$stds)<true_value & (results$means+1.96*results$stds)>true_value)
coverage_df_t=sum((results$means-qt(0.975,df_t)*results$stds)<true_value & (results$means+qt(0.975,df_t)*results$stds)>true_value)
coverage_df_adj_t=sum((results$means-qt(0.975,df_adj_t)*results$stds)<true_value & (results$means+qt(0.975,df_adj_t)*results$stds)>true_value)
coverage # 98.2%cd
coverage_df_t
coverage_df_adj_t
