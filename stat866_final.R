##################################
##### ---- PRELIMINARIES ---- ####
##################################

# Limits printing output
options(max.print = 100)

# Lad in packages
ll <- c('tidyverse','magrittr','cowplot','xtable','stargazer','stringr','forcats','broom',
        'survival','survminer','survMisc','MASS','glmnet','caret')
for (k in ll) {library(k,character.only=T) }

rm(list=ls())

# Set up directory
setwd('C:/Users/erikinwest/Documents/Courses/STAT886/final/')

##############################
##### ---- DATA LOAD ---- ####
##############################

# Color paletter
gg_color_hue <- function(n) { hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n] }

# Note: http://cancerguide.org/trials_glossary.html
# Evaluable Disease: A tumor or tumors which you can tell are present but the size of which cannot be 
# measured accurately. For example, for technical reasons, bone metastases are hard to measure exactly 
# and are usually counted as Evaluable Disease.

# Load in the data
dat <- read.table('ma8_surv.txt',header=T) %>% tbl_df
# Tidy it up a bit 
dat.clean <- 
  dat %>% mutate(arm=fct_recode(factor(arm),'DOX+VNB'='1','DOX'='0'), # arm==1:DOX+VNB, arm==2:DOX
               perform=factor(perform,levels=0:2,
                              labels=c('Active','Restricted','Self-care')), # naturally ordered but probably not right for "linear" effect
               meno=fct_recode(meno,'post'='1','pre'='0'), # pre/post menaupausal
               measure=fct_recode(factor(measure),'yes'='0','no'='1'),
               bone=fct_recode(factor(bone),'yes'='1','no'='0'), # At a bonesite
               nodal=fct_recode(nodal,'pos'='1','neg'='0'),
               er=fct_recode(er,'pos'='1','neg'='0'),
               global=as.numeric(as.character(global)),
               adv=fct_recode(factor(adv),'yes'='1','no'='0'),
               chemmet=fct_recode(factor(chemmet),'yes'='1','no'='0'),
               regimen=fct_recode(factor(regimen),'yes'='1','no'='0')) %>% 
        rename(bonesite=bone,measurable=measure,adverse=adv)


################################################
##### ------------ QUESTION 1  ------------ ####
################################################

# State clearly your results and draw your conclusion based on your results.

# Get the surv object
dat.Surv <- Surv(time=dat$survival,event=dat$dead)

# (a) generate Kaplan-Meier survival curves for patients in each treatment group
dat.survfit <- survfit(dat.Surv~arm,data=dat)
col.survfit <- c("#F8766D","#00BFC4")

# Plot it
plot1a <-
  ggsurvplot(dat.survfit,palette=col.survfit,surv.scale='percent', 
           legend.title='Treatment arm',legend='top',xlab='Time (days)',
           legend.labs=c('DOX','DOX+VNB'),risk.table = T)
           # surv.median.line='hv')
plot1a.plot <- plot1a$plot
plot1a.table <- plot1a$table

# (b) estimate the median survival and associated 95% confidence intervals for each treatment group
# Median is the smallest t such that S<=0.5
tidy.survfit <- tidy(dat.survfit) %>% mutate(treatment=ifelse(strata=='arm=0','DOX','DOX+VNB'))

# Define a function
kernel.ci <- function(tsf,p=0.5,h=0.05,alpha=0.05) {
  tsf.slice <- tsf %>% group_by(treatment) %>% filter(estimate <= (1-p)) %>% filter(time==min(time))
  tp <- tsf.slice$time
  VarSt50 <- tsf.slice$std.error^2
  # The the l/u
  l50 <- tsf %>% group_by(treatment) %>% filter(estimate <= (1-p)-h) %>% filter(time==min(time))
  u50 <- tsf %>% group_by(treatment) %>% filter(estimate >= (1-p)+h) %>% filter(time==max(time))
  # Get the estimate of f.hat
  f.hat <- (u50$estimate-l50$estimate)/(l50$time-u50$time)
  Vart50 <- VarSt50/(f.hat^2)
  lb <- tp-qnorm(0.975)*sqrt(Vart50)
  ub <- tp+qnorm(0.975)*sqrt(Vart50)
  # table
  lbub.tbl <- tibble('Treatment arm'=tsf.slice$treatment,'Lower bound'=lb,'Median'=tp,'Upper bound'=ub)
  return(lbub.tbl)
}

# Save table
ci.tab <- cbind(N=dat.survfit$n,kernel.ci(tsf=tidy.survfit,p=0.5,h=0.05,alpha=0.05)) %>% tbl_df
# Save the avearge lb range
lb.list <- list(half=kernel.ci(tsf=tidy.survfit,p=0.5,h=0.025,alpha=0.05)$`Lower bound` %>% mean,
                curr=kernel.ci(tsf=tidy.survfit,p=0.5,h=0.05,alpha=0.05)$`Lower bound` %>% mean,
                double=kernel.ci(tsf=tidy.survfit,p=0.5,h=0.1,alpha=0.05)$`Lower bound` %>% mean)


# (c) test whether there is any difference in survival between two treatment groups based on log-rank test
# Start writing to an output file
sink('analysis-output.txt')
survMisc::comp(ten(survfit(dat.Surv~arm,data=dat.clean)),p=c(0,1,0,1), q=c(0,0,1,1))
# Stop writing to the file
sink()
# Append to the file
sink('analysis-output.txt', append=TRUE)
sink()
# p=q=0 log-rank test
# p=1,q=0 Mann-Whitney-Wilcoxon test
# Peto-Peto (Andersen) p=0;q=1
# Fleming-Harrington p=1;q=1

# Call in saved output
comp.dat <- read_table('analysis-output.txt') 
comp.clean <-
comp.dat[grep('FH',comp.dat$X1)[1:4],] %>% 
  mutate(Chi2=as.numeric(Z)^2,Pchi=pchisq(Chi2,df=1,lower.tail=F),X1=gsub('FH_','',X1)) %>% 
  rename(Test=X1,Statistic=Chi2,Pvalue=Pchi) %>% dplyr::select(Test,Statistic,Pvalue) %>% 
  mutate(Test=fct_recode(Test,'Log-rank'='p=0_q=0','Wilcoxon'='p=1_q=0',
                         'Peto-Peto'='p=0_q=1','Fleming-Harrington'='p=1_q=1'))


# (d) estimate the hazard ratio (assume it is constant) between two treatment groups and associated 
# 95% confidence interval.
q4.cox <- coxph(dat.Surv~arm,data=dat.clean)
q4.coef <- summary(q4.cox)$coef %>% data.frame
# Calculate the 95% for beta and then exponentiate 
beta.ci <- with(q4.coef,coef + c(-1,1)*qnorm(0.975)*se.coef.)
hr.ci <- exp(beta.ci)
# Make the table 
ci.tbl <- 
tibble(Term=c('$\\hat\\beta$','Hazard ratio ($\\exp\\hat\\beta$)'),'Lower bound'=c(beta.ci[1],hr.ci[1]),
           'Estimate'=c(q4.coef$coef,q4.coef$exp.coef.),'Upper bound'=c(beta.ci[2],hr.ci[2]))


######################################################################
##### ------------ QUESTION 2: EXPLORATORY ANALYSIS  ------------ ####
######################################################################

# For the numeric data
dat.num <- dat.clean %>% dplyr::select_if(is.numeric) %>% dplyr::select(-c(survival,dead)) %>% data.frame
dat.num.names <-  fct_recode(colnames(dat.num),'Age'='age','# of sites'='numsites',
                             'Disease free (days)'='disfree','QLQ'='global','WBC count'='WBC',
                             'Granulocyte count'='gran','Platelet count'='platelet',
                             'Time to diagnosis'='timdiag') %>% as.character()
colnames(dat.num) <- dat.num.names
dat.num.df <-
  dat.num %>% summarise_each(funs(N=length(.)-sum(is.na(.)),Mean=mean(.,na.rm=T),Median=median(.,na.rm=T),
                                  Stdev=sd(.,na.rm=T),Min=min(.,na.rm=T),Max=max(.,na.rm=T))) %>% 
  gather(Variable,Value) %>% separate(Variable,c('Variable','Measure'),'\\_') %>% 
  spread(Measure,Value) %>% dplyr::select(Variable,Mean,Median,Stdev,Min,Max,N)

# For the factor data
dat.fac <-
  dat.clean %>% dplyr::select_if(is.factor) %>% dplyr::select(-c(ID,arm)) %>% 
  mutate(meno=fct_recode(meno,'yes'='post','no'='pre'),
         nodal=fct_recode(nodal,'yes'='pos','no'='neg'),
         er=fct_recode(er,'yes'='pos','no'='neg')) %>% 
  mutate_each(funs(ifelse(as.character(.)=='.',NA,as.character(.)))) %>%
  mutate_each(funs(gsub('yes','Yes',.))) %>%
  mutate_each(funs(gsub('no','No',.))) %>% 
  rename_('Physical ability'='perform','Menopausal'='meno',
          'Measurable tumor'='measurable','Bone tumor'='bonesite',
          'Nodal involvement'='nodal','ER positive'='er',
          'Adverse event'='adverse','Chemotherapy'='chemmet','Prior regimen'='regimen')
# Count it up! 
dat.fac.yn <- dat.fac[,-1] %>% gather(Variable,Value) %>% group_by(Value) %>% 
  count(Variable) %>% na.omit
dat.fac.df <- rbind(dat.fac.yn %>% data.frame,
                    dat.fac.yn %>% group_by(Variable) %>% summarise(n=sum(n)) %>% cbind(Value='Total',.)) %>%
  spread(Value,n) %>% dplyr::select(c(1,4,2,3))
# --- Combine for nicer printing --- #
dat.num.fac.df <- cbind(dat.num.df %>% rename_('Continuous'='Variable'),
      dat.fac.df %>% rename_('Categorical'='Variable')) %>% 
      mutate(N=as.integer(N),Min=paste('(',paste(Min,Max,sep='-'),')',sep='')) %>% 
    rename_('Min-Max'='Min','N '='Total') %>% dplyr::select(-Max)

# For PA variable
pa.var <- dat.fac$`Physical ability` %>% table(deparse.level=0)

# ----- EXPLORATORY DATA ANALYSIS ------ #

# For the continuous variables we want to show the median survival times for features 
#    being/above below the average values
# Mutate each column so it's above/below mean
dat.num.mean <- dat.num %>% tbl_df %>% mutate_each(funs(ifelse(.>mean(.,na.rm=T),T,F))) 
# Store
cont.names <- colnames(dat.num.mean)
cont.n <- length(cont.names)
# Now run the survfit to get the KM estimate of the median for each
dat.num.med <- NULL
  # tibble(Variable=cont.names,'< Mean'=rep(NA,cont.n),'> Mean'=rep(NA,cont.n))
# Loop
for (k in 1:cont.n) {
  vn <- names(dat.num.mean[,k])
  dat.num.med <- 
    rbind(dat.num.med,
  survfit(dat.Surv~vv,data=dat.num.mean[,k] %>% set_colnames('vv')) %>% 
    tidy %>% rename(treatment=strata) %>% kernel.ci(p=0.5,h=0.05,alpha=0.05) %>%
    rename(Measure=`Treatment arm`,LB=`Lower bound`,UB=`Upper bound`) %>% 
    mutate(Measure=ifelse(gsub('vv=','',Measure)=='FALSE','< Mean','> Mean'),Variable=vn) 
    )
}
# Get the median for the whole survival curve
all.median <- survfit(dat.Surv~1) %>% tidy %>% filter(estimate<1/2) %>% filter(time==min(time))
# Clean up the data for better plotting
dat.num.med.gg <- dat.num.med %>% mutate(hack='Continuous variables',diff=(UB-Median)/qnorm(0.975),
                       LB2=Median-diff,UB2=Median+diff)

# gg it
gg.continuous <-
ggplot(dat.num.med.gg,aes(y=Median,x=Variable)) + 
  geom_point(aes(color=Measure),size=3,position=position_dodge(0.2)) + 
  geom_linerange(aes(ymin=LB2,ymax=UB2,color=Measure),linetype=2,position=position_dodge(0.2)) + 
  coord_flip() +
  geom_hline(yintercept = all.median$time,linetype=2) +
  labs(y='Median survival time',caption='SE uses bandwidth h=0.05',
       subtitle='Vertical line shows overall median\nDotted lines show ± SE') + 
  theme(axis.title.y=element_blank()) +
  scale_color_discrete(name='Covariate') + 
  facet_wrap(~hack)

#  ---- Factorial variables ---- #
fac.names <- colnames(dat.fac)
fac.n <- length(fac.names)
# Now run the survfit to get the KM estimate of the median for each
dat.fac.med <- NULL
# tibble(Variable=cont.names,'< Mean'=rep(NA,cont.n),'> Mean'=rep(NA,cont.n))
# Loop
for (k in 2:fac.n) { # Start at 2 to skip Physical activity
  vn <- names(dat.fac[,k])
  dat.fac.med <- 
    rbind(dat.fac.med,survfit(dat.Surv~vv,data=dat.fac[,k] %>% set_colnames('vv')) %>% 
            tidy %>% rename(treatment=strata) %>% kernel.ci(p=0.5,h=0.05,alpha=0.05) %>%
            rename(Measure=`Treatment arm`,LB=`Lower bound`,UB=`Upper bound`) %>% 
            mutate(Measure=gsub('vv=','',Measure),Variable=vn))
}
# Clean up the data for better plotting
dat.fac.med.gg <- dat.fac.med %>% mutate(hack='Categorical variables',diff=(UB-Median)/qnorm(0.975),
                                         LB2=Median-diff,UB2=Median+diff)

# gg it
gg.factorial <-
  ggplot(dat.fac.med.gg,aes(y=Median,x=Variable)) + 
  geom_point(aes(color=Measure),size=3,position=position_dodge(0.2)) + 
  geom_linerange(aes(ymin=LB2,ymax=UB2,color=Measure),linetype=2,position=position_dodge(0.2)) + 
  coord_flip() +
  geom_hline(yintercept = all.median$time,linetype=2) +
  labs(y='Median survival time',subtitle='Vertical line shows overall median\nDotted lines show ± SE') + 
  theme(axis.title.y=element_blank(),legend.position = c(0.65,0.65)) +
  scale_color_manual(name='Covariate',values=c('#00BA38','#C77CFF')) + 
  facet_wrap(~hack)

# Now repeat for the PA variable
pa.dat <- survfit(dat.Surv~PA,data=set_colnames(dat.fac[,'Physical ability'],'PA')) %>% 
  tidy %>% rename(treatment=strata) %>% kernel.ci(p=0.5,h=0.05,alpha=0.05) %>%
  rename(Measure=`Treatment arm`,LB=`Lower bound`,UB=`Upper bound`) %>% 
  mutate(Measure=gsub('PA=','',Measure),Variable=vn) %>%
  mutate(hack='Physical ability',diff=(UB-Median)/qnorm(0.975),
         LB2=Median-diff,UB2=Median+diff)
# gg it
gg.pa <-
ggplot(pa.dat,aes(y=Median,x=Measure)) + 
  geom_point(aes(color=Measure),size=3) + 
  geom_linerange(aes(ymin=LB2,ymax=UB2,color=Measure),linetype=2) + 
  coord_flip() + facet_wrap(~hack) +
  geom_hline(yintercept = all.median$time,linetype=2) +
  labs(y='Median survival time',subtitle='\n') + 
  theme(axis.title.y=element_blank(),legend.position='none') 

######################################################################
##### ------------ QUESTION 2: PARAMETRIC MODELLING  ------------ ####
######################################################################

# ----------------- DROP VARIABLES FOR REGRESSION ANALYSIS ---------------- #

# Three variables have missing values nodal/er/global
# I opt to drop nodal/er and drop rows with missing global values
dat.sub1 <- dat.clean %>% dplyr::select(-c(nodal,er,regimen,arm)) %>%
                filter(!(global=='.' | meno=='.')) %>% # Condence perform into either Self-care not self-care
                mutate(perform=ifelse(perform=='Self-care','Self-care','Active/Restricted'),
                       meno=fct_drop(meno,only='.'))

# dat.sub2 <- dat.clean %>% filter(!(nodal=='.' | er=='.' | global=='.' | meno=='.')) %>%
#                     mutate_if(is.factor,funs(fct_drop(.,only='.'))) %>%
#                     dplyr::select(-c(arm,regimen))

# Get some nice names to match!
clean.sub1 <- 
fct_recode(colnames(dat.sub1[,-c(1:3)]),'Age'='age','# of sites'='numsites', #'Treatment'='arm',
           'Disease free (days)'='disfree','QLQ'='global','WBC count'='WBC',
           'Granulocyte count'='gran','Platelet count'='platelet',
           'Time to diagnosis'='timdiag','Physical ability'='perform','Menopausal'='meno',
           'Measurable tumor'='measurable','Bone tumor'='bonesite',
           'Adverse event'='adverse','Chemotherapy'='chemmet') %>% as.character()

# Create the Surv object
sub1.Surv <- Surv(time=dat.sub1$survival,event=dat.sub1$dead)
# Get the formula names
big.formula <- as.formula(str_c('sub1.Surv',paste(colnames(dat.sub1[,-c(1:3)]),collapse='+'),sep='~'))
# Different dists
dist.list <- c('exponential','weibull','lognormal')

# Fit each model
para.exp <- survreg(big.formula,data=dat.sub1,dist='exponential')
para.wei <- survreg(big.formula,data=dat.sub1,dist='weibull')
para.logn <- survreg(big.formula,data=dat.sub1,dist='lognormal')
# Put in list
para.list <- list(para.exp,para.wei,para.logn)
# Do stepwise selection
# step.aic <- list()
step.bic <- list()
for (k in 1:length(para.list)) {
  # Get the model
  modk <- para.list[[k]]
  # aick <- stepAIC(modk,scope=list(lower=~1,upper=~.),trace=F,direction='backward',k=2)
  bick <- stepAIC(modk,scope=list(lower=~1,upper=~.),trace=F,direction='backward',k=log(nrow(modk$y)))
  step.bic[[k]] <- bick; print(k) #step.aic[[k]] <- aick ; 
}
# Combine the list with the other para
para.master <- list()
for (k in 1:length(para.list)) { 
  para.master[[2*k-1]] <- para.list[[k]]
  para.master[[2*k]] <- step.bic[[k]] 
}


# Custom stats
dof.stat <- c('D.o.F',unlist(lapply(para.master,function(ll) ll$df.residual)))
dev.stat <- c('Deviance',round(unlist(lapply(para.master,function(ll) 2*(ll$loglik[2]-ll$loglik[1]))),1))
crit.stat <- c('Crit. Dev.',rep(qchisq(0.95,df=1),length(para.master)) %>% round(1))
scale.stat <- c('Scale',round(unlist(lapply(para.master,function(ll) ll$scale)),2))
lines.stat <- list(dof.stat,dev.stat,crit.stat,scale.stat)

# See if we can reject the null of no model difference for stepwise vs non-stepwise
nmod <- length(para.list)
distname <- c('Exponential','Weibull','Lognormal')
df.lrtest <- tibble(Distribution=distname,'$\\Delta p$'=rep(NA,nmod),
                    Deviance=rep(NA,nmod),'Critical Value'=rep(NA,nmod))
# loop it
for (k in 1:nmod) { 
  print(k)
  # Get difference in p
  dp <- step.bic[[k]]$df.residual-para.list[[k]]$df.residual
  # And critical value
  crit <- qchisq(0.95,dp,lower.tail = T)
  # Diff if LL
  dLL <- 2*(para.list[[k]]$loglik[2]-step.bic[[k]]$loglik[2])
  # Store
  df.lrtest[k,2:4] <- c(dp,dLL,crit)
}
df.lrtest <- df.lrtest %>% mutate('Reject Null?'=ifelse(Deviance>`Critical Value`,'Yes','No'),
                                     `$\\Delta p$`=as.integer(`$\\Delta p$`))

# Define the winning candidate models
cand.mdls <- list(Exponential=para.list[[1]],Weibull=step.bic[[2]],Lognormal=step.bic[[3]])

# ----- RESIDUAL ANALYSIS ------ #

# Define the cumulative hazard functions for each
xb.exp <- cand.mdls$Exponential$linear.predictors
xb.wei <- cand.mdls$Weibull$linear.predictors
xb.logn <- cand.mdls$Lognormal$linear.predictors
# Get the scales for the two models
scale.wei <- cand.mdls$Weibull$scale
scale.logn <- cand.mdls$Lognormal$scale
# Get the event and time
sub1.time <- dat.sub1$survival
sub1.delta <- dat.sub1$dead

# Calculate H(t) for each
sub1.Ht <- tibble(res.exp=sub1.time*exp(-xb.exp),
                   res.weibull=exp( (log(sub1.time) - xb.wei)/scale.wei ),
                   res.logn=-log(1-pnorm((log(sub1.time) - xb.logn)/scale.logn)),
                   time=sub1.time,
                   event=sub1.delta)
# Put in the line form
sub1.Ht <- sub1.Ht %>% gather(dist,res,-time,-event) %>%
  group_by(dist) %>% # mutate(res.cw=res + (1-event)) %>% 
  arrange(dist,res) %>% mutate(idx=str_c(dist,1:length(dist),sep='.'))
# Estimate a KM survial reslationship between the residuals and the censoring events
sub1.rHt <-
  sub1.Ht %>% group_by(dist) %>% 
  do(SO=survfit(Surv(time=.$res,event=.$event)~1)) %>% tidy(SO) %>% 
  dplyr::select(c(dist,time,estimate)) %>% 
  mutate(Ht=-log(estimate),idx=str_c(dist,1:length(dist),sep='.')) %>% 
  rename(ri=time,St=estimate) %>% ungroup %>% dplyr::select(-dist)
# Combine
sub1.res <- sub1.Ht %>% left_join(sub1.rHt,by='idx') %>% ungroup() %>% 
                  mutate(dist=lvls_reorder(dist,c(1,3,2)))

# Plot ri against -logSt(ri)=Ht(ri)
res.labeller <- c('res.exp'='Exponential','res.weibull'='Weibull','res.logn'='Lognormal')
gg.coxsnell <- ggplot(sub1.res,aes(x=res,y=Ht)) + 
  geom_point(aes(color=dist),show.legend = F) + 
  facet_wrap(~dist,nrow=1,labeller=labeller(dist=res.labeller)) + 
  geom_abline(intercept=0,slope=1,color='black',linetype=2) + 
  labs(x=expression(r[i]),y=expression(hat(H) * '(' * r[i] * ')'))

# ----- TESTING THE WEIBULL AFT ASSUMPTIONS ------- #

# Save data for later
dat.wei <- with(dat.sub1,survfit(Surv(time=survival,event=dead)~1)) %>% tidy %>% 
  mutate(logt=log(time),logHt=log(-log(estimate)))
# Get the linear regression model
mdl.wei <- lm(logHt~logt,data=dat.wei) %>% coef %>% round(1)
str.mdl.wei <- paste(mdl.wei[1],'+',mdl.wei[2],'* log(t)',sep=' ')
str.mdl.df <- data.frame(x=4,y=0.5,label=paste('Regression formula',str.mdl.wei,sep='\n'))
# gg it
gg.km.wei <-
ggplot(dat.wei,aes(x=logt,y=logHt)) + 
  geom_point(size=1) + 
  geom_smooth(method='lm',se=F,color='red') + 
  geom_text(data=str.mdl.df,aes(x=x,y=y,label=label),color='red',size=5) + 
  labs(x=expression(log(t)),y=expression(log(H(t))),subtitle='From KM estimates')

# ---- CHECK PARALLEL LINES FOR ONE OF THE FEATURES ----- #

dat.numsites <- dat.sub1 %>% dplyr::select(survival,dead,numsites) %>% 
  mutate(fac.numsites=fct_lump(factor(numsites),n=4)) %>% 
  mutate(fac.numsites=fct_recode(fac.numsites,'5+'='Other'))

# Update
dat.surv.numsites <- 
dat.numsites %>% group_by(fac.numsites) %>% 
  do(SO=survfit(Surv(time=.$survival,event=.$dead)~1)) %>% 
  tidy(SO) %>% mutate(logt=log(time),logHt=log(-log(estimate)))

# Plot it
gg.parallel <-
ggplot(dat.surv.numsites,aes(x=logt,y=logHt)) + 
  geom_point(aes(color=fac.numsites),size=1) + 
  scale_color_discrete(name='# of sites') + 
  labs(x=expression(log(t)),y=expression(log(H(t))),
       subtitle='From KM estimates') + 
  theme(legend.position = c(0.2,0.7))

###########################################################################
##### ------------ QUESTION 3: SEMI-PARAMETRIC MODELLING  ------------ ####
###########################################################################

# Run the full Cox model (note we need to recode perform for naming issues)
cox.full <- survival::coxph(big.formula,data=dat.sub1 %>% 
                              mutate(perform=ifelse(perform=='Self-care','Selfcare',perform)))
# Run the backward selection Cox with BIC
cox.bw <- stepAIC(cox.full,scope=list(lower=~1,upper=~.),trace=F,direction='backward',k=log(nrow(cox.full$y)))
# Run the COX Lasso with CV to find the best "lambda"
y.cox <- dat.sub1[,c('survival','dead')] %>% set_colnames(c('time','status')) %>% as.matrix
x.cox <- model.matrix(~.,data=dat.sub1[,-c(1:3)])[,-1]
# Use cross validation
cvlambda.range.cox <- cv.glmnet(x=x.cox,y=y.cox,family='cox',alpha=1,nlambda = 500,lambda.min.ratio = 0.00001)
lambda.1se <- cvlambda.range.cox$lambda.1se
lambda.min <- cvlambda.range.cox$lambda.min
# Get the CV 
cox.lasso <- glmnet(x=x.cox,y=y.cox,family='cox',alpha=1,lambda=lambda.1se)

# Create a dummy variable version for latex
dummy.X <- cbind(y=dat.sub1$survival,model.matrix(~.,data=dat.sub1[,-c(1:3)])[,-1]) %>% data.frame %>%
    set_colnames(gsub('\\.','',names(.)))
dummy.lm <- lm(y~-1+.,data=dummy.X)
# Now get the coefficient and standard error lists
coef.lasso <- coef(cox.lasso) %>% as.numeric
names(coef.lasso) <- coef(cox.lasso) %>% rownames
se.lasso <- rep(NA,length(coef.lasso))
names(se.lasso) <- coef(cox.lasso) %>% rownames
# Store in a list
coef.list <- list(coef(cox.full),coef(cox.bw),coef.lasso[which(coef.lasso!=0)])
se.list <- list(summary(cox.full)$coef[,3],summary(cox.bw)$coef[,3],se.lasso[which(coef.lasso!=0)])
# Now create the add lines list
add.lines.cox <- list(obs=c('N',rep(length(dummy.lm$residuals),3)),
                      p=c('\\# of params',list(cox.full,cox.bw,cox.lasso) %>% 
                            lapply(.,function(ll) ll %>% coef %>% equals(0) %>% `n'est pas`() %>% sum) %>% unlist),
                      LL=c('LogLik',round(unlist(lapply(list(cox.full,cox.bw),function(ll) ll$loglik[2])),1),NA))


# ----- MODEL CHECKING ----- #

# Have a cleaning dictionary 
clean.dict <- 
c('Age'='age','# of sites'='numsites',
'Disease free (days)'='disfree','QLQ'='global','WBC count'='WBC',
'Granulocyte count'='gran','Platelet count'='platelet',
'Time to diagnosis'='timdiag','Physical ability'='perform','Menopausal'='meno',
'Measurable tumor'='measurable','Bone tumor'='bonesite',
'Adverse event'='adverse','Chemotherapy'='chemmet')

# Use the Schoenfeld residuals and test developed by Grambsch and Therneau
test.zph <- cox.zph(cox.bw,transform='km')

# Put the residuals and time in a data.frame
schoen.dat <- data.frame(time=rownames(test.zph$y),tt=test.zph$x,schoen=test.zph$y) %>% tbl_df %>%
  mutate(time=time %>% as.character %>% as.numeric) %>% arrange(time) %>% 
    gather(var,val,-time,-tt) %>% separate(var,c('drop','var'),'[.]') %>% dplyr::select(-drop)
# Clean up the variables..
clean.nam <- sapply(unique(schoen.dat$var),function(ss) grep(str_sub(ss,1,5),clean.dict,value=T)) %>%
  names %>% str_split_fixed(.,'[.]',2)  %>% extract(,2)
dirty.nam <- unique(schoen.dat$var)
# Update
schoen.dat$var <- factor(schoen.dat$var,levels=dirty.nam,labels=clean.nam)

# Plot it
gg.schoen <- ggplot(schoen.dat,aes(x=time,y=val)) + 
  geom_point(aes(color=var),show.legend = F,size=1) + 
  facet_wrap(~var,nrow=2,scales='free_y') + 
  labs(x='Time (days)',y='Schoenfeld residuals',subtitle='Black line shows LOESS smoother') + 
  geom_smooth(se=F,method = 'loess',color='black',linetype=2)

# Make a pretty table
pretty.zph <-
test.zph$table %>% data.frame %>% rownames_to_column('Variable') %>% 
  mutate(Variable=factor(Variable,levels=c(dirty.nam,'GLOBAL'),labels=c(clean.nam,'Global test'))) %>% 
  rename(Statistic=chisq,Pvalue=p) %>% dplyr::select(-rho)

######################################
##### ---- SAVE FOR MARKDOWN ---- ####
######################################

# Put in list
rmd.list <- list(dat.num.fac.df=dat.num.fac.df,
                 plot1a.plot=plot1a.plot,
                 plot1a.table=plot1a.table,
                 pa.var=pa.var,
                 ci.tab=ci.tab,
                 lb.list=lb.list,
                 comp.clean=comp.clean,
                 ci.tbl=ci.tbl,
                 all.median=all.median,
                 gg.continuous=gg.continuous,
                 gg.factorial=gg.factorial,
                 gg.pa=gg.pa,
                 para.master=para.master,
                 lines.stat=lines.stat,
                 clean.sub1=clean.sub1,
                 df.lrtest=df.lrtest,
                 gg.coxsnell=gg.coxsnell,
                 str.mdl.df=str.mdl.df,
                 gg.km.wei=gg.km.wei,
                 gg.parallel=gg.parallel,
                 dummy.lm=dummy.lm,
                 coef.list=coef.list,
                 se.list=se.list,
                 add.lines.cox=add.lines.cox,
                 gg.schoen=gg.schoen,
                 pretty.zph=pretty.zph)
# Save
save(rmd.list,file = 'rmd_data.RData')
