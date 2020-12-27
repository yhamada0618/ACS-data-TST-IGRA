
d<-read.delim("../acs_clinical_DB.txt", header = T, sep = ",", dec = ".")

#manipulate to match variable names and coding with Rishi's data dictionary

d[d=="999"|d=="998"]<-NA
d[d==""]<-NA
library(tidyr)
library(dplyr)
d<-as_tibble(d)
d<-filter(d,!is.na(SubjectID))
#correct error in QFT result
d$QFT_TBAg_Nil[d$QFT_TBAg_Nil=="positive"]<-12.6
d$QuantiferonInTubeResult[d$QuantiferonInTubeResult=="12.6"]<-"positive"
d$QFT_TBAg_Nil<-as.numeric(d$QFT_TBAg_Nil)
d<-rename(d,studyid=SubjectID,
          age=AgeAtLastBirthDay,
          contact=LivedInHouseWithTBPatient,
          prevbcg=VaccinatedWithBCG,
          prevtbdiag=PreviousDiagnosisOfTB,
          prevtbdiagyear=YearPreviousDiagnosisTB,
          ethinicity=Ethnicity,
          hivpos=HIVPositive,
          height=Height,
          weight=Weight,
          recruitmentdate=EnrolmentDate,
          censordate=DateCensor,
          studytime=DateSurvExit,
          dateofdeath=DateDeath,
          ltfu_date=DateAWOL,
          mantoux_date=DateTSTRead,
          mantoux_result=TSTResult,
          qfngit_date=DateOfQuantiferon,
          igrapos=QuantiferonInTubeResult,
          qfngit_tbag_nil=QFT_TBAg_Nil,
          sex=Sex)
d$mantoux_result<-as.numeric(d$mantoux_result)
#study ID 40371 has invalid TST result. No numeric value"
d[,c("studytime","dateofdeath","ltfu_date","censordate","recruitmentdate","qfngit_date","mantoux_date")]<-lapply(d[,c("studytime","dateofdeath","ltfu_date","censordate","recruitmentdate","qfngit_date","mantoux_date")],as.Date.character)
d$igrapos<-as.factor(d$igrapos)
d$igrapos<-factor(d$igrapos,levels=c("positive","negative"),labels = c("positive","negative"))          

d<-mutate(d,diabetes=ifelse(ChronicIllnessDetails=="diabetes",1,0),
          bmi=weight/(height^2),
          died=ifelse(!is.na(dateofdeath),1,0),
          ltfu=ifelse(!is.na(ltfu_date),1,0),
          activetb=ifelse(TBCaseType=="prevalent"|is.na(TBCaseType),0,1),
          activetbdate=ifelse(activetb==1,studytime,NA),
          mantoux_bin=ifelse(mantoux_result>=5,1,0),
          qfngit_result=igrapos)
d$sex<-as.factor(d$sex)
d$sex<-factor(d$sex,levels=c("male","female"))
d$contact<-as.factor(d$contact)
d$BCGScar<-as.factor(d$BCGScar)
d$prevtbdiag<-as.factor(d$prevtbdiag)
d$hivpos<-as.factor(d$hivpos)
d$igrapos<-as.factor(d$igrapos)
d$igrapos<-factor(d$igrapos,levels=c("negative","positive"))
d$mantoux_bin<-as.factor(d$mantoux_bin)
d$mantoux_bin<-factor(d$mantoux_bin,levels=c("0","1"),labels=c("negative","positive"))


#estimate TB incidence. Include only baseline visit
d2<-filter(d,VisitType=="D0")
#exclude participants with prevalent TB and one HIV+ participant
d2<-filter(d2,TBCaseType=="incident"|is.na(TBCaseType))
d2<-filter(d2,hivpos!=1)
d2<-mutate(d2,followup=as.numeric(studytime-recruitmentdate))

#exclude participants without followup
d2<-filter(d2,followup>0)

#overall incidence rate per 100p-y
sum(d2$activetb==1)/(sum(d2$followup)/365)*100
inc<-d2 %>% select(followup,mantoux_bin,activetb,igrapos) %>% group_by(mantoux_bin,igrapos) %>% summarise_at(c("activetb","followup"),sum,na.rm=T)
inc_tst<-inc %>% group_by(mantoux_bin) %>% summarise_at(c("activetb","followup"),sum, na.rm=T)
inc_tst<-mutate(inc_tst,ir=activetb/(followup/365)*100)
inc_qft<-inc %>% group_by(igrapos) %>% summarise_at(c("activetb","followup"),sum, na.rm=T)
inc_qft<-mutate(inc_qft,ir=activetb/(followup/365)*100)

library(gtsummary)
library("survival")

#binary cox

library("survminer")
fit<-coxph(Surv(followup, activetb) ~ mantoux_bin, data = d2)
tbl_regression(fit,exponentiate = T)

fit<-coxph(Surv(followup, activetb) ~ igrapos, data = d2)
tbl_regression(fit,exponentiate = T)

#explore different TST cut-off values

d2<-mutate(d2,tst10=ifelse(mantoux_result>=10,1,0),
           tst15=ifelse(mantoux_result>=15,1,0),
           tst_stratify=ifelse(mantoux_result>=15&prevbcg=="yes",1,
                               ifelse(mantoux_result>=5&prevbcg=="no",1,0)
                               )
           )
fit<-coxph(Surv(followup, activetb) ~ tst10, data = d2)
tbl_regression(fit,exponentiate = T)
fit<-coxph(Surv(followup, activetb) ~ tst15, data = d2)
tbl_regression(fit,exponentiate = T)
fit<-coxph(Surv(followup, activetb) ~ tst_stratify, data = d2)
tbl_regression(fit,exponentiate = T)

#binary poisson

fit<-glm(activetb ~ mantoux_bin+offset(log(followup)),family = poisson(link = "log"), data = d2)
tbl_regression(fit,exponentiate = T)

fit<-glm(activetb ~ igrapos+offset(log(followup)),family = poisson(link = "log"), data = d2)
tbl_regression(fit,exponentiate = T)

fit<-glm(activetb ~ tst10+offset(log(followup)),family = poisson(link = "log"), data = d2)
tbl_regression(fit,exponentiate = T)
fit<-glm(activetb ~ tst15+offset(log(followup)),family = poisson(link = "log"), data = d2)
tbl_regression(fit,exponentiate = T)
fit<-glm(activetb ~ tst_stratify+offset(log(followup)),family = poisson(link = "log"), data = d2)
tbl_regression(fit,exponentiate = T)

#use quantitative value
#normalize quantiative QFT and tst results
library(BBmisc)
#d2<-mutate(d2,normqft=cut(qfngit_tbag_nil,breaks=c(-Inf,unique(quantile(qfngit_tbag_nil,seq(0,1,0.01),na.rm=T))
  #                                               )
 #                        )
#)
#d2$normqft<-as.numeric(d2$normqft)
d2<-mutate(d2,normqft=normalize(qfngit_tbag_nil,method="range",range=c(0,1)))

d2<-mutate(d2,normtst=normalize(mantoux_result,method="range",range=c(0,1)))


#d2$normqft<-as.numeric(d2$normqft)
#d2<-mutate(d2,normtst=cut(qfngit_tbag_nil,breaks=c(-Inf,unique(quantile(mantoux_result,seq(0,1,0.01),na.rm=T))
 #                                                )
  #                      )
   #    )
#d2$normtst<-as.numeric(d2$normtst)

#fit without normalization

library(rms)
fit_tst<-coxph(Surv(followup, activetb) ~ rcs(mantoux_result,3), data = d2)
d2$pred_tst[!is.na(d2$mantoux_result)] <- predict(fit_tst, data = d2, type="risk")
fit_qft<-coxph(Surv(followup, activetb) ~ rcs(qfngit_tbag_nil,3), data = d2)
d2$pred_qft[!is.na(d2$qfngit_tbag_nil)] <- predict(fit_qft, data = d2, type="risk")
ggplot(d2) +
  geom_smooth(aes(x=qfngit_tbag_nil, y=pred_qft)) 
ggplot(d2) +
  geom_smooth(aes(x=mantoux_result, y=pred_tst)) 


#prediction using normalized data

if(FALSE){
fit_tst<-coxph(Surv(followup, activetb) ~ rcs(normtst,3), data = d2)
d2$pred_tst[!is.na(d2$mantoux_result)] <- predict(fit_tst, data = d2, type="risk")
fit_qft<-coxph(Surv(followup, activetb) ~ rcs(normqft,3), data = d2)
d2$pred_qft[!is.na(d2$normqft)] <- predict(fit_qft, data = d2, type="risk")
d3<-select(d2,pred_tst,pred_qft,normtst,normqft)
d3<-gather(d3, test, norm_value, normtst:normqft)
d3<-mutate(d3,test=ifelse(test=="normtst","TST","IGRA"))
d3<-mutate(d3,pred=ifelse(test=="TST",pred_tst,pred_qft))
d3<-select(d3,norm_value,test,pred)
library(ggplot2)
ggplot(d3) +
geom_smooth(aes(x=norm_value, y=pred,color=test)) 
}

#sub-group-by contact


fit<-glm(activetb ~ mantoux_bin*contact+offset(log(followup)),family = poisson(link = "log"), data = d2)
tbl_regression(fit,exponentiate = T)

fit<-glm(activetb ~ igrapos*contact+offset(log(followup)),family = poisson(link = "log"), data = d2)
tbl_regression(fit,exponentiate = T)

#predictors
d4<-select(d2,activetb,tst15,igrapos,contact,age,sex,prevtbdiag,prevbcg,bmi,followup)
gg_miss_var(d4)
res<-summary(aggr(d4, sortVar=TRUE))$combinations
d4<-d4[complete.cases(d4),]
library(glmnet)
full<- glm(activetb~tst15+igrapos+contact+age+sex+prevtbdiag+prevbcg+bmi+offset(log(followup)), d4,family = poisson(link="log"))
step.model <- stepAIC(full, direction = "both", 
                      trace = FALSE)
summary(step.model)
