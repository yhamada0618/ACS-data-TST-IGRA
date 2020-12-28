
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
d$activetb<-as.factor(d$activetb)
d$activetb<-factor(d$activetb, levels=c(0,1),labels=c("no","yes"))
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

library(table1)
table1(~activetb+mantoux_bin+mantoux_result+igrapos+qfngit_tbag_nil+contact+age+sex+prevtbdiag+prevbcg+bmi+MotherGrossMonthlyIncome +FatherGrossMonthlyIncome+GuardianGrossMonthlyIncome,data=d2)

#overall incidence rate per 100p-y
d2$activetb<-as.numeric(d2$activetb)
d2$activetb<-d2$activetb-1
sum(d2$activetb==1)/(sum(d2$followup)/365)*100
inc<-d2 %>% select(followup,mantoux_bin,activetb,igrapos) %>% group_by(mantoux_bin,igrapos) %>% summarise_at(c("activetb","followup"),sum,na.rm=T)
inc_tst<-inc %>% group_by(mantoux_bin) %>% summarise_at(c("activetb","followup"),sum, na.rm=T)
inc_tst<-mutate(inc_tst,ir=activetb/(followup/365)*100)
inc_qft<-inc %>% group_by(igrapos) %>% summarise_at(c("activetb","followup"),sum, na.rm=T)
inc_qft<-mutate(inc_qft,ir=activetb/(followup/365)*100)



#binary cox
library("survival")
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

#fit without normalization

#replace 0< with 0 for QFT results
d2<-mutate(d2,qfngit_tbag_nil=ifelse(qfngit_tbag_nil<0,0,qfngit_tbag_nil))

library(rms)
fit_tst<-coxph(Surv(followup, activetb) ~ rcs(mantoux_result,3), data = d2)
d2$pred_tst[!is.na(d2$mantoux_result)] <- predict(fit_tst, data = d2, type="risk")
fit_qft<-coxph(Surv(followup, activetb) ~ rcs(qfngit_tbag_nil,3), data = d2)
d2$pred_qft[!is.na(d2$qfngit_tbag_nil)] <- predict(fit_qft, data = d2, type="risk")
ggplot(d2) +
  geom_smooth(aes(x=qfngit_tbag_nil, y=pred_qft)) 
ggplot(d2) +
  geom_smooth(aes(x=mantoux_result, y=pred_tst)) 

#normalize quantiative QFT and tst results
library(BBmisc)
d2<-mutate(d2,mantoux_result=round(mantoux_result,0))
qft<-quantile(d2$qfngit_tbag_nil,seq(0.01,1,0.01),na.rm=T)
tst<-quantile(d2$mantoux_result,seq(0.01,1,0.01),na.rm=T)
tile<-as.data.frame(cbind(qft,tst))
tile$tst<-round(tile$tst)
tile_tst<-data.frame(x=cumsum(table(tile$tst)),y=table(tile$tst))
tile_tst<-mutate(tile_tst,p=ifelse(y.Freq>1,x-y.Freq+1,x))
tile_tst$p[1]<-0
tile<-mutate(tile,p_tile_tst=tile_tst$p[match(tst,tile_tst$y.Var1)])
tile_qft<-data.frame(x=cumsum(table(tile$qft)),y=table(tile$qft))
tile_qft<-mutate(tile_qft,p=ifelse(y.Freq>1,x-y.Freq+1,x))
tile_qft$p[1]<-0
tile_qft$y.Var1<-as.numeric(as.character(tile_qft$y.Var1))
tile<-mutate(tile,p_tile_qft=tile_qft$p[match(round(qft,2),round(tile_qft$y.Var1,2))])


d2<-mutate(d2,normqft=cut(qfngit_tbag_nil,breaks=c(-Inf,unique(quantile(qfngit_tbag_nil,seq(0,1,0.01),na.rm=T))
                                             )
                    )
)
tst<-distinct(tile,tst,p_tile_tst)
qft<-distinct(tile,qft,p_tile_qft)



d2<-mutate(d2,normtst=cut(mantoux_result,breaks=c(-Inf,unique(round(quantile(mantoux_result,seq(0,1,0.01),na.rm=T),0)
)
)
,labels=tst$p_tile_tst)
)

d2<-mutate(d2,normqft=cut(qfngit_tbag_nil,breaks=c(-Inf,unique(quantile(qfngit_tbag_nil,seq(0,1,0.01),na.rm=T))
                                                )
                      ,labels=qft$p_tile_qft)
    )
d2$normtst<-as.numeric(as.character(d2$normtst))
d2$normqft<-as.numeric(as.character(d2$normqft))

#d2<-mutate(d2,normqft=normalize(qfngit_tbag_nil,method="range",range=c(0,1)))

#d2<-mutate(d2,normtst=normalize(mantoux_result,method="range",range=c(0,1)))


#prediction using normalized data


fit_tst<-coxph(Surv(followup, activetb) ~ rcs(normtst,c(10,50,90)), data = d2)
d2$pred_tst[!is.na(d2$mantoux_result)] <- predict(fit_tst, data = d2, type="risk")
fit_qft<-coxph(Surv(followup, activetb) ~ rcs(normqft,c(10,50,90)), data = d2)
d2$pred_qft[!is.na(d2$normqft)] <- predict(fit_qft, data = d2, type="risk")
d3<-select(d2,pred_tst,pred_qft,normtst,normqft)
d3<-gather(d3, test, norm_value, normtst:normqft)
d3<-mutate(d3,test=ifelse(test=="normtst","TST","IGRA"))
d3<-mutate(d3,pred=ifelse(test=="TST",pred_tst,pred_qft))
d3<-select(d3,norm_value,test,pred)
library(ggplot2)
ggplot(d3) +
geom_smooth(aes(x=norm_value, y=pred,color=test)) 

#Normalization following the PRESKOPE-TB percentile
#prepare breaks and labels for cut()
qft_cat<-c(0,0.001,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.13,0.15,0.17,0.2,0.22,0.26,0.3,0.35,0.4,0.45,0.52,0.6,0.69,0.81,0.92,1.07,1.27,1.52,1.77,2.06,2.48,2.94,3.55,4.28,5.1,6.29,7.59,9.18,10,+Inf)
qft_label<-c(1,39,40,49,54,57,60,62,63,65,66,67,69:98)
tst_cat<-c(0,2:23,25,27,31,+Inf)
tst_label<-c(1,52,53,55,57,60,62,64,66,68,72,74,77,79,81,84,86,88,89,91,94:98,99)

d2<-mutate(d2,normqft=cut(qfngit_tbag_nil,breaks=qft_cat,labels=qft_label,right=F))
d2<-mutate(d2,normtst=cut(mantoux_result,breaks=tst_cat,labels=tst_label,right=F))
d2$normtst<-as.numeric(as.character(d2$normtst))
d2$normqft<-as.numeric(as.character(d2$normqft))



fit_tst<-coxph(Surv(followup, activetb) ~ rcs(normtst,c(10,50,90)), data = d2)
d2$pred_tst[!is.na(d2$mantoux_result)] <- predict(fit_tst, data = d2, type="risk")
fit_qft<-coxph(Surv(followup, activetb) ~ rcs(normqft,c(10,50,90)), data = d2)
d2$pred_qft[!is.na(d2$normqft)] <- predict(fit_qft, data = d2, type="risk")
d3<-select(d2,pred_tst,pred_qft,normtst,normqft)
d3<-gather(d3, test, norm_value, normtst:normqft)
d3<-mutate(d3,test=ifelse(test=="normtst","TST","IGRA"))
d3<-mutate(d3,pred=ifelse(test=="TST",pred_tst,pred_qft))
d3<-select(d3,norm_value,test,pred)
library(ggplot2)
ggplot(d3) +
  geom_smooth(aes(x=norm_value, y=pred,color=test)) 
#sub-group-by contact


fit<-glm(activetb ~ mantoux_bin*contact+offset(log(followup)),family = poisson(link = "log"), data = d2)
tbl_regression(fit,exponentiate = T)

fit<-glm(activetb ~ igrapos*contact+offset(log(followup)),family = poisson(link = "log"), data = d2)
tbl_regression(fit,exponentiate = T)

#predictors
library(MASS)
d4<-select(d2,activetb,tst15,igrapos,contact,age,sex,prevtbdiag,prevbcg,bmi,followup)
gg_miss_var(d4)
res<-summary(aggr(d4, sortVar=TRUE))$combinations
d4<-d4[complete.cases(d4),]
full<- glm(activetb~tst15+igrapos+contact+age+sex+prevtbdiag+prevbcg+bmi+offset(log(followup)), d4,family = poisson(link="log"))
step.model <- stepAIC(full, direction = "both", 
                      trace = FALSE)
summary(step.model)

