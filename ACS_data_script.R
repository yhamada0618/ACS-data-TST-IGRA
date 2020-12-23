
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
          sex=ifelse(sex=="male",0,1),
          qfngit_result=igrapos)
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
#exclude participants with prevalent TB
d2<-filter(d2,TBCaseType=="incident"|is.na(TBCaseType))
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
library("survminer")
fit<-coxph(Surv(followup, activetb) ~ mantoux_bin, data = d2)
tbl_regression(fit,exponentiate = T)

fit<-coxph(Surv(followup, activetb) ~ igrapos, data = d2)
tbl_regression(fit,exponentiate = T)

explore