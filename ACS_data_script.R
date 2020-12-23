
d<-read.delim("../acs_clinical_DB.txt", header = T, sep = ",", dec = ".")

#manipulate to match variable names and coding with Rishi's data dictionary

d[d=="999"|d=="998"]<-NA
d[d==""]<-NA
library(tidyr)
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
          qfngit_result=QuantiferonInTubeResult,
          qfngit_tbag_nil=QFT_TBAg_Nil,
          igrapos=QuantiferonInTubeResult)
d[,c("studytime","dateofdeath","ltfu_date","censordate","recruitmentdate")]<-lapply(d[,c("studytime","dateofdeath","ltfu_date","censordate","recruitmentdate")],as.Date.character)
d<-mutate(d,diabetes=ifelse(ChronicIllnessDetails=="diabetes",1,0),
          bmi=weight/(height^2),
          died=ifelse(!is.na(dateofdeath),1,0),
          ltfu=ifelse(!is.na(ltfu_date),1,0),
          activetb=ifelse(TBCaseType==incident,1,0),
          activetbdate=ifelse(activetb==1,studytime,NA),
          mantoux_bin=iflese(mantoux_result>=5,1,0),
          
