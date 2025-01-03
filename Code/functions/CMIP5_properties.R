
## Functions
CMIP5_properties = function(CMIP.name) {
  
  proper=strsplit(CMIP.name,"_")[[1]]
  
  stats.model1 = data.frame()
  stats.model1[1,"File.name"] = CMIP.name
  
  stats.model1[1,c("variable","relam","model")]= proper[1:3]
  stats.model1[1,c("scenario","varivant","resolution")]= proper[4:6]
  
  temp=strsplit(proper[6],split=c(".n"))[[1]]
  A.date =as.numeric(strsplit(temp,split=c("-"))[[1]])
  stats.model1[1,c("startdate","enddate")]=A.date
  stats.model1[1,c("startyear","endyear")]=c(year(ymd(A.date[1]*100+1)),year(ymd(A.date[2]*100+1)))
  stats.model1[1,c("startmonth","endmonth")]=c(month(ymd(A.date[1]*100+1)),month(ymd(A.date[2]*100+1)))
  
  int <- interval(ymd(A.date[1]*100+1),ymd(A.date[2]*100+1))
  stats.model1[1,"totalmonth"]=time_length(int, "month")+1
  
  return(stats.model1)
}
