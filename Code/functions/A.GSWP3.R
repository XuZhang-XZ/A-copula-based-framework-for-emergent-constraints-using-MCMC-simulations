
## A.GSWP3
A.GSWP3 = function(i.precip_1,i.time1) {
  
  ## Date
  i.date = data.frame("Date" = i.time1) %>%
    dplyr::mutate(year = year(Date), month = month(Date))
  i.month = i.date %>%
    dplyr::group_by(year, month) %>%
    dplyr::summarise()
  
  ## Return Data
  i.precip_re = array(NA,c(dim(i.precip_1)[c(1,2)],NROW(i.month)))
  
  ## Monthly
  i = 1
  for(i in 1:NROW(i.month)) {
    sel.ge = which(i.date$year == i.month$year[i] & i.date$month == i.month$month[i])
    i.precip_re[,,i] = apply(i.precip_1[,,sel.ge],c(1,2),mean)
  }
  
  ## Return
  i.precip_re
}
