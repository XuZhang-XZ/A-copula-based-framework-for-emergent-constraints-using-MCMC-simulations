

## Var Functions
Running.var = function(x,n = 30){
  x.var = array(NA,length(x))
  for(i in (n+1):(length(x))){
    x.var[i] = sd(x[(i-n):(i)])
  }
  return(x.var)
}



## Var Functions
Running.trend = function(x,n = 30){
  x.var = array(NA,length(x))
  for(i in (n+1):(length(x))){
    x.var[i] = linear.trend(x[(i-n):(i)])
  }
  return(x.var)
}
