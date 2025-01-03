


## Functions
EC.fun = function(His,Future,Obs.mean,Int.Var = NA,Type) {
  
  ## Internal Variability
  if(is.na(Int.Var)) {
    return("Errors")
  }
  
  ## Variance
  Obs_mean = Obs.mean
  Int.Var = Int.Var
  His_mean = mean(His)
  His_Var = var(His)
  Future_mean = mean(Future)
  Future_Var = var(Future)
  
  ## Properties
  Cor_EC = cor10.1(His,Future)[1]

  ## Slope and Relative reduction in Variance
  k_EC =  Cor_EC * std(His) * std(Future) / (var(His) + Int.Var)
  RRV = (Cor_EC)^2 / (1 +  Int.Var / var(His))
  
  ## Future EC_Mean and Var
  Future_mean_EC = mean(Future) + k_EC * (Obs.mean - mean(His))
  Future_Var_EC = var(Future) * (1-RRV)
  
  ## Return
  HEC_Summary = data.frame(
    "Type" = Type,
    "Obs.mean" = Obs_mean,
    "Int.Var" = Int.Var,
    "His_mean" = His_mean,
    "His_Var" = His_Var,
    "Future_mean" = Future_mean,
    "Future_Var" = Future_Var,
    "k_EC" = k_EC,
    "RRV" = RRV,
    "Future_mean_EC" = Future_mean_EC,
    "Future_Var_EC" = Future_Var_EC
  )
  
  ## Return
  return(HEC_Summary)

}
