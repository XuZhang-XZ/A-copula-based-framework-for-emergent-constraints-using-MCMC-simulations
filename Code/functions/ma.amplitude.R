

# x = rnorm(100)
# j = 20
# ave.length = 20

ma.amplitude = function(x,ave.length = 20,var = "sd") {
  
  x0 = array(NA,length(x))
  
  for(j in ave.length:length(x)) {
    sel.ges = c((j - ave.length + 1) : j)
    if(var == "sd") {
      x0[j] = sd(x[sel.ges])
    }
    if(var == "var") {
      x0[j] = var(x[sel.ges])
    }
  }
  
  x0
  
  
}
