

## Marginal Functions
Marignal.distribution = function(x = x,Type = "d", Num = "x_1") {
  
  ## Data
  Univariate.best_S = Univariate.best %>% dplyr::filter(Variable == Num)
  paras = c(Univariate.best_S$para1,Univariate.best_S$para2)
  
  if(Num == "x_1") {
    if(Type == "d") {
      re = dgamma(x,paras[1],paras[2])
    }
    if(Type == "p") {
      re = pgamma(x,paras[1],paras[2])
    }
  }
  if(Num == "x_2") {
    if(Type == "d") {
      re = dgamma(x,paras[1],paras[2])
    }
    if(Type == "p") {
      re = pgamma(x,paras[1],paras[2])
    }
  }
  if(Num == "y_1") {
    if(Type == "d") {
      re = dlnorm(x,paras[1],paras[2])
    }
    if(Type == "p") {
      re = plnorm(x,paras[1],paras[2])
    }
  }
  
  ## Return
  return(re)
  
}





# 
# 
# ## Functions
# Marignal.distribution = function(x = x,Type = "d", Num = "x_1") {
#   
#   ## Data
#   Univariate.best = read.csv("Output Data/CMIP6/tas/Univariate.best.csv",row.names = 1) %>%
#     dplyr::filter(Variable == Num)
#   paras = c(Univariate.best$para1,Univariate.best$para2)
#   
#   if(Num == "x_1") {
#     if(Type == "d") {
#       re = dgamma(x,paras[1],paras[2])
#     }
#     if(Type == "p") {
#       re = pgamma(x,paras[1],paras[2])
#     }
#   }
#   if(Num == "x_2") {
#     if(Type == "d") {
#       re = dgamma(x,paras[1],paras[2])
#     }
#     if(Type == "p") {
#       re = pgamma(x,paras[1],paras[2])
#     }
#   }
#   if(Num == "y_1") {
#     if(Type == "d") {
#       re = dlnorm(x,paras[1],paras[2])
#     }
#     if(Type == "p") {
#       re = plnorm(x,paras[1],paras[2])
#     }
#   }
# 
#   ## Return
#   return(re)
#   
# }


# 
# Marignal.distribution = function(x = x,Type = "d", Num = 1,Univariate.paras = Univariate.paras) {
#   
#   paras = Univariate.paras[1:2,Num]
#   
#   if(Type == "d") {
#     re = dcauchy(x,paras[1],paras[2])
#   }
#   if(Type == "p") {
#     re = pcauchy(x,paras[1],paras[2])
#   }
#   
#   # if(Num ==1 & Type == "d") {
#   #   re = dlnorm(x,-2.81549131,0.39166272)
#   # }
#   # 
#   # if(Num ==1 & Type == "p") {
#   #   paras = Univariate.paras[1:2,i]
#   #   re = plnorm(x,-2.81549131,0.39166272)
#   # }
#   # 
#   # if(Num ==2 & Type == "d") {
#   #   re = dlnorm(x,0.04573080,0.20125698)
#   # }
#   # 
#   # if(Num ==2 & Type == "p") {
#   #   re = plnorm(x,0.04573080,0.20125698)
#   # }
#   # 
#   # 
#   # if(Num ==3 & Type == "d") {
#   #   re = dlnorm(x,-0.25967926,0.55409044)
#   # }
#   # 
#   # if(Num ==3 & Type == "p") {
#   #   re = plnorm(x,-0.25967926,0.55409044)
#   # }
#   
#   ## Return
#   return(re)
# 
# }


# ## Test
# x = seq(0,3,0.1)
# plot(x,dlnorm(x,0.04573080,0.20125698))
