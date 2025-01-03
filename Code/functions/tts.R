


## Group_by
tts = function(x) {
  n.len = length(x)
  if(n.len == 1) {
    re.x = x
  } else {
    re.x = x[1]
    for(i in 2:n.len){
      re.x = paste0(re.x,"_",x[i])
    }
  }
  return(re.x)
}
