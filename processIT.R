
source("configIt.R")
configs <- setConfig_f()

slist<-c("E10.5_Mnd_D","E10.5_Mnd_P","E10.5_Max_D","E10.5_Max_P","E11.5_Mnd_D","E11.5_Mnd_P","E11.5_Max_D","E11.5_Max_P","E12.5_Mnd_D","E12.5_Mnd_P","E12.5_Max_D","E12.5_Max_P","E13.5_Mnd_D","E13.5_Mnd_P","E13.5_Max_D","E13.5_Max_P","E14.5_Mnd_D","E14.5_Mnd_P","E14.5_Max_D","E14.5_Max_P")

c <- 1
cc <- 2
while ( cc < 20 ) { 
  sel <-c(slist[c:cc])
  c <- c+2
  cc <- cc+2

  configs$sel <-sel
  print(configs$sel)
  source("processCEL.R")
  source("process2JSON.R")
  processCELdata_f(configs)
}


