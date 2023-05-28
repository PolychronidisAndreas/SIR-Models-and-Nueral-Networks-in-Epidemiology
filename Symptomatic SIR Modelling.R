#########################Symptomatic SIR Modelling#########################

#Load required libraries
library(ggplot2)
library(deSolve)
library(reshape2)

# Model inputs

initial_state_values=c(S=999999,I=1,R=0)
parameters=c(gamma=0.20,beta=0.50)

# Time points

time=seq(from=1,to=100,by=1)

# SIR model function 

sir_model <- function(time,state,parameters){
  with(as.list(c(state,parameters)),{
    N=S+I+R
    lambda=beta*(I/N) 
    dS=-lambda*S
    dI=lambda*S-gamma*I
    dR=gamma*I
    
    return(list(c(dS,dI,dR)))
  }
  )
}


#Solving the differential equations
output<-as.data.frame(ode(y=initial_state_values,func = sir_model,parms=parameters,times = time))


out_long=melt(output,id="time")
# To plot the proportion of susceptible, infected and recovered individuals over time
ggplot(data = out_long,          
       aes(x = time, y = value/1000000, colour = variable, group = variable)) +  
  geom_line() +xlab("Χρόνος (Μέρες)")+ylab("Ποσοστό του πληθυσμού")+scale_color_discrete(name="State")






