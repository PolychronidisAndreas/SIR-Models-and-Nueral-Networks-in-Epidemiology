#########################Asymptomatic SIR Modelling############################

#Load required libraries
library(ggplot2)
library(deSolve)
library(reshape2)

# Model inputs

initial_state_values=c(S=999999,I=1,J=0,R=0,U=0)
parameters=c(gamma=0.8,beta=0.4,v=0,ksi=0.2,mu=0.8)

# Time points

time=seq(from=1,to=15,by=1)

# SIR model function 

sir_model <- function(time,state,parameters){
  with(as.list(c(state,parameters)),{
    N=S+I+J+R+U
    dS=-beta*(I+J)*S-v*S
    dI=beta*ksi*(I+J)*S-gamma*I
    dJ=beta*(1-ksi)*(I+J)*S-mu*J
    dR=gamma*I+v*ksi*S
    dU=mu*J+v*(1-ksi)*S
    return(list(c(dS,dI,dJ,dR,dU)))
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








