# Environmental Modelling
# Project 3 

# Note: migration (= dispersal parameter) between sub-populations should be ignored

##### a) #####

# Define Parameters:

s1 <- 0.27
s2 <- 0.272
s3 <- 0.305
g1 <- 0.398
g2 <- 0.366
f3 <- 9.74
d <- 1 

# Lefkovitch matrix for cherry guava crop: 

cherry_MTR = rbind(c(s1*d, 0.000, f3*d),c(g1*d, s2*d, 0.000),c(0.000, g2*d, s3))
NAMES = c('Seedlings', 'Juveniles', 'Matures')
colnames(cherry_MTR) <-NAMES             
rownames(cherry_MTR) <- NAMES

# Starting population

ini_N <- matrix(c(1000,100,10), ncol=1)         
rownames(ini_N) <- NAMES            
colnames(ini_N) <- "Abundance"

# Matrix for the population development

YEARs <- 25
pop_MTR <- matrix(0,nrow=3,ncol=YEARs+1)
rownames(pop_MTR) <- NAMES
colnames(pop_MTR) <- seq(0,YEARs)
pop_MTR[,1] <- ini_N                   

# Calculation of the population development over time via for-loop

for(t in 2:(YEARs+1)){  
  
  pop_MTR[,t] <-  cherry_MTR %*% pop_MTR[,t-1]     
  
}

plot(colnames(pop_MTR),pop_MTR[1,], main = "No Removal", xlab="Time [Years]", ylab="Population size [#]",col="orange", type = "l", lwd=2, lty=1, ylim = c(0,max(pop_MTR)*1.05))
lines(colnames(pop_MTR),pop_MTR[2,],col="green", type = "l", lwd=2, lty=1)
lines(colnames(pop_MTR),pop_MTR[3,],col="brown", type = "l", lwd=2, lty=1)
legend("topleft", NAMES, col = c("orange","green","brown"), lty=1)


##### b) #####

# Finite rate of exponential growth

Lambda <- as.numeric(eigen(cherry_MTR)$values[1])

# Stable stage structure

SAD <- as.numeric(eigen(cherry_MTR)$vectors[,1])
SSS <- SAD/sum(SAD)


##### c) #####

# Introduce removal

Removal <- 0.35 # 0.1

for(t in 2:(YEARs+1)){  
  
  pop_MTR[,t] <-  cherry_MTR %*% pop_MTR[,t-1] %*% (1-Removal)  
  
}

plot(colnames(pop_MTR),pop_MTR[1,], main = "35% Removal"  , xlab="", ylab="",col="orange", type = "l", lwd=2, lty=1,ylim = c(0,max(pop_MTR)*1.05))
lines(colnames(pop_MTR),pop_MTR[2,],col="green", type = "l", lwd=2, lty=1)
lines(colnames(pop_MTR),pop_MTR[3,],col="brown", type = "l", lwd=2, lty=1)
legend("topright", NAMES, col = c("orange","green","brown"), lty=1)


##### d) #####

# Add parameters for density dependency (d)

gammaS <- 8.76*10^(-5)
gammaJ <- 1.02*10^(-4)
gammaM <- 3.65*10^(-4)
A <- 9 # [ha]

# run for-loop with density dependency build in

for(t in 2:(YEARs+1)){  
  
  d <- 1 - ((gammaS*pop_MTR[1,t-1]+gammaJ*pop_MTR[2,t-1]+gammaM*pop_MTR[3,t-1])/A)
  cherry_MTR = rbind(c(s1*d, 0.000, f3*d),c(g1*d, s2*d, 0.000),c(0.000, g2*d, s3))
  pop_MTR[,t] <-  cherry_MTR %*% pop_MTR[,t-1]
  
}

plot(colnames(pop_MTR),pop_MTR[1,],xlab="Time [Years]", ylab="Population size [#]", main = "Density dependent", col="orange", type = "l", lwd=2, lty=1,ylim = c(0,max(pop_MTR)*1.05))
lines(colnames(pop_MTR),pop_MTR[2,],col="green", type = "l", lwd=2, lty=1)
lines(colnames(pop_MTR),pop_MTR[3,],col="brown", type = "l", lwd=2, lty=1)
legend("topleft", NAMES, col=c("orange","green","brown"), lty=1)


##### extra #####

# finding the threshold (equilibrium) removal rate

# Reset the parameter d and the Lefkovitch matrix

d <- 1
cherry_MTR = rbind(c(s1*d, 0.000, f3*d),c(g1*d, s2*d, 0.000),c(0.000, g2*d, s3))

# Run the model with initial abundances consistent with SSS and Removal balancing the gowth rate (Lambda)

pop_MTR[,1] <- SSS*1000         
Removal <- 1-(1/Lambda)


for(t in 2:(YEARs+1)){  
  
  pop_MTR[,t] <-  cherry_MTR %*% pop_MTR[,t-1] %*% (1-Removal)  
  
}

plot(colnames(pop_MTR),pop_MTR[1,], main="Extra 1", xlab="Time [Years]", ylab="Population size [#]",col="orange", type = "l", lwd=2, lty=1,ylim = c(0,max(pop_MTR)*1.05))
lines(colnames(pop_MTR),pop_MTR[2,],col="green", type = "l", lwd=2, lty=1)
lines(colnames(pop_MTR),pop_MTR[3,],col="brown", type = "l", lwd=2, lty=1)
legend("topright", NAMES, col = c("orange","green","brown"), lty=1)


##### extra 2 #####

# effect of d on Lambda

# reset initial abundances

pop_MTR[,1] <- ini_N  

# create a matrix for the storage of the parameter values

parameters <- matrix(0,nrow=2,ncol=YEARs+1)
rownames(parameters) <- c("d","Lambda")
colnames(parameters) <- seq(0,YEARs)
parameters[,1] <- c(d,Lambda) 

# run for-loop for density dependency, save values of d and Lambda in matrix

for(t in 2:(YEARs+1)){  
  
  d <- 1 - ((gammaS*pop_MTR[1,t-1]+gammaJ*pop_MTR[2,t-1]+gammaM*pop_MTR[3,t-1])/A)
  cherry_MTR = rbind(c(s1*d, 0.000, f3*d),c(g1*d, s2*d, 0.000),c(0.000, g2*d, s3))
  pop_MTR[,t] <-  cherry_MTR %*% pop_MTR[,t-1]
  
  parameters[1,t] <- d
  
  Lambda <- as.numeric(eigen(cherry_MTR)$values[1])
  parameters[2,t] <- Lambda
  
}

plot(colnames(parameters),parameters[2,], main="Extra 2", xlab="Time [Years]", ylab="Parameter Values",col="blue", type = "l", lwd=2, lty=1,ylim = c(0.5,max(parameters)*1.05))
lines(colnames(parameters),parameters[1,],col="green", type = "l", lwd=2, lty=1)
legend("bottomleft", c("Lambda","d"), col = c("blue","green"), lty=1)
abline(1,0,lty=2)
