######################################################
#Data Generation Functions
########################################################


dataGenerationPairedT <- function(normality, n, hypothesis){
  Ha <- 4
  sd <- 10
  if(normality){
    data <- rnorm(n = n, mean = ifelse(hypothesis == 'H0', 0, Ha), sd = sd)
  }else{
    sampleMeanAndSD <- rnorm(n = n, mean = ifelse(hypothesis == 'H0', 0, Ha), sd = sd)
    mean <- mean(sampleMeanAndSD)
    SD <- sd(sampleMeanAndSD)
    data <- sn::rst(n = n, xi = 0, omega = 1, alpha = -30)
    data <- mean + SD*scale(data)
  }
  return(data)
}


dataGenerationRANOVA <- function(normality, k, sphericity,  n, hypothesis){
  muHa <- seq(0,6.5, 6.5/(k-1))
  muh0 <- rep(0, k)
  sd <- 100
  
  if(hypothesis == 'H0'){
    mu <- muh0
  }
  else{
    mu <- muHa
  }
  
  if(normality && sphericity){
    if(k == 3){
      sigma <- matrix(c(sd,2,2,
                        2,sd,2,
                        2,2,sd),3,3)
      
    }else{
      sigma <- matrix(c(sd,2,2,2,2,
                        2,sd,2,2,2,
                        2,2,sd,2,2,
                        2,2,2,sd,2,
                        2,2,2,2,sd),5,5)
      
    }
    data <- MASS::mvrnorm(n = n, mu = mu, Sigma = sigma)
  }else if(normality && !sphericity){
    if(k == 3){
      sigma <- matrix(c(sd,30,40,
                        30,sd,60,
                        40,60,sd),3,3)
      
    }else{
      sigma <- matrix(c(sd,20,40,60,80,
                        20,sd,40,40,60,
                        40,40,sd,60,40,
                        60,40,60,sd,80,
                        80,60,40,80,sd),5,5)
      
    }
    data <- MASS::mvrnorm(n = n, mu = mu, Sigma = sigma)
  }else if(sphericity && !normality){
    if(k == 3){
      sigma <- matrix(c(sd,2,2,
                        2,sd,2,
                        2,2,sd),3,3)
      
    }else{
      sigma <- matrix(c(sd,2,2,2,2,
                        2,sd,2,2,2,
                        2,2,sd,2,2,
                        2,2,2,sd,2,
                        2,2,2,2,sd),5,5)
      
    }
    D <- diag(rep(-350,k)) # skew parameter
    
    sampleMeansAndSDs <- MASS::mvrnorm(n = n, mu = mu, Sigma = sigma)
    means <- apply(sampleMeansAndSDs, 2, mean)
    SDs <- apply(sampleMeansAndSDs, 2, sd)
    
    data <- miscF::rmvsn(n = n, D = D, Mu = mu, Sigma = sigma)
    
    for(i in 1:ncol(data)){
      data[,i] <- means[i] + SDs[i]*scale(data[,i])
    }
    
  }
  
  dataframe <- data.frame(ID = 1:nrow(data))
  for(i in 1:ncol(data)){
    dataframe <- cbind(dataframe, data[,i])
  }
  
  colnames(dataframe) <- c('ID', paste('v', 1:k, sep = ''))
  
  return(dataframe)
}

################################################################################
#Paired t-test Simulation
################################################################################


repetitions <- 1:100
n <- c(15, 40, 150)
normality <- c(TRUE, FALSE)
hypothesis <- c('H0', 'Ha')


options <- expand.grid(n = n, normality = normality, hypothesis = hypothesis, repetitions = repetitions, stringsAsFactors = F)
results <- data.frame(n = NA, hypothesis = NA, normality = NA, p = NA, mean = NA, BF = NA, Shapiro = NA)

for(i in 1:nrow(options)){
  set.seed(i)
  data <- dataGenerationPairedT(normality = options[["normality"]][i], n = options[["n"]][i], 
                                hypothesis = options[["hypothesis"]][i])
  p.value <- t.test(data, mu = 0)$p.value
  Bttest <- BayesFactor::ttestBF(x = data, mu = 0)
  BF <- BayesFactor::extractBF(Bttest)$bf
  results[i,] <- list(n = options[["n"]][i], hypothesis = options[["hypothesis"]][i], normality = options[["normality"]][i],
                      p = p.value, mean = mean(data), BF = BF, Shapiro = shapiro.test(data)$p.value)
}


########################################################################
#RESULTS PAIRED T TEST
########################################################################
results <- read.csv("pairedTresults.csv")

par(mfrow=c(2,3))


#Normal Data


Normal <- subset(results, results$normality == TRUE)
length(Normal$n[Normal$Shapiro < 0.05])/nrow(Normal) #positive shapiro tests


#Skewed Data
Skewed <- subset(results, results$normality == FALSE)
length(Skewed$n[Skewed$Shapiro < 0.05])/nrow(Skewed) #positive shapiro tests

#Skewed H0
H0s <- subset(Skewed, Skewed$hypothesis == 'H0')

#Normal H0
H0n <- subset(Normal, Normal$hypothesis == 'H0')

#Normal Ha
Han <- subset(Normal, Normal$hypothesis == 'Ha')


#Skewed Ha
Has <- subset(Skewed, Skewed$hypothesis == 'Ha')


#H0 Normal vs Skewed Frequentist
#plotting
#pdf("pairedTtestH0.png")

m <- matrix(c(1,1,1,2,3,4,5,6,7),nrow = 3,ncol = 3,byrow = TRUE)

layout(mat = m)

plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x = "bottom",inset = 0,
       legend = c("Normal", "Skewed"),lty = c(1,2),lwd=2.5, cex=1.5, horiz = TRUE)

type1errors <- data.frame()
for(i in 1:length(n)){
  nH0n <- subset(H0n, H0n$n == n[i])
  type1 <- length(nH0n$n[nH0n$p < 0.05])/nrow(nH0n) # type 1 errors
  plot(density((nH0n$p), from = 0 ), lty = 1,lwd = 1, main = paste('n =',n[i]), xlab = "", xlim = c(0,1))
  abline(v = 0.05)
  type1errors <- rbind(type1errors, list(n = n[i], shape = "Normal", type1 = type1))

  nH0s <- subset(H0s, H0s$n == n[i])
  type1 <- length(nH0s$n[nH0s$p < 0.05])/nrow(nH0s) # type 1 errors
  lines(density((nH0s$p), from = 0), lty = 2 ,lwd = 1,)
  abline(v = 0.05)
  type1errors <- rbind(type1errors, list(n = n[i], shape = "Skewed", type1 = type1))
  
}
type1errors

#H0 Normal vs Skewed Bayesian
type1errors <- data.frame()
for(i in 1:length(n)){
  nH0n <- subset(H0n, H0n$n == n[i])
  type1 <- length(nH0n$n[nH0n$BF > 3])/nrow(H0n) #type 1 errors
  plot(density(log(nH0n$BF)), lty = 1,lwd = 1, main = paste('n =', n[i]), xlab = "")
  abline(v = -log(3))
  abline(v = log(3))
  type1errors <- rbind(type1errors, list(n = n[i], shape = "Normal", type1 = type1))

  nH0s <- subset(H0s, H0s$n == n[i])
  type1 <- length(nH0s$n[nH0s$BF > 3])/nrow(nH0s) #type I errors
  lines(density(log(nH0s$BF)), lty = 2,lwd = 1)
  abline(v = -log(3))
  abline(v = log(3))
  type1errors <- rbind(type1errors, list(n = n[i], shape = "Skewed", type1 = type1))
  
}
type1errors

#dev.off()

#Ha Normal vs Skewed Frequentist
#plotting
#pdf("pairedTtestHa.pdf")
m <- matrix(c(1,1,1,2,3,4,5,6,7),nrow = 3,ncol = 3,byrow = TRUE)

layout(mat = m)

plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x = "bottom",inset = 0,
       legend = c("Normal", "Skewed"),lty = c(1,2),lwd=2.5, cex=1.5, horiz = TRUE)

powerFrame <- data.frame()
for(i in 1:length(n)){
  nHan <- subset(Han, Han$n == n[i])
  power <- 1 - length(nHan$n[nHan$p > 0.05])/nrow(nHan) # power: 1 - type 2 error
  plot(density(nHan$p, from = 0), lty = 1,lwd = 1, main = paste('n =',n[i]), xlab = "", xlim = c(0,0.5))
  abline(v = 0.05)
  powerFrame <- rbind(powerFrame, list(n = n[i], shape = "Normal", power = power))

  nHas <- subset(Has, Has$n == n[i])
  power <- 1 - length(nHas$n[nHas$p > 0.05])/nrow(nHas) # power: 1 - type 2 error
  lines(density(nHas$p, from = 0), lty = 2,lwd = 1,)
  abline(v = 0.05)
  powerFrame <- rbind(powerFrame, list(n = n[i], shape = "Skewed", power = power))
}
powerFrame

#Ha Normal vs Skewed Bayesian
powerFrame <- data.frame()

for(i in 1:length(n)){
  nHan <- subset(Han, Han$n == n[i])
  power <- 1 - length(nHan$n[nHan$BF < 3])/nrow(nHan) # power: 1 - type 2 error
  plot(density(log(nHan$BF)), lty = 1,lwd = 1, main = paste('n =',n[i]), xlab = "", xlim = c(-3, 15))
  abline(v = log(3))
  abline(v = -log(3))
  powerFrame <- rbind(powerFrame, list(n = n[i], shape = "Normal", power = power))

  nHas <- subset(Has, Has$n == n[i])
  power <- 1- length(nHas$n[nHas$BF < 3])/nrow(nHas) #true positive
  lines(density(log(nHas$BF)), lty = 2,lwd = 1, main = paste('n =',n[i]))
  abline(v = log(3))
  abline(v = -log(3))
  powerFrame <- rbind(powerFrame, list(n = n[i], shape = "Skewed", power = power))
}
powerFrame
#dev.off()

#######################################################################
#Repeated Measures ANOVA
#######################################################################



repetitions <- 1:100
n <- c(15, 45, 150)
k <- c(3,5)
sphericity <- c(TRUE,FALSE)
normality <- c(TRUE,FALSE)
hypothesis <- c('H0', 'Ha')

options <- expand.grid(n = n, normality = normality, sphericity = sphericity, hypothesis = hypothesis, repetitions = repetitions, k = k, stringsAsFactors = F)
options <- subset(options, !(options$normality == F & options$sphericity == F))
results <- data.frame(n = NA, k = NA, normality = NA, sphericity = NA, 
                      hypothesis = NA, p = NA, BF = NA, mauchly = NA, mardia = NA, eta = NA)

for(i in 1:nrow(options)){
  set.seed(i)
  print(paste(round((i/nrow(options)*100),digits = 2),'%')) #percent complete
  
  data <- dataGenerationRANOVA(normality = options[["normality"]][i], k = options[["k"]][i],
                               sphericity = options[["sphericity"]][i], n = options[["n"]][i], hypothesis = options[["hypothesis"]][i])
  
  mardia <- psych::mardia(data, plot = F)$p.skew
  reshapedData <- tidyr::gather(data = data, key = "key", value = "value", 2:ncol(data))
  reshapedData$ID <- as.factor(reshapedData$ID)
  reshapedData$key <- as.factor(reshapedData$key)
  freqANOVA <- rstatix::anova_test(data = reshapedData, formula = value ~ key + Error(ID/key))
  pvalue <- freqANOVA$ANOVA$p
  mauchly <- freqANOVA$`Mauchly's Test for Sphericity`$p
  eta <- freqANOVA$ANOVA$ges
  bANOVA <- BayesFactor::anovaBF(value ~ key + ID, data = reshapedData, whichRandom = "ID", progress = F)
  BF <- BayesFactor::extractBF(bANOVA)$bf
  
  results[i,] <- list(n = options[["n"]][i], k = options[["k"]][i], normality = options[["normality"]][i],
                      sphericity = options[["sphericity"]][i], hypothesis = options[["hypothesis"]][i], 
                      p = pvalue, BF = BF, mauchly = mauchly, mardia = mardia, eta = eta)
}




###################################################
#RESULTS RANOVA
###################################################





results <- read.csv("RANOVAsimulation.csv")

par(mfrow=c(3,2))

#Manipulation check
#NothingViolated
nonV <- subset(results, results$sphericity == TRUE & results$normality == TRUE)
length(nonV$n[nonV$mauchly < 0.05])/nrow(nonV) # positive mauchly test
length(nonV$n[nonV$mardia < 0.05])/nrow(nonV)  # positive mardia test

#SphericityViolated
spherV <- subset(results, results$sphericity == FALSE)
length(spherV$n[spherV$mauchly < 0.05])/nrow(spherV) #positive mauchly test
length(spherV$n[spherV$mardia < 0.05])/nrow(spherV)  # positive mardia test

#NormalityViolated
normV <- subset(results, results$normality == FALSE)
length(normV$n[normV$mauchly < 0.05])/nrow(normV) #positive mauchly test
length(normV$n[normV$mardia < 0.05])/nrow(normV)  # positive mardia test

#n150
n150 <- subset(results, results$n == 150)
normV150 <- subset(n150, n150$normality == FALSE)
length(normV150$n[normV150$mardia < 0.05])/nrow(normV150)  # positive mardia test




#Under H0
H0data <- subset(results, results$hypothesis == 'H0')


##Nothing Violated
nonViolatedH0 <- subset(H0data, H0data$sphericity == TRUE & H0data$normality == TRUE)

##Sphericity Violated
spherViolatedH0 <- subset(H0data, H0data$sphericity == FALSE)

##Normality Violated
normViolatedH0 <- subset(H0data, H0data$normality == FALSE)


#H0 Nothing violated vs Sphericity Violated vs Normality Violated Frequentist

#pdf('ANOVAh0Freq.pdf')


m <- matrix(c(1,1,1,2,4,6,3,5,7),nrow = 3,ncol = 3,byrow = TRUE)

layout(mat = m)

plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x = "bottom",inset = 0,
       legend = c("Ideal", "Sphericity Violated", "Normality Violated" ),lty = c(1,2,3),lwd=2.5, cex=1.5, horiz = TRUE)

type1errors <- data.frame()

for(i in 1:length(n)){
  sampleSize <- n[i]
  for(j in 1:length(k)){
    subgroupSize <- k[j]
    nnonViolatedH0 <- subset(nonViolatedH0, nonViolatedH0$n == sampleSize & nonViolatedH0$k == subgroupSize)
    type1 <- length(nnonViolatedH0$n[nnonViolatedH0$p < 0.05])/nrow(nnonViolatedH0) #type I errors
    plot(density((nnonViolatedH0$p), from = 0, to = 1), lty = 1,lwd = 1, main = paste("n = ", n[i], 'k =', k[j]), xlab="")
    abline(v = 0.05)
    type1errors <- rbind(type1errors, list(n = n[i], k = k[j], violations = "None", type1 = type1))
    
    nspherViolatedH0 <- subset(spherViolatedH0, spherViolatedH0$n == sampleSize & spherViolatedH0$k == subgroupSize)
    type1 <- length(nspherViolatedH0$n[nspherViolatedH0$p < 0.05])/nrow(nspherViolatedH0) #type I errors
    lines(density((nspherViolatedH0$p)), lty = 2, lwd = 1)
    type1errors <- rbind(type1errors, list(n = n[i], k = k[j], violations = "Sphericity", type1 = type1))
    
    nnormViolatedH0 <- subset(normViolatedH0, normViolatedH0$n == sampleSize & normViolatedH0$k == subgroupSize)
    type1 <- length(nnormViolatedH0$n[nnormViolatedH0$p < 0.05])/nrow(nnormViolatedH0) #type I errors
    lines(density((nnormViolatedH0$p)), lty = 3, lwd = 1)
    type1errors <- rbind(type1errors, list(n = n[i], k = k[j], violations = "Normality", type1 = type1))
  }
}
type1errors

#dev.off()

#H0 Nothing violated vs Sphericity Violated vs Normality Violated Bayesian

#pdf('ANOVAh0Bayes.pdf')

m <- matrix(c(1,1,1,2,4,6,3,5,7),nrow = 3,ncol = 3,byrow = TRUE)

layout(mat = m)


plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x = "bottom",inset = 0,
       legend = c("Ideal", "Sphericity Violated", "Normality Violated" ),lty = c(1,2,3),lwd=2.5, cex=1.5, horiz = TRUE)


type1errors <- data.frame()

for(i in 1:length(n)){
  sampleSize <- n[i]
  for(j in 1:length(k)){
    subgroupSize <- k[j]
    nnonViolatedH0 <- subset(nonViolatedH0, nonViolatedH0$n == sampleSize & nonViolatedH0$k == subgroupSize)
    type1 <- length(nnonViolatedH0$n[nnonViolatedH0$BF > 3])/nrow(nnonViolatedH0) #type I errors
    plot(density(log(nnonViolatedH0$BF)), lty = 1,lwd = 1, main = paste("n = ", n[i], 'k =', k[j]), xlab = '')
    abline(v = -log(3))
    abline(v = log(3))
    type1errors <- rbind(type1errors, list(n = n[i], k = k[j], violations = "None", type1 = type1))
    
    nspherViolatedH0 <- subset(spherViolatedH0, spherViolatedH0$n == sampleSize & spherViolatedH0$k == subgroupSize)
    type1 <- length(nspherViolatedH0$n[nspherViolatedH0$BF > 3])/nrow(nspherViolatedH0) #type I errors
    lines(density(log(nspherViolatedH0$BF)), lty = 2,lwd = 1, main = paste("n = ", n[i], 'k =', k[j]))
    abline(v = -log(3))
    abline(v = log(3))
    type1errors <- rbind(type1errors, list(n = n[i], k = k[j], violations = "Sphericity", type1 = type1))
    
    nnormViolatedH0 <- subset(normViolatedH0, normViolatedH0$n == sampleSize & normViolatedH0$k == subgroupSize)
    type1 <- length(nnormViolatedH0$n[nnormViolatedH0$BF > 3])/nrow(nnormViolatedH0) #type I errors
    lines(density(log(nnormViolatedH0$BF)), lty = 3,lwd = 1, main = paste("n = ", n[i], 'k =', k[j]))
    abline(v = -log(3))
    abline(v = log(3))
    type1errors <- rbind(type1errors, list(n = n[i], k = k[j], violations = "normality", type1 = type1))
  }
}
type1errors

#dev.off()

#Under Ha
Hadata <- subset(results, results$hypothesis == 'Ha')
mean(Hadata$eta)


##Nothing Violated
nonViolatedHa <- subset(Hadata, Hadata$sphericity == TRUE & Hadata$normality == TRUE)

##Sphericity Violated
spherViolatedHa <- subset(Hadata, Hadata$sphericity == FALSE)

##Normality Violated
normViolatedHa <- subset(Hadata, Hadata$normality == FALSE)


#Ha Nothing violated vs Sphericity Violated vs Normality Violated Frequentist

#pdf('ANOVAhaFreq.pdf')

m <- matrix(c(1,1,1,2,4,6,3,5,7),nrow = 3,ncol = 3,byrow = TRUE)

layout(mat = m)


plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x = "bottom",inset = 0,
       legend = c("Ideal", "Sphericity Violated", "Normality Violated" ),lty = c(1,2,3),lwd=2.5, cex=1.5, horiz = TRUE)


powerFrame <- data.frame()

for(i in 1:length(n)){
  sampleSize <- n[i]
  for(j in 1:length(k)){
    subgroupSize <- k[j]
    nnonViolatedHa <- subset(nonViolatedHa, nonViolatedHa$n == sampleSize & nonViolatedHa$k == subgroupSize)
    power <- 1 - length(nnonViolatedHa$n[nnonViolatedHa$p > 0.05])/nrow(nnonViolatedHa) #power
    plot(density((nnonViolatedHa$p), from = 0, to = 1), lty = 1,lwd = 1, main = paste("n = ", n[i], 'k =', k[j]), xlab = '')
    abline(v = 0.05)
    powerFrame <- rbind(powerFrame, list(n = n[i], k = k[j], violations = "None", power = power))
    
    nspherViolatedHa <- subset(spherViolatedHa, spherViolatedHa$n == sampleSize & spherViolatedHa$k == subgroupSize)
    power <- 1 -length(nspherViolatedHa$n[nspherViolatedHa$p > 0.05])/nrow(nspherViolatedHa) #power
    lines(density((nspherViolatedHa$p)), lty = 2, lwd = 1)
    powerFrame <- rbind(powerFrame, list(n = n[i], k = k[j], violations = "Sphericity", power = power))
    
    nnormViolatedHa <- subset(normViolatedHa, normViolatedHa$n == sampleSize & normViolatedHa$k == subgroupSize)
    power <- 1 -length(nnormViolatedHa$n[nnormViolatedHa$p > 0.05])/nrow(nnormViolatedHa) #power
    lines(density((nnormViolatedHa$p)), lty = 3, lwd = 1)
    powerFrame <- rbind(powerFrame, list(n = n[i], k = k[j], violations = "Normality", power = power))
  }
}
powerFrame

#dev.off()

#Ha Nothing violated vs Sphericity Violated vs Normality Violated Bayesian

#pdf('ANOVAhaBayes.pdf')

m <- matrix(c(1,1,1,2,4,6,3,5,7),nrow = 3,ncol = 3,byrow = TRUE)

layout(mat = m)


plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x = "bottom",inset = 0,
       legend = c("Ideal", "Sphericity Violated", "Normality Violated" ),lty = c(1,2,3),lwd=2.5, cex=1.5, horiz = TRUE)


powerFrame <- data.frame()

for(i in 1:length(n)){
  sampleSize <- n[i]
  for(j in 1:length(k)){
    subgroupSize <- k[j]
    nnonViolatedHa <- subset(nonViolatedHa, nonViolatedHa$n == sampleSize & nonViolatedHa$k == subgroupSize)
    power <- 1 - length(nnonViolatedHa$n[nnonViolatedHa$BF < 3])/nrow(nnonViolatedHa) #power
    plot(density(log(nnonViolatedHa$BF)), lty = 1,lwd = 1, main = paste("n = ", n[i], 'k =', k[j]), xlab = '')
    abline(v = -log(3))
    abline(v = log(3))
    powerFrame <- rbind(powerFrame, list(n = n[i], k = k[j], violations = "None", power = power))
    
    nspherViolatedHa <- subset(spherViolatedHa, spherViolatedHa$n == sampleSize & spherViolatedHa$k == subgroupSize)
    power <- 1 - length(nspherViolatedHa$n[nspherViolatedHa$BF < 3])/nrow(nspherViolatedHa) #power
    lines(density(log(nspherViolatedHa$BF)), lty = 2,lwd = 1, main = paste("n = ", n[i], 'k =', k[j]))
    abline(v = -log(3))
    abline(v = log(3))
    powerFrame <- rbind(powerFrame, list(n = n[i], k = k[j], violations = "Sphericity", power = power))
    
    nnormViolatedHa <- subset(normViolatedHa, normViolatedHa$n == sampleSize & normViolatedHa$k == subgroupSize)
    power <- 1 - length(nnormViolatedHa$n[nnormViolatedHa$BF < 3])/nrow(nnormViolatedHa) #power
    lines(density(log(nnormViolatedHa$BF)), lty = 3,lwd = 1, main = paste("n = ", n[i], 'k =', k[j]))
    abline(v = -log(3))
    abline(v = log(3))
    powerFrame <- rbind(powerFrame, list(n = n[i], k = k[j], violations = "Normality", power = power))
  }
}
powerFrame

#dev.off()
