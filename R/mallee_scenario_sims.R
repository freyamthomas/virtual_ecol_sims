
# Stochastic parameters

# Our modelled height values

Hmax <- c(145.42035, 103.65811, 115.94892, 89.81047, 459.25922, 64.70510, 48.75292, 81.85831, 258.23389, 267.63473, 81.73318, 42.15995, 211.01529, 62.83146, 86.18343, 84.96240, 57.87685, 41.01645, 143.17702, 59.73425)
a <- c(4.226643, 4.543715, 3.849143, 3.226788, 4.580450, 4.489764, 4.770788, 3.483222, 6.622986, 6.644306, 4.597916, 4.892076, 4.385591, 4.410173, 4.462377, 4.700410, 4.362824, 5.073423, 4.689330, 3.371535)
b <- c(1.720544, 1.555347, 1.818385, 2.122570, 1.818912, 1.441581, 1.265560, 1.927080, 1.133050, 1.136036, 1.551051, 1.203531, 1.670703, 1.481089, 1.515627, 1.431935, 1.479525, 1.149516, 1.570982, 1.903793)
Time <- c(0,1,2,3,4,6,8,13,15,26,28,33,36,41,86)
mu <- matrix(NA, nrow=length(Hmax), ncol=length(Time))
row.names(mu) <- c('Acabra', 'Acamon', 'Acawil', 'Beyopa', 'Codon', 'Dodbur', 'Erecrass', 'Eregla', 'EucBlue', 'EucGreen', 'Grehue', 'Halcya', 'Mellan', 'Olemul', 'Olepim', 'Olesub', 'Phebalium', 'ProstantheraGreen', 'Senart', 'Wesrig')
TSFvalues <- c(0,1,2,3,4,6,8,13,15,26,28,33,36,41,86)
colnames(mu) <- TSFvalues
for (j in 1:length(Hmax)) { 
  for (x in 1:length(Time))   {
    mu[j, x] <- Hmax[j] / (1 + exp(-a[j] * (Time[x] - b[j])))  
  }  
}


# Our detection time matrix

N.species <- dim(mu)[1]
TSFvalues <- c(0,1,2,3,4,6,8,13,15,26,28,33,36,41,86)
p_matrix <- matrix(NA, nrow=N.species, ncol=length(TSFvalues))
row.names(p_matrix) <- c('Acabra', 'Acamon', 'Acawil', 'Beyopa', 'Codon', 'Dodbur', 'Erecrass', 'Eregla', 'EucBlue', 'EucGreen', 'Grehue', 'Halcya', 'Mellan', 'Olemul', 'Olepim', 'Olesub', 'Phebalium', 'ProstantheraGreen', 'Senart', 'Wesrig')
TSFvalues <- c(0,1,2,3,4,6,8,13,15,26,28,33,36,41,86)
colnames(p_matrix) <- TSFvalues
for (r in 1:length(TSFvalues)) {
  
  p_matrix[,r] <- c(2, 10, 20, 5, 5, 2, 2, 2, 2, 2, 10, 3, 4, 5, 6, 7, 3, 5, 3, 5) 
  
}


# Our occupancy matrix

TSFvalues <- c(0,1,2,3,4,6,8,13,15,26,28,33,36,41,86)
Occupancy_meta_matrix <- matrix(NA, nrow=N.species, ncol=length(TSFvalues))
row.names(Occupancy_meta_matrix) <- c('Acabra', 'Acamon', 'Acawil', 'Beyopa', 'Codon', 'Dodbur', 'Erecrass', 'Eregla', 'EucBlue', 'EucGreen', 'Grehue', 'Halcya', 'Mellan', 'Olemul', 'Olepim', 'Olesub', 'Phebalium', 'ProstantheraGreen', 'Senart', 'Wesrig')
TSFvalues <- c(0,1,2,3,4,6,8,13,15,26,28,33,36,41,86)
colnames(Occupancy_meta_matrix) <- TSFvalues
for (w in 1:length(TSFvalues)) {
  
  Occupancy_meta_matrix[,w] <- c(0.9, 0.9, 0.9, 0.9, 0.9,0.9,0.9, 0.4, 0.4,0.4,0.4, 0.4, 0.4, 0.4, 0.4, 0.2, 0.2, 0.1, 0.1, 0.1) 
  
}

# Strategic sampling with a PriorityList of sites to visit

PriorityList <- matrix(NA, 1400, 2)
PriorityList[,1]  <- rep(c(2, 8, 26, 41, 4, 15, 33, 86, 1,6, 13, 3, 36, 28), 100)
PriorityList[,2] <- rep(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1), 100)
colnames(PriorityList) <- c('TSF', 'Rep')



######################################################################################

# Deterministic parameters

HrsPerDay <- 10 # long days
FieldMinsDay <- HrsPerDay * 60   
mu
sd.obs <- 0.395
PriorityList 
Occupancy_meta_matrix
p_matrix
N.fieldtrips <- 100  # replicate fieldtrips 100 times
TSFvalues <- c(0,1,2,3,4,6,8,13,15,26,28,33,36,41,86) 
HometoSiteTravel <- (30 * 2) # 30 minutes to get to site and home
setuptime <- 15 # 15 minutes to set up
mt <- 3.5 # 3.5 minutes to measure 
TravelbwSites <- 30 # 30 minutes to travel between sites
Maxtime <- 60 * 30   # 30 hours at one replicate
N.reps <- 1 # one replicate
N.TSF <- length(TSFvalues)


#####################################################################################

# Field day values set up the probability of achieving x number of real field days given 'unworkable' events

FieldDays7 <- 7
FieldDayvalues7 <- matrix(0, nrow=2, ncol=FieldDays7*2+1)
FieldDayvalues7[2,] <- c(rep(0,8), rep(1,3), rep(2,1),rep(3,1), rep(5,1), rep(10,1))
for (q in 2:length(FieldDayvalues7[1,])) { 
  
  FieldDayvalues7[1,q] <- FieldDayvalues7[1, q-1] + 0.5
  
} 
FieldDayvalues7[2,] <- FieldDayvalues7[2,]/sum(FieldDayvalues7[2,]) 


for (q in 2:length(FieldDayvalues7[2,])) {
  
  FieldDayvalues7[2,q] <- FieldDayvalues7[2, q-1] + FieldDayvalues7[2,q]
  
}

FieldDays14 <- 14
FieldDayvalues14 <- matrix(0, nrow=2, ncol=FieldDays14*2+1)
FieldDayvalues14[2,] <- c(rep(0,16), rep(1,5), rep(2,4), rep(3,2), rep(4,1), rep(5,1))

for (q in 2:length(FieldDayvalues14[1,])) { 
  
  FieldDayvalues14[1,q] <- FieldDayvalues14[1, q-1] + 0.5
  
} 
FieldDayvalues14[2,] <- FieldDayvalues14[2,]/sum(FieldDayvalues14[2,]) 
for (q in 2:length(FieldDayvalues14[2,])) {
  
  FieldDayvalues14[2,q] <- FieldDayvalues14[2, q-1] + FieldDayvalues14[2,q]
  
}


FieldDays21 <- 21
FieldDayvalues21 <- matrix(0, nrow=2, ncol=FieldDays21*2+1)
FieldDayvalues21[2,] <- c(rep(0,26), rep(1,8), rep(2,5), rep(3,2), rep(4,1), rep(5,1))
for (q in 2:length(FieldDayvalues21[1,])) { 
  
  FieldDayvalues21[1,q] <- FieldDayvalues21[1, q-1] + 0.5
  
} 
FieldDayvalues21[2,] <- FieldDayvalues21[2,]/sum(FieldDayvalues21[2,]) 
for (q in 2:length(FieldDayvalues21[2,])) {
  
  FieldDayvalues21[2,q] <- FieldDayvalues21[2, q-1] + FieldDayvalues21[2,q]
  
}



FieldDays28 <- 28
FieldDayvalues28 <- matrix(0, nrow=2, ncol=FieldDays28*2+1)
FieldDayvalues28[2,] <- c(rep(0,34), rep(1,12), rep(2,6),rep(3,2), rep(4,1), rep(5,1), rep(10,1))
for (q in 2:length(FieldDayvalues28[1,])) { 
  
  FieldDayvalues28[1,q] <- FieldDayvalues28[1, q-1] + 0.5
  
} 
FieldDayvalues28[2,] <- FieldDayvalues28[2,]/sum(FieldDayvalues28[2,]) 
for (q in 2:length(FieldDayvalues28[2,])) {
  
  FieldDayvalues28[2,q] <- FieldDayvalues28[2, q-1] + FieldDayvalues28[2,q]
  
}


FieldDays35 <- 35
FieldDayvalues35 <- matrix(0, nrow=2, ncol=FieldDays35*2+1)
FieldDayvalues35[2,] <- c(rep(0,43), rep(1,14), rep(2,7),rep(3,3), rep(4,1), rep(5,1), rep(6,1), rep(10,1))
for (q in 2:length(FieldDayvalues35[1,])) { 
  
  FieldDayvalues35[1,q] <- FieldDayvalues35[1, q-1] + 0.5
  
} 
FieldDayvalues35[2,] <- FieldDayvalues35[2,]/sum(FieldDayvalues35[2,]) 
for (q in 2:length(FieldDayvalues35[2,])) {
  
  FieldDayvalues35[2,q] <- FieldDayvalues35[2, q-1] + FieldDayvalues35[2,q]
  
}



FieldDays42 <- 42
FieldDayvalues42 <- matrix(0, nrow=2, ncol=FieldDays42*2+1)
FieldDayvalues42[2,] <- c(rep(0,51), rep(1,17), rep(2,9),rep(3,4), rep(4,1), rep(5,1), rep(6,1), rep(10,1))
for (q in 2:length(FieldDayvalues42[1,])) { 
  
  FieldDayvalues42[1,q] <- FieldDayvalues42[1, q-1] + 0.5
  
} 
FieldDayvalues42[2,] <- FieldDayvalues42[2,]/sum(FieldDayvalues42[2,]) 
for (q in 2:length(FieldDayvalues42[2,])) {
  
  FieldDayvalues42[2,q] <- FieldDayvalues42[2, q-1] + FieldDayvalues42[2,q]
  
}


FieldDays49 <- 49
FieldDayvalues49 <- matrix(0, nrow=2, ncol=FieldDays49*2+1)
FieldDayvalues49[2,] <- c(rep(0,60), rep(1,20), rep(2,10), rep(3,3), rep(4,1), rep(5,1), rep(6,1), rep(7,1), rep(8,1), rep(10,1))
for (q in 2:length(FieldDayvalues49[1,])) { 
  
  FieldDayvalues49[1,q] <- FieldDayvalues49[1, q-1] + 0.5
  
} 
FieldDayvalues49[2,] <- FieldDayvalues49[2,]/sum(FieldDayvalues49[2,]) 
for (q in 2:length(FieldDayvalues49[2,])) {
  
  FieldDayvalues49[2,q] <- FieldDayvalues49[2, q-1] + FieldDayvalues49[2,q]
  
}



FieldDays56 <- 56
FieldDayvalues56 <- matrix(0, nrow=2, ncol=FieldDays56*2+1)
FieldDayvalues56[2,] <- c(rep(0,70), rep(1,20), rep(2,11), rep(3,4), rep(4,2), rep(5,2), rep(6,1), rep(7,1), rep(8,1), rep(10,1))
for (q in 2:length(FieldDayvalues56[1,])) { 
  
  FieldDayvalues56[1,q] <- FieldDayvalues56[1, q-1] + 0.5
  
} 
FieldDayvalues56[2,] <- FieldDayvalues56[2,]/sum(FieldDayvalues56[2,]) 
for (q in 2:length(FieldDayvalues56[2,])) {
  
  FieldDayvalues56[2,q] <- FieldDayvalues56[2, q-1] + FieldDayvalues56[2,q]
  
}


FieldDays63 <- 63
FieldDayvalues63 <- matrix(0, nrow=2, ncol=FieldDays63*2+1)
FieldDayvalues63[2,] <- c(rep(0,77), rep(1,25), rep(2,13), rep(3,4), rep(4,3), rep(5,1), rep(6,1), rep(7,1), rep(8,1), rep(10,1))
for (q in 2:length(FieldDayvalues63[1,])) { 
  
  FieldDayvalues63[1,q] <- FieldDayvalues63[1, q-1] + 0.5
  
} 
FieldDayvalues63[2,] <- FieldDayvalues63[2,]/sum(FieldDayvalues63[2,])
for (q in 2:length(FieldDayvalues63[2,])) {
  
  FieldDayvalues63[2,q] <- FieldDayvalues63[2, q-1] + FieldDayvalues63[2,q]
  
}



FieldDays70 <- 70
FieldDayvalues70 <- matrix(0, nrow=2, ncol=FieldDays70*2+1)
FieldDayvalues70[2,] <- c(rep(0,85), rep(1,28), rep(2,14), rep(3,5), rep(4,4), rep(5,1), rep(6,1), rep(7,1), rep(8,1), rep(10,1))
for (q in 2:length(FieldDayvalues70[1,])) { 
  
  FieldDayvalues70[1,q] <- FieldDayvalues70[1, q-1] + 0.5
  
} 
FieldDayvalues70[2,] <- FieldDayvalues70[2,]/sum(FieldDayvalues70[2,]) 
for (q in 2:length(FieldDayvalues70[2,])) {
  
  FieldDayvalues70[2,q] <- FieldDayvalues70[2, q-1] + FieldDayvalues70[2,q]
  
}


FieldDays77 <- 77
FieldDayvalues77 <- matrix(0, nrow=2, ncol=FieldDays77*2+1)
FieldDayvalues77[2,] <- c(rep(0,93), rep(1,31), rep(2,16), rep(3,6), rep(4,3), rep(5,2), rep(6,1), rep(7,1), rep(8,1), rep(10,1))
for (q in 2:length(FieldDayvalues77[1,])) { 
  
  FieldDayvalues77[1,q] <- FieldDayvalues77[1, q-1] + 0.5
  
} 
FieldDayvalues77[2,] <- FieldDayvalues77[2,]/sum(FieldDayvalues77[2,]) 
for (q in 2:length(FieldDayvalues77[2,])) {
  
  FieldDayvalues77[2,q] <- FieldDayvalues77[2, q-1] + FieldDayvalues77[2,q]
  
}


FieldDays84 <- 84
FieldDayvalues84 <- matrix(0, nrow=2, ncol=FieldDays84*2+1)
FieldDayvalues84[2,] <- c(rep(0,102), rep(1,34), rep(2,17), rep(3,6), rep(4,4), rep(5,2), rep(6,1), rep(7,1), rep(8,1), rep(10,1))
for (q in 2:length(FieldDayvalues84[1,])) { 
  
  FieldDayvalues84[1,q] <- FieldDayvalues84[1, q-1] + 0.5
  
}
FieldDayvalues84[2,] <- FieldDayvalues84[2,]/sum(FieldDayvalues84[2,]) 
for (q in 2:length(FieldDayvalues84[2,])) {
  
  FieldDayvalues84[2,q] <- FieldDayvalues84[2, q-1] + FieldDayvalues84[2,q]
  
}


FieldDays91 <- 91
FieldDayvalues91 <- matrix(0, nrow=2, ncol=FieldDays91*2+1)
FieldDayvalues91[2,] <- c(rep(0,110), rep(1,37), rep(2,18), rep(3,6), rep(4,4), rep(5,3), rep(6,2), rep(7,1), rep(8,1), rep(10,1))

for (q in 2:length(FieldDayvalues91[1,])) {
  
  FieldDayvalues91[1,q] <- FieldDayvalues91[1, q-1] + 0.5
  
} 
FieldDayvalues91[2,] <- FieldDayvalues91[2,]/sum(FieldDayvalues91[2,]) 
for (q in 2:length(FieldDayvalues91[2,])) {
  
  FieldDayvalues91[2,q] <- FieldDayvalues91[2, q-1] + FieldDayvalues91[2,q]
  
}


FieldDays98 <- 98
FieldDayvalues98 <- matrix(0, nrow=2, ncol=FieldDays98*2+1)
FieldDayvalues98[2,] <- c(rep(0,118), rep(1,40), rep(2,20), rep(3,7), rep(4,5), rep(5,3), rep(6,1), rep(7,1), rep(8,1), rep(10,1))
for (q in 2:length(FieldDayvalues98[1,])) { 
  
  FieldDayvalues98[1,q] <- FieldDayvalues98[1, q-1] + 0.5
  
}
FieldDayvalues98[2,] <- FieldDayvalues98[2,]/sum(FieldDayvalues98[2,]) 
for (q in 2:length(FieldDayvalues98[2,])) {
  
  FieldDayvalues98[2,q] <- FieldDayvalues98[2, q-1] + FieldDayvalues98[2,q]
  
}


FieldDays105 <- 105
FieldDayvalues105 <- matrix(0, nrow=2, ncol=FieldDays105*2+1)
FieldDayvalues105[2,] <- c(rep(0,127), rep(1,42), rep(2,21), rep(3,8), rep(4,6), rep(5,3), rep(6,1), rep(7,1), rep(8,1), rep(10,1))
for (q in 2:length(FieldDayvalues105[1,])) { 
  FieldDayvalues105[1,q] <- FieldDayvalues105[1, q-1] + 0.5
  
} 
FieldDayvalues105[2,] <- FieldDayvalues105[2,]/sum(FieldDayvalues105[2,])
for (q in 2:length(FieldDayvalues105[2,])) {
  
  FieldDayvalues105[2,q] <- FieldDayvalues105[2, q-1] + FieldDayvalues105[2,q]
  
}


FieldDays111 <- 111
FieldDayvalues111 <- matrix(0, nrow=2, ncol=FieldDays111*2+1)
FieldDayvalues111[2,] <- c(rep(0,134), rep(1,45), rep(2,22), rep(3,9), rep(4,6), rep(5,3), rep(6,1), rep(7,1), rep(8,1), rep(10,1))
for (q in 2:length(FieldDayvalues111[1,])) { 
  
  FieldDayvalues111[1,q] <- FieldDayvalues111[1, q-1] + 0.5
  
} 
FieldDayvalues111[2,] <- FieldDayvalues111[2,]/sum(FieldDayvalues111[2,])
for (q in 2:length(FieldDayvalues111[2,])) {
  
  FieldDayvalues111[2,q] <- FieldDayvalues111[2, q-1] + FieldDayvalues111[2,q]
  
}

FieldDays150 <- 150
FieldDayvalues150 <- matrix(0, nrow=2, ncol=FieldDays150*2+1)
FieldDayvalues150[2,] <- c(rep(0,180), rep(1,60), rep(2,30), rep(3,12), rep(4,7), rep(5,4), rep(6,2), rep(7,2), rep(8,2), rep(10,2))
for (q in 2:length(FieldDayvalues150[1,])) { 
  
  FieldDayvalues150[1,q] <- FieldDayvalues150[1, q-1] + 0.5
  
} 
FieldDayvalues150[2,] <- FieldDayvalues150[2,]/sum(FieldDayvalues150[2,]) 
for (q in 2:length(FieldDayvalues150[2,])) {
  
  FieldDayvalues150[2,q] <- FieldDayvalues150[2, q-1] + FieldDayvalues150[2,q]
  
}


FieldDays222 <- 222
FieldDayvalues222 <- matrix(0, nrow=2, ncol=FieldDays222*2+1)
FieldDayvalues222[2,] <- c(rep(0,268), rep(1,90), rep(2,47), rep(3,18), rep(4,10), rep(5,5), rep(6,3), rep(7,2), rep(8,1), rep(10,1))
for (q in 2:length(FieldDayvalues222[1,])) { 
  
  FieldDayvalues222[1,q] <- FieldDayvalues222[1, q-1] + 0.5
  
} 
FieldDayvalues222[2,] <- FieldDayvalues222[2,]/sum(FieldDayvalues222[2,]) 
for (q in 2:length(FieldDayvalues222[2,])) {
  
  FieldDayvalues222[2,q] <- FieldDayvalues222[2, q-1] + FieldDayvalues222[2,q]
  
}


##########################################################################

# Set up lists of the mallee scenario over the varying time steps. 

Scenario_One <- list(mu, FieldDayvalues7, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=7, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Two <- list(mu, FieldDayvalues14, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=14, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Three <- list(mu, FieldDayvalues21, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=21, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Four <- list(mu, FieldDayvalues28, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=28, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Five <- list(mu, FieldDayvalues35, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=35, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Six <- list(mu, FieldDayvalues42, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=42, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Seven <- list(mu, FieldDayvalues49, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=49, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Eight <- list(mu, FieldDayvalues56, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=56, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Nine <- list(mu, FieldDayvalues63, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=63, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Ten <- list(mu, FieldDayvalues70, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=70, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Eleven <- list(mu, FieldDayvalues77, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=77, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Twelve <- list(mu, FieldDayvalues84, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=84, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Thirteen <- list(mu, FieldDayvalues91, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=91, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Fourteen <- list(mu, FieldDayvalues98, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=98, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Fifteen <- list(mu, FieldDayvalues105, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=105, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Sixteen <- list(mu, FieldDayvalues111, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=111, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Seventeen <- list(mu, FieldDayvalues150, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=150, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenario_Eighteen <- list(mu, FieldDayvalues222, FieldMinsDay, HometoSiteTravel, Maxtime, TravelbwSites, PriorityList, N.fieldtrips, TotalFieldDays=222, TSFvalues, N.reps, setuptime, mt, HrsPerDay, p_matrix, Occupancy_meta_matrix, sd.obs, N.TSF)
Scenarios_all <- list(Scenario_One, Scenario_Two, Scenario_Three, Scenario_Four, Scenario_Five, Scenario_Six, Scenario_Seven, Scenario_Eight,Scenario_Nine,Scenario_Ten,Scenario_Eleven,Scenario_Twelve,Scenario_Thirteen,Scenario_Fourteen,Scenario_Fifteen,Scenario_Sixteen,Scenario_Seventeen,Scenario_Eighteen)


##################################

# Fieldtrip function - collecting between 2 and 5 individuals per species

Fieldtrip <- function(Scenario) { 
  
  as.numeric(Sys.time())-> g; set.seed((g - floor(g)) * 1e8 -> seed); print(seed)
  
  N.species <- dim(mu)[1]
  
  species.ID <- c(rep(0,20))
  
  Rando <- runif(1, 0, 1)
  
  RealIndex <- min(which(Rando <= Scenario[[2]][2,])) 
  
  RealFieldDays <- Scenario[[2]][1, RealIndex] 
  
  RealFieldDaysMins <- RealFieldDays * Scenario[[3]] 
  
  TotalSiteTime <- RealFieldDaysMins - Scenario[[4]] * RealFieldDays  
  
  TotalSites <- round(TotalSiteTime/ (Scenario[[5]] + Scenario[[6]]))
  
  rm(Rando, RealIndex, RealFieldDays, RealFieldDaysMins, TotalSiteTime) 
  
  
  FittingData <- array(NA, c(40000, 4))
  colnames(FittingData) <- c('Sp', 'Height', 'TSF', 'Rep')
  
  N.plants <- 0
  
  Findsthrutime <- array(NA, c(40000, 6, N.reps, N.TSF))
  colnames(Findsthrutime) <- c('SpeciesID', ' ', ' ', ' ', 'Time(mins)','Height(cm)')
  
  SummaryData <- array(NA, c(N.species, 6, N.reps, N.TSF))
  colnames(SummaryData) <- c('indivs_counts', 'Presence', 'Occup', 'Searchable', 'mean H(cm)', 'TSF')
  
  
  for (t in 1:TotalSites) {
    
    TSF <- Scenario[[7]][t, 1]
    
    TSFIndex <- which(Scenario[[10]] == TSF)
    
    Rep <- Scenario[[7]][t, 2]
    
    matrix_ind_counts <- matrix(0, N.species, 1)
    matrix_test <- matrix(NA, 40000, 6)   
    colnames(matrix_test) <- c('ChosenSpecies','indivID',  'TSF',  'Rep', 'time_mins','H')
    
    Presence <- rbinom(N.species, 1, Scenario[[16]][,TSFIndex])
    
    Searchable <- Presence
    
    timeatsite <- 0
    
    setuptime <- Scenario[[12]] 
    timeatsite <- timeatsite + setuptime
    
    datarowcounter <- 0
    
    PossDetectiontime <- rexp(N.species, 1/Scenario[[15]][,TSFIndex])
    PossMinDetectiontime <- min(PossDetectiontime[Searchable==1]) 
    
    
    while((1 - prod(matrix_ind_counts[Presence==1] >= 2 )) & (timeatsite + PossMinDetectiontime + Scenario[[13]] < Scenario[[5]])) {  # this is where the minimum number of species is defined
      
      Detectiontime <- PossMinDetectiontime
      
      ChosenSpecies <- which(Detectiontime == PossDetectiontime)
      species.ID[ChosenSpecies] <- 1
      
      datarowcounter <- datarowcounter + 1
      
      N.plants <- N.plants + 1
      
      matrix_ind_counts[ChosenSpecies, 1] <- matrix_ind_counts[ChosenSpecies, 1] + 1
      
      matrix_test[datarowcounter, 1] <- ChosenSpecies
      matrix_test[datarowcounter, 2] <- N.plants
      matrix_test[datarowcounter, 3] <-  TSF
      matrix_test[datarowcounter, 4] <-  Rep
      
      k <- log(mu[ChosenSpecies, TSFIndex]) - Scenario[[17]]^2
      
      H <- rlnorm(1, k, Scenario[[17]])
      
      FittingData[N.plants, 1] <- ChosenSpecies
      FittingData[N.plants, 2] <- H
      FittingData[N.plants, 3] <- TSFvalues[TSFIndex]
      FittingData[N.plants, 4] <- Rep
      
      
      matrix_test[datarowcounter, 6] <- H
      
      timeatsite <- timeatsite + Detectiontime + Scenario[[13]]
      
      matrix_test[datarowcounter, 5] <- timeatsite
      
      
      if  (matrix_ind_counts[ChosenSpecies, 1] == 5 ) { # this is where the maximum number of species is defined
        
        Searchable[ChosenSpecies] <- 0
        
      }
      
      PossDetectiontime <- rexp(N.species, 1/ Scenario[[15]][,TSFIndex])
      
      PossMinDetectiontime <- min(PossDetectiontime[Searchable == 1]) 
      
      rm(k,H, ChosenSpecies)
      
    }
    
    rm(timeatsite,setuptime,datarowcounter, PossDetectiontime, PossMinDetectiontime, Detectiontime)
    
    Findsthrutime[ ,  , Rep, TSFIndex] <- matrix_test
    
    SummaryData[ , 1, Rep, TSFIndex] <- matrix_ind_counts
    SummaryData[ , 2, Rep, TSFIndex] <- Presence     
    SummaryData[ , 3, Rep, TSFIndex] <- Scenario[[16]][,TSFIndex]
    SummaryData[ , 4, Rep, TSFIndex] <- Searchable
    SummaryData[ , 5, Rep, TSFIndex] <- Scenario[[1]][,6]
    SummaryData[ , 6, Rep, TSFIndex] <- TSF
    
    rm(TSF, TSFIndex, Rep, matrix_ind_counts, matrix_test, Presence, Searchable)
    
  }
  
  rm(TotalSites)
  
  return(list(FittingData, N.plants, Findsthrutime, SummaryData, species.ID))
  
  rm(FittingData, N.plants, Findsthrutime, SummaryData, species.ID)  
  
}




#########################################################################

# simulate 100 fieldtrips for dataset simulation - perhaps model less?

N.species <- dim(mu)[1]
N.scen <- 17
N.fieldtrips <- 100
All_FittingData <- array(NA, c(40000, 4, N.scen, N.fieldtrips))
#colnames(All_FittingData) <- c('Sp', 'Height', 'TSF', 'Reps')
Sample_sizes <- array(NA, c(N.scen, N.fieldtrips))
Sample_species <- array(NA, c(1,N.scen, N.fieldtrips))
All_Findsthrutime <- array(NA, c(40000, 6, N.reps, N.TSF, N.scen, N.fieldtrips))
#colnames(All_Findsthrutime) <- c('ChosenSpecies','N.plants',  'TSF',  'Rep', 'time_mins','H')
All_SummaryData <- array(NA, c(N.species, 6, N.reps, N.TSF, N.scen, N.fieldtrips))
#colnames(All_SummaryData) <- c('indivs_counts', 'Presence', 'Occup', 'Searchable', 'mean H(cm)', 'TSF')
All_species.ID <- array(NA, c(N.species, N.scen, N.fieldtrips))


##########################################################################


# Running these models over many fieldtrips many scenarios may take a long time.

for (i in 1:N.scen) {
  
  for (f in 1:N.fieldtrips) {
    
    temp.list <- Fieldtrip(Scenarios_all[[i]])
    All_FittingData[ , ,i,f] <-  temp.list[[1]]
    Sample_sizes[i,f] <- temp.list[[2]]
    All_Findsthrutime[, , , ,i,f] <- temp.list[[3]]
    All_SummaryData[, , , ,i,f] <- temp.list[[4]]
    All_species.ID[ ,i,f] <- temp.list[[5]]
    
    Sample_species[ , i,f] <- length(na.omit(unique(temp.list[[1]][,1])))
    temp.species <- as.numeric(length(unique(All_FittingData[ ,1 ,i,f])))
    
    #All_Results[1:(Sample_species[ ,i,f] * 3 + 7), ,i,f] <- Fitmodel(All_FittingData[1:Sample_sizes[i,f], ,i,f])
    
    rm(temp.list)
    
    
  }
  
  
}


############################################################################



