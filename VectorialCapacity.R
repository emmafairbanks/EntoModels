library(ggplot2)
library(matrixStats)

set.seed(135)

# load climate data
Temp = c(rnorm(50, 15, 5), rnorm(90, 20, 5), rnorm(90, 25, 5),
         rnorm(90, 20, 5), rnorm(45, 25, 5)) # example data
ts = 1:length(Temp)

VectorCompetence_hv = 0.52
VectorCompetence_vh = 0.77
BloodIndex = 0.56
BitesPerCycle = 1
      
#### Vector control
Coverage = 0
Usage = 0
BiteReduction = 0
RelativePreKill = 0
RelativeDisarm = 0
PostKill = 0
        
if(RelativePreKill>BiteReduction){
  print('Warning: The relative rate of preprandial killing compared to
          the biting rate without the tool is larger than the reduction
          in the rate of biting with the tool. Check this is the correct 
          assumption for the tool modelled.')
}
        
#### Climate dependent rates
OviRate =  0.0002*Temp*(Temp - 3.7)*(41.9 - Temp)^0.37
OviRate[which(Temp<3.7)] = 0
MortRate =  0.0009*exp(0.16*Temp)
IncubRate =  0.017*(Temp - 12.6)
IncubRate[which(Temp<12.6)] = 0
        
ProbStopHostSeeking = 1 - exp(-(BitesPerCycle*OviRate*((
    (Coverage*(1-Usage*(BiteReduction - RelativePreKill - RelativeDisarm)) +
        (1-Coverage))*BloodIndex) + 
        (1 -BloodIndex))))

ProbAlive = function(t, tau, MortRate){
  return(exp(-sum(MortRate[t:tau])))
}

ProbInf = function(t, tau, IncubRate){
  k = 1
  scale = 1
  return(pgamma(sum(IncubRate[(t+1):tau]),
                shape = k, scale = scale/k))
}

ProbInfectiousBite_hv = function(t, BloodIndex, Coverage, Usage,
                                 BiteReduction,
                                 RelativePreKill, RelativeDisarm,
                                 VectorCompetence_hv, ProbStopHostSeeking){
  return(VectorCompetence_hv*ProbStopHostSeeking[t]*
           BloodIndex*(Coverage*(1 - Usage*BiteReduction) + (1-Coverage))/
           (BloodIndex*(Coverage*(1 - Usage*(BiteReduction - RelativePreKill - RelativeDisarm)) +(1-Coverage)) + (1 - BloodIndex)))
}

ProbInfectiousBite_vh = function(BloodIndex, Coverage, Usage,
                                 BiteReduction,
                                 RelativePreKill, RelativeDisarm, 
                                 VectorCompetence_vh,
                                 tau, IncubRate,
                                 ProbStopHostSeeking){
  return(ProbStopHostSeeking[tau]*
           (BloodIndex*((Coverage*(1 - Usage*BiteReduction)) + (1-Coverage))/
           (BloodIndex*(Coverage*(1 - Usage*(BiteReduction - RelativePreKill - RelativeDisarm)) + (1-Coverage)) + (1 - BloodIndex)))*
           VectorCompetence_vh*ProbInf(t, tau, IncubRate))
}

ProbSurviveTool = function(t, tau, RatePreKill, RatePostKill){
  k = 1
  mean = 1
  ProbPreKilled = pgamma(sum(RatePreKill[t:tau]),
                         shape = k, scale = mean/k)
  ProbPostKilled = pgamma(sum(RatePostKill[(t+1):tau]),
                          shape = k, scale = mean/k)
  return((1-ProbPreKilled)*(1-ProbPostKilled))
}

RatePreKill = data.frame(matrix(NA, nrow = 1, ncol = length(ts)))
RatePostKill = data.frame(matrix(NA, nrow = 1, ncol = length(ts)))
for(t in ts){
  RatePreKill[t] = as.numeric(ProbStopHostSeeking[t]*
                     BloodIndex*Coverage*Usage*RelativePreKill/
                     (BloodIndex*(Coverage*(1 - Usage*(BiteReduction - RelativePreKill - RelativeDisarm)) +
                       (1-Coverage)) + (1 - BloodIndex)))
}
for(t in ts){
  RatePostKill[t] = as.numeric(ProbStopHostSeeking[t]*
                      (BloodIndex*Coverage*(1 - Usage*BiteReduction)  + (1-Coverage))*PostKill/
                      (BloodIndex*(Coverage*(1 - Usage*(BiteReduction - RelativePreKill - RelativeDisarm)) + (1-Coverage)) + (1 - BloodIndex)))
}
      
#####BASELINE####
      
VC_baseline = data.frame(matrix(NA, nrow = 0, ncol = 2))
names(VC_baseline) = c('Time', 'Value')
        
for(t in ts){
  Output1 = ProbInfectiousBite_hv(t, BloodIndex, Coverage, Usage,
                                  BiteReduction,
                                  RelativePreKill, RelativeDisarm,
                                  VectorCompetence_hv, ProbStopHostSeeking)*
              ProbAlive(t,t,MortRate)
  Output2 = 0
  for(tau in (t+1):(max(ts)-1)){
    if(tau<max(ts) & t<tau){
    Output2 = Output2 + 
      ProbInfectiousBite_vh(BloodIndex, Coverage, Usage,
                            BiteReduction,
                            RelativePreKill, RelativeDisarm, 
                            VectorCompetence_vh,
                            tau, IncubRate,
                            ProbStopHostSeeking)*
        ProbAlive(t,tau,MortRate)
    }
  }
  VC_baseline[dim(VC_baseline)[1] + 1,] = c(t, Output1*Output2)
}

ggplot() +
  geom_line(data = VC_baseline, aes(x = Time, y = Value), color = "blue") +
  theme_minimal()
