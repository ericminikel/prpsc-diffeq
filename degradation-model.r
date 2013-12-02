# Eric Minikel
# CureFFI.org
# 2013-12-02
# Comparison of numerically simulated and analytically derived PrP degradation models

# half lives from literature
mrna_thalf  = 7  # Pfeiffer 1993  http://www.ncbi.nlm.nih.gov/pubmed/8095862
prpc_thalf  = 5 # Borcheldt 1990 http://www.ncbi.nlm.nih.gov/pubmed/1968466/
prpsc_thalf = 30 # Peretz 2001    http://www.ncbi.nlm.nih.gov/pubmed/11507642

# values at initial steady state
l0 = log(2) / mrna_thalf  # lambda_0
m0 = log(2) / prpc_thalf  # mu_0
v0 = log(2) / prpsc_thalf # nu_0

hrs = 8*24 # number of hours to model and plot. here, set to 8 days

# values of parameters, which you can tweak to explore the effects of diff. compounds
r = l0 * 0.50 # try setting to r0 * 0.50 to reduce transcription rate by half
l = l0 * 1.00 # try setting to l0 * 2.00 to double the mRNA degradation rate
a = m0 * 1.00
m = m0 * 1.00 
b = v0 * 1.00
v = v0 * 1.00

# initial steady state
R0 = 1
S0 = 1
T0 = 1

########### simulation #################

mrna_produced = numeric(hrs)
mrna_degraded = numeric(hrs)
mrna_present = numeric(hrs)
mrna_present[1] = R0
mrna_produced = rep(r, hrs)
for (i in 2:hrs) {
   mrna_degraded[i] = l * mrna_present[i-1]
   mrna_present[i] = mrna_present[i-1] + mrna_produced[i] - mrna_degraded[i]
}

prpc_produced = numeric(hrs)
prpc_degraded = numeric(hrs)
prpc_present = numeric(hrs)
prpc_present[1] = S0
prpc_produced = a * mrna_present
for (i in 2:hrs) {
   prpc_degraded[i] = m * prpc_present[i-1]
   prpc_present[i] = prpc_present[i-1] + prpc_produced[i] - prpc_degraded[i]
}

prpsc_produced = numeric(hrs)
prpsc_degraded = numeric(hrs)
prpsc_present = numeric(hrs)
prpsc_present[1] = T0
prpsc_produced = b * prpc_present
for (i in 2:hrs) {
   prpsc_degraded[i] = v * prpsc_present[i-1]
   prpsc_present[i] = prpsc_present[i-1] + prpsc_produced[i] - prpsc_degraded[i]
}

plottitle = 'Degradation of PrP mRNA, PrPC and PrPSc over time'

png('degradation.model.png',width=600,height=400)

# plot mRNA, PrPC and PrPSc values from simulation
plot(1:hrs,mrna_present,pch=1,col='red',ylim=c(0,1),ylab='normalized level',xlab='hours of incubation',main=plottitle,cex.main=.8,xaxt='n')
points(1:hrs,prpc_present,pch=1,col='orange')
points(1:hrs,prpsc_present,pch=1,col='purple')

######### analytical solutions #############

# to compare to simulated solutions, plot 1:hrs vs. 0:(hrs-1) because
# in simulation, knockdown begins at t = 1 hr, while in analytical
# solution, knockdown begins at t = 0 hr

mrna = function(t) {
  return ( r/l + (R0 - r/l)*exp(-l*t) )
}
points(1:hrs, mrna(0:(hrs-1)), type='l', lwd=2, col='red')

prpc = function(t) {
  return ( a*r/(m*l) + a*((R0 - r/l)/(m-l))*exp(-l*t) + (S0 - a*r/(m*l) - a*(R0 - r/l)/(m-l))*exp(-m*t) )
}
points(1:hrs, prpc(0:(hrs-1)), type='l', lwd=2, col='orange')

prpsc = function(t) {
  return ( b*a*r/(v*m*l) + b*a*(R0-r/l)/((v-l)*(m-l))*exp(-l*t) + b*(S0 - a*r/(m*l) - a*(R0-r/l)/(m-l))/(v-m)*exp(-m*t) + 
         (T0 - b*a*r/(v*m*l) - b*a*(R0-r/l)/((v-l)*(m-l)) - b*(S0 - a*r/(m*l) - a*(R0-r/l)/(m-l))/(v-m))*exp(-v*t) )
}
points(1:hrs, prpsc(0:(hrs-1)), type='l', lwd=2, col='purple')

legend('topright',c('PrPSc simulated','PrPC simulated','mRNA simulated','PrPSc analytical','PrPC analytical','mRNA analytical'),
       col=c('purple','orange','red','purple','orange','red'),
       pch=c(1,1,1,NA,NA,NA), lwd=c(NA,NA,NA,2,2,2))

dev.off()
