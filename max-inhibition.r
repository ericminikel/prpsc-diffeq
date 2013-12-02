########### what percent of possible degradation is achieved by 1 day, 5 days and 6 days?

maxint = .Machine$integer.max # longest possible time in future to model

# transcription inhibition
r = l0 * 0.50 # 
l = l0 * 1.00 # 
a = m0 * 1.00
m = m0 * 1.00 
b = v0 * 1.00
v = v0 * 1.00
# percent of maximum inhibition achived after various time periods
(1 - prpc(1*24))/(1-prpc(maxint))   # .76
(1 - prpsc(5*24))/(1-prpsc(maxint)) # .90
(1 - prpsc(6*24))/(1-prpsc(maxint)) # .94
