# This R script saves the matching list between events and their stochastic event sets as an .RData file.

# For matching events with SES:
load(file = "D:/Documents/Proj_Damage_Spatial_Corr/Data/master.events_OQ.RData")
# master.events.

event_ses_list <- master.events

save(event_ses_list, file = "F:/Documents/RPackages/fragilitycurves/data/event_ses_list.RData")
