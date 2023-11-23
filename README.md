# Bird_Detectability
 Development of new approaches for landbird detectability modeling from human point counts.
 
 - Goal is to critically evaluate 'goodness-of-fit' of standard distance and time-removal modeling, and explore alternative parameterizations that may provide better fit.


**Notes**
*2023-11-22*
- "current" analysis is in Analysis folder
- goodness of fit of standard models seems poor (many extra birds are detected in the first minute compared to predictions, and then subsequently less)

- model_fns_new.R contains latest attempts at fitting mixture models, where detected birds can be of two types (fast/slow cue rates, and high/low EDRs)

- models are fit to empirical data

- appears to be difficult to fit the models; likelihood surface very rough, and optimizers arrive at local optima, so need to be run multiple times with many initial values.

- next steps may be to ask David Hope to implement a revised version in CPP.  However, need to convince ourselves that this would be worth it.