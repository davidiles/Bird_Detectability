# Bird_Detectability
 Development of new approaches for landbird detectability modeling from human point counts.
 
Goal is to critically evaluate 'goodness-of-fit' of standard distance and time-removal modeling, and explore alternative parameterizations that may provide better fit.


**Notes**

*2024-01-24*
1)	Elly found that ARUs tend to result in much higher cue rate estimates, compared to human point counts.  This could be due to humans being unable to bin detections in time properly (e.g., being overwhelmed by a busy dawn chorus).
a.	Action: Check whether this occurs in paired human/ARU data, based on Steve’s SK dataset

2)	Heterogeneity in bird behaviour would result in biased estimates of offsets and density.  E.g., if some birds call much less frequently than others, or are much less detectable.  Could also occur if there is unmodeled heterogeneity (e.g., in habitat)
a.	Action: Evaluate if standard QPAD (which does not attempt to estimate a mixture of bird types) fits the data adequately, and if not, whether a mixture model is more appropriate.  Do mixture models fit better than non-mixture models?

3)	Joint vs independent – are first cues always detected?


*2023-11-22*
- "current" analysis is in Analysis folder
- goodness of fit of standard models seems poor (many extra birds are detected in the first minute compared to predictions, and then subsequently less)

- model_fns_new.R contains latest attempts at fitting mixture models, where detected birds can be of two types (fast/slow cue rates, and high/low EDRs)

- models are fit to empirical data

- appears to be difficult to fit the models; likelihood surface very rough, and optimizers arrive at local optima, so need to be run multiple times with many initial values.

- next steps may be to ask David Hope to implement a revised version in CPP.  However, need to convince ourselves that this would be worth it.