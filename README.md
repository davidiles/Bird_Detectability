# Bird_Detectability
This repo contains code to test approaches for landbird detectability modeling.
 
Over-arching goal is to critically evaluate 'goodness-of-fit' of standard distance and time-removal modeling, and the sensitivity of model outputs to violations of model assumptions.  For example, 

Secondary goal: explore alternative parameterizations that may provide better fit.

## Investigation 1: Do cue rate estimates differ when estimated from ARUs and humans




**Notes**

*Analysis ideas*
 1) In recordings, is the time of first detection associated with sound pressure level (volume)?
 2) From acoustic arrays, does time of first detection depend on distance of the bird?
 3) From acoustic recordings, can we estimate cue intervals (using more than just the first detection)?
 4) From acoustic recordings, do cue rate estimates based on first detection correspond to estimates from first two (or more) cues?
 5) Does it really make sense to talk about "cue rate"?  This doesn't take into account the duration of each cue.
 6) Are cues distributed as an exponential?

*Meeting with Steve VW on 2024-01-29*
- confirmed that ARUs could provide better estimates of cue rates
- farnsworth mixture distribution


*2024-01-24*
1)	Elly found that ARUs tend to result in much higher cue rate estimates, compared to human point counts.  This could be due to humans being unable to bin detections in time properly (e.g., being overwhelmed by a busy dawn chorus).
a.	Action: Check whether this occurs in paired human/ARU data, based on Steve