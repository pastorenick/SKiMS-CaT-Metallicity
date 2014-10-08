#### README FILE FOR THE CORRECT ORDER OF THE OPERATIONS
#### v 1.2 - September 2014
####

1) CaTindexMeasurement (or _modG2 to run it in parallel)

2) fittingOffset (to measure the offset respect to SAURON value)
2b) fittingMetOffsetN4365 (to measure the offset respect ATLAS3d value

3) MetallicityConversion (CaT->metallicity and sigma correction)

4) CleanOutput 

5) KrigingMapping (v2 or in parallel)

6) radialProfilesAndGradients

#### If the datasets have been updated, for some galaxies it could be 
#### necessary to update the offsets with SAURON/ATLAS3d.
#### Run fittingMetOffset.py and fittingMetOffsetN4365.py after 4) and 
#### then again 3), 4), 5) and 6)