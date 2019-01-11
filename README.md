# Surface-water availability governs survival of an amphibian in arid mountain streams

### Erin R. Zylstra, Don E. Swann, and Robert J. Steidl

### Freshwater Biology 64:164-174 

### Manuscript DOI: [10.1111/fwb.13204](https://onlinelibrary.wiley.com/doi/10.1111/fwb.13204)

### Please contact the first author for questions about the code or data: Erin R. Zylstra (ezylstra@email.arizona.edu)
_____________________________________________________________________________

## Abstract:
1. Declines in amphibian populations worldwide have been linked to multiple factors, including recent changes in climate. Changes in precipitation, for example, can alter the availability of aquatic resources that are required by many amphibians for successful reproduction and larval development.  
2. Changes in climate have the potential to affect other important demographic parameters, particularly for species that are active year-round. Therefore, we studied survival of post-metamorphic (snout-vent length ≥50 mm) lowland leopard frogs (*Lithobates yavapaiensis*) in arid mountain canyons, where surface water is limited and its availability can vary markedly within and among years.  
3. Between 2013 and 2015, we surveyed frogs 33-74 times in each of six stream reaches distributed across two watersheds in southern Arizona and used capture-recapture methods based on in-situ photographs to identify individuals. We used Cormack-Jolly-Seber models to explore how surface-water availability, weather, and vegetation influenced seasonal variation in apparent survival of post-metamorphic individuals.  
4. Overall, mean annual apparent survival in this dynamic, arid environment was low (0.11, 95% CI = 0.07-0.14). Survival varied with ambient temperature, dew point, perimeter groundcover, and year, but especially with changes in surface-water availability. When water levels were at or near 100% of maximum pool depths, mean monthly survival was high (≥0.88); when water levels were at 50%, survival decreased modestly (0.81) and when at 20%, survival dropped sharply (0.36).  
5. A decrease in survival of post-metamorphic frogs in response to severe drought almost certainly contributed to the extirpation of frogs from one watershed in 2015. We anticipate that predicted increases in frequency and severity of drought will decrease the probability that lowland leopard frogs persist in this region over the long-term, as droughts are expected to increase local extirpations and limit the ability of individuals to disperse through an increasingly arid landscape.  

## Main files
We ran analyses in R using the RMark package (which calls Program MARK).  The data and code needed to run analyses described in the paper are contained in the following three files:

1. [CensoredData_InputFile.inp](CensoredData_InputFile.inp): A text file formatted for Program MARK that contains capture histories (after censoring the initial observation of each individual) for post-metamorphic lowland leopard frogs in the Rincon Mountains between 2013 and 2015.
2. [CensoredData_Covariates.csv](CensoredData_Covariates.csv): Covariates associated with each survey occasion or survey interval (after censoring the initial observation of each individual). The “duration” column lists the length of each interval between survey occasions, in days. Covariates are described in detail in the R code.
3. [CensoredData_RMark.R](CensoredData_RMark.R): R code to import and format data for analyses in MARK, run Cormack-Jolly-Seber models, and summarize and plot results.

## Supplementary files
In the paper, we based inferences on capture-recapture data from frogs that were photographed in-situ. We identified individuals based on their unique spot patterns. Because we were not able to photograph the entire spot pattern for each individual, it was possible that the number of encounter histories exceeded the number of individuals observed. To reduce potential bias associated with these “false rejection” errors, we censored the first observation of each individual. To assess how censoring methods may have affected parameter estimates, we estimated survival of leopard frogs based on capture-histories that included all observations (uncensored) for 1) all frogs, and 2) only those frogs whose spot patterns were mapped completely. More information about these comparisons can be found in [Appendix S1](Zylstra_etal_2019_AppendixS1.docx). The data files and code needed to run these analyses are described below:

1. [UncensoredData_AllFrogs_InputFile.inp](UncensoredData/UncensoredData_AllFrogs_InputFile.inp): A text file formatted for Program MARK that contains capture histories with all observations of all post-metamorphic lowland leopard frogs in the Rincon Mountains between 2013 and 2015.
2. [UncensoredData_CompleteSpotMaps_InputFile.inp](UncensoredData/UncensoredData_CompleteSpotMaps_InputFile.inp): A text file formatted for Program MARK that contains capture histories with all observations of only those leopard frogs for which we were able to document their entire spot pattern.
3. [UncensoredData_Covariates.csv](UncensoredData_Covariates.csv): Covariates associated with each survey occasion or survey interval. The “duration” column lists the length of each interval between survey occasions, in days. Covariates are described in detail in the R code.
4. [UncensoredData_AllFrogs_RMark.R](UncensoredData_AllFrogs_RMark.R): R code to import and format data for analyses in MARK, run Cormack-Jolly-Seber models, and summarize results for all observations of all frogs.
5. [UncensoredData_CompleteSpotMaps_RMark.R](UncensoredData_CompleteSpotMaps_RMark.R): R code to import and format data for analyses in MARK, run Cormack-Jolly-Seber models, and summarize results for all observations of frogs with complete spot maps.
