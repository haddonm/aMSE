## aMSE Recent Activity

* Added a `NEWS.md` file to track/log the development of the package.

* Currently can generate a full region, characterize its productivity, and if run with no catch or variability maintains an equilibrium, even in the presence of larval drift.

* 16-04-2020 worked on _setuphtml_ and _addfilename_ to simplify the use of make-html when adding a plot or table to the results for a simulation.

* 17-04-2020 0.0.0.8920: Cleaned up the make_html code to produce valid HTML5 and CSS3, separated the CSS code into a separate file that is linked into the different html files. Cleaned up the imports and fixed a bunch of the examples, so that there are now no errors, no warnings, and no notes. Although the examples for the read[ctrl][data][region1] functions all currently have a 'dontrun' status. An issue to fix. 

* 01-05-2020 0.0.0.8600: Continued with the development of Running_aMSE.Rmd. This has the great advantage of forcing me to clarify the text and the relationships between the R functions and the tasks they do. It even suggests ways to improve on the names of the functions. Naming things continues to be one of the harder things (if I want to keep these things meaningful, nothing wrong with a1, a2, x1, x2, etc). Have discovered how to include some CSS code to customize the vignette format. But, more importantly, I have improved the text, changed _addfilename_ to _logfilename_ and automated recording the file type (which is used by make_html to process the file correctly on the local webpage).

* 10-05-2020: 0.0.0.8500 Continued the development of Running_aMSE.Rmd. Clarified the intent of this document. Cleaned up the makehtml code and added the dodepletion section. Now need to implement at least one HCR function and associated file.

* 11-05-2020 0.0.0.8450 Put the results generation into internal functions so their production can be automated as new ones are developed. Ensured all internal functions were operating consistently with regard the various name changes I have made. e.g. fixed the emergence calculations, which was really messing with the results! Set up a better worked example in the readme in the github version so it now produces plots and tables and can generate the local website should the user wish.

* 12-05-2020 0.0.0.8400 Simplified the intended directory structure to a single results directory _resdir_. Cleaned out old, now unused code, and tidied the help files.

* 15-05-2020 0.0.0.8350 Started including estimates relating to SMUs as well as populations and whole of region. This entailed making changes to the design of regionD, to ensure internal consistency and ease of use. Began adding details of the fleet dynamics to Running-aMSE.Rmd.

* 17-05-2020 0.0.0.8300 Some big changes. I have now put a larval dispersal movement matrix into the global glb object. This is now applied to the predicted recruitments before they settle so that the equilibrium and dynamic conditions now take this potential movement into account within the usual dynamics. All other examples and functions now have been modified to account for this. In addition, I have included a src directory to contain the Rcpp and RcppArmadillo code, which I hoped to use to speed things up. However, potential speedups with inverting matrices and multiplying vectors by matrices appear to provide no benefit or make things worse! Further explorations are needed. However, I have introduced two small functions to replace the use of sapply when accessing components of regionC, and this has sped thing up by about 10%. A strange error where the example code from the vignettes and their output are put into files in the aMSE.check directory and then that is all deleted but 'check' finds these files and complains, even though it put them there? Again, further investigation required. 

* 18-05-2020 0.0.0.8250 Now implemented the movement matrix into oneyearC and oneyearcat, which use catches rather than harvest rates to impose a fishery onto the simulation. The larval dispersal matrix is now applied after the calculation of the expected recruitment levels using _oneyearrec_. That way if no fishing is applied everything remains static (when there is no variation on anything). But now, we can apply catches rather than harvest rates so all harvest control rule outputs (in catches) will be valid, even when they fail miserably. The check note about unknown directories is continuing and I think it is due to my explicit use of tempdir in a few examples, especially in makehtml.R. I will keep looking into it, even though it does not affect functionality in any way.

* 22-05-2020 0.0.0.8200 Added a fishery dataset from East coast block 13 that includes cpue, catch, and effort. This is for use with exploring the operation of the estimators of the performance measures. I have also now included two new get functions (among others) one to calculate the grad1 values and then one to calculate the Tasmanian harvest strategy gradient1 scores (as yet unweighted).

* 25-05-2020 0.0.0.8150 Added getgrad4 and modified getgrad1 and getscore1 on the path to having an experimental harvest control rule that approximates the Tasmanian HS. Started a Exploring_Tas_HS.Rmd file.

* 26-05-2020 0.0.0.8100 Modified getgrad4 to match current Tasmanian practice (divide each cpue chunk by the cpue in the first year of the chunk; basically I forgot but will document it thoroughly), changed getscore1 to getscore as it now works for both grad1 and grad4, and will work for gradX. Now need to add a targCE and a way to combine them in a weighted fashion and the HS and HCR are ready to go.

* 10-06-2020 0.0.0.8000 Altered all cases of region and reg into zone. Also altered all cases of SMU and smu to SAU and sau. So the data files have been altered accordingly.

* 10/08/2020 0.0.0.7950 Added to some of the vignettes and made minor improvements to some functions.

* 01/09/2020 0.0.0.7900 Added more to the Running the MSE vignette, which included the addition of a getzoneprod to summarize the productivity of the whole zone.

* 02/09/2020 0.0.0.7850 Modified generate_results.r bare-bones.r, devel.r, and README.Rmd to reflect changes I have made both in makehtml and in aMSE.

* 04/09/2020 0.0.0.7800 Have now stabilized makehtml and the places in aMSE where it is used (generate_results.R). Now it behaves as it should and it will be simple to add further plots and tables to the characterization.

* 13/09/2020 0.0.0.7750 Have removed 'diagrams' from the requirements of running the MSE and started adding functions and data-sets to faciliate the full blown operating model conditioning.

* 22/09/2020 0.0.0.7700 Have started adding files to begin conditioning the operating model on teh west coast data. This is designed to bring out the issues and problems raised when including such things insode the package. So far this has entailed changes to the readfiles for the control, zone, and data files. This probably means I will need now to fix the bare-bones example in teh readme.md file on github.

* 29/09/2020 0.0.0.7650 Have added functions to assist with the conditioning of the OM. Need to update the vignette. All examples should work, but some of them too need modification. 

* 01/10/2020 0.0.0.7600 Am working on conditioning the model and am now using the observed CPUE in an attempt to match it to the predicted CPUE from teh OM after initial depletion levels are set and the historical catches area applied. I have added a new candidatefunctions.R file to the data-raw folder and will expect to make changes pretty much every day, some small some large.

* 09/10/2020 0.0.0.7550 Have begun making changes to the characterization of zones, both unfished and depleted. Making functions internally consistent with the requirements of getting an MSZE up and running. Added more characterization functions.

* 11/10/2020 0.0.0.7500 Have transferred work to my travel laptop. Today have started to implement setting up the harvest strategy to be used in the projections after having generated a generic MSE operating model. This entailed modifying or adding to the zone definition .csv file as well as the readzonefile function and others to ensure internal consistency. I may eventually scrap the separation of the control file and the zone files and have it all in the control file, which currently has little remaining.

* 14/10/2020 0.0.0.7450 Have started to iron out inconsistencies between different ways of running the MSE. It needs an initial LML in case there is no conditioning and no projection planned. The initial LML will influecen estimates of exlpoitable biomass, MSY, and other aspects of production. So have modified the template generating files and the read data file functions.

* 26/10/2020 0.0.0.7400 Have altered the input control file to contain the zone wide information that was previously read in seperately. This entailed changing the ctrlfiletemplate and the readctrlzone functions. I have also been modifying the functions that assist with the automatic documentation of functions and the aMSE functions they reference. Not surprisingly, the makezone and oneyear related functions are referencing the most aMSE functions, while many others reference none.

* 01/11/2020 0.0.0.7350 Now have the route through the generation of the pre-projection zone clarified and consistent. Using the generic MSE, once depleted to a pre-determined level the use of the function addrepvar literally adds variation to the beginning of each individual replicate prior to application of the chosen HCR.

* 08/11/2020 0.0.0.7300 Implemented a constant Catch (constant TAC) harvest strategy so as to test the run time of a typical run. 1000 replicates in about 73 seconds.

* 09/11/2020 0.0.0.7250 Included new test data zone1, ctrl, and constants, with which to begin testing the MSE framework

* 22/11/2020 0.0.0.7000 Now includes an operational mcdahcr and applymcda functions. Currently 1000 replicates should take about 2.5 minutes as I have been optimizing again. There is much still to wrap up inside functions and plot convenient outputs, but it is already interesting.

* 24/11/2020 0.0.0.6950 Started a the_data_files.Rmds vignette, that attempts to describe the structure and contents of the data files (being the control.csv and populatiom.csv). Sone other minor adjustments, such as implementing the use of the randomseed.

* 25/11/2020 0.0.0.6900 Added first zone wide summary of outputs, still need to include length-composition of catch

* 03/12/2020 0.0.0.6800 Added separate randomseeds for population definition and projections. Added wtedmean and worked on vignettes. Modified zone.RData.

* 05/12/2020 0.0.0.6700 Exploring the initiation of projections.

* 06/12/2020 0.0.0.6600 Started developing a method to calibrate the MCDA using fishery data, without needing to condition the model completely. Fixed up examples so they now work with available data-sets.

* 10/12/2020 0.0.0.6500 Generalized the selection of the applicable HS but now need to define applyHS before running doprojection.

* 24/12/2020 0.0.0.6400 Strangely readctrlfile seemed to have reverted to readctrlzone!!? Fixed that and included applyHS as an explicit pointer to the harvest strategy function in doprojection.

* 05/01/21 0.0.0.6300 Added the preferred option of using a time-series of historical catches to deplete the unfished zone. Also added 'qmult' so as to adjust the cpue output, but may revert to just using the 'MaxCEpar'

* 12/01/2021 0.0.0.6200 Increasing generality now have calibrateMCDA and dohistoricC, plus extra HS related functions. More to come.

* 16/02/2021 0.0.0.6000 Now all components are currently operational.
