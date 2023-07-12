## aMSE Recent Activity.

See the ReadMe for the latest few events

aMSE v0.0.19 (Release date: 2022-10-26)
==============

Changes:

* added more functions to facilitate comparing scenarios.


aMSE v0.0.18 (Release date: 2022-10-18)
==============

Changes:

* In the projections replicates across years has now been changed to years across replicates to allow for inclusion of FIS performance measures. Also, further comparison functions included.

* 2022-10-26 aMSE 0.0.19 
* 2022-10-18 aMSE 0.0.18 In the projections replicates across years has now been changed to years across replicates to allow for inclusion of FIS performance measures. Also, further comparison functions included.

* 2022-09-26 aMSE 0.0.17 Added finalcondyeardepletion to aid in comparisons.

* 2022-08-24 aMSE 0.0.16 Started to add depensation and other refinements.

* 2022-08-19 aMSE 0.0.15 Added do_comparison as a wrapper function for comparing scenarios.

* 2022-08-05 aMSE 0.0.14 Added more functions to assist with comparing scenarios. 

* 2022-08-05 aMSE 0.0.12 Added more functions to assist with comparing scenarios. Still sorting out references to other packages.

* 2022-08-01 aMSE 0.0.11 Added a bunch of functions to assist with comparing scenarios. No EWN.

* 2022-07-10 aMSE 0.0.10 Switched dependency from hutils to codeutils. re-added zone1 for S21 to the data sets. No EWN.

* 2022-07-10 aMSE 0.0.9 Added data sets to allow for examples and documentation to illustrate the conditioning and running of the MSE. Added adjustavrec as a function, nade a wrapper rnormz and rlnormz to enable zero variation while still using up a random number so the sequence was maintained. Ready for production work. No EWN.  

* 2022-06-29 aMSE 0.0.8 Corrected a large mistake where the estimated SAU harvest rates were being biased high because the SAU's end of year exploitable biomass was being used instead of the mid-year exploitable biomass. 

* 2022-06-22 aMSE 0.0.7 Now have rewritecontrolfile and rewritesaudatafile and the option to include a table of parameters from sizemod into the inputs of aMSE. Have begun development of functions to aid comparing the output from different scenarios, including tables and plots. No EWN.

* 2022-06-08 aMSE 0.0.6 Have revised the ctrltemplate and datatemplate generators so that they reflect a lambda = 0.75 and correctly use the hyperstability. No EWN

* 2022-06-06 aMSE 0.0.5. Have reviewed the dynamics in detail so that the outputs from the aMSE for each SAU more closely match the outputs from sizemod. The original disjunction was partly due to a faulty implementation of the inclusion of the option of hyper-stability, but that has now been amended and the operatnig model is now behaving as it should. Breaking the dynamics within each SAU into multiple sub-populations with variation also affects the equivalence between aMSE and sizemod butmostly influences the productivity rather than the recruitment dynamics and the trajectories mainly follow those from sizemod but with the cpue and biomasses being projected incorrectly. It appears (still working on a why) that by changing the productivity and sub-dividing that productivity among sub-populations, the overall recruitment levels required to mimic the dynamics in all cases appears to be less in aMSE than in sizemod. Otheriwse nothing now needs to change. One uses parameters obtained from sizemod, especially the recruitment deviates, and then apply the functions that search for more optimal values for AvRec, the unfished recruitment levles in each SAU. The fits obtained should be close to those found in sizemod. In 0.0.5 there remain issues relating to the internal data sets, which include one called 'zone', which is not designed to work with the option of hyper-stability (lambda != 1.0). Otherwise, No EWN.

* 2022-05-31 aMSE 0.0.4 Added functions that assist in conditioning the operating model by optimizing the fit to CPUE and size-composition data after opimizing the AvRec values for each SAU. No EWN.

* 2021-11-29 aMSE 0.0.1 Removed datadir. Started implementing hyper-stability while keeping the interface the same, so that the results from sizemod can be used.

* 2021-11-23 aMSE 0.0.0.1 Added functions to simplify the optimization of AvRec relative to the CPUE data after initial parameter estimation for each SAU using the sizemod package. Also added a comparison of size-composition data to the do_condition function. Next step will be to include the optimization of the recdevs using both CPUE and size-comp data. 

* 2021-11-21 aMSE 0.0.0.020 So as to use any observed size-composition data in the conditioning, a number of functions have been added to convert the predicted numbers-at-size in the catch by population into NAS x SAU. These can now be plotted against the observed NAS for each SAU and this can be either to the console or to a new tab in the webpage called 'predictedcatchN'. This work derived from the progress made with the sizemod package used to fit size-based assessment models to individual SAU. Because the MSE sub-divides each SAU into a number of populations each with somewhat different MaxDL, L95, etc, even after fitting the models, when using that data it is then necessary to refit the AvRec values, which, because of the sub-division into populations each having different biology can lead to the MSE AvRec being either smaller or larger than the estimated R0 for a SAU. Similar arguments apply to teh recruitment residuals but generally the pattern provides a reasonable fit.

* 2021-10-17 aMSE 0.0.0.010 As cleanslate, within setuphtml() has been deprecated in the package makehtml, I have removed reference to it throughout aMSE. Now, to clean a directory of aMSE files use cleanrundir.

* 2021-09-22 aMSE 0.0.0.025 Modifications to do_condition and other functions to facilitate the ease of changing values within an array of control files and saudata files.

* 2021-09-17 aMSE 0.0.0.050 Finalized details of the optional initial depletion. Amended the logistic function so it now works when using the knifeedge option!
 
* 2021-09-13 aMSE 0.0.0.075 Implemented the option of having an initial depletion as well as application of historical catches.

* 2021-09-10 aMSE 0.0.0.100 Many small, but significant, changes. Modified getdata which interacts with the ###HS.R file and added the
zoneDP$NAS data input. Added a summary recdev plot to condition tab.

* 2021-09-03 aMSE 0.0.0.200 Lots of tidying of help pages, added numbersatsizeSAU. Added a fixed TasmanianHS.R file to the data-raw directory.

* 2021-08-23 aMSE 0.0.0.300 Using text as SAU names in the control file now labels all plots and those labels are not included in the tables. glb$saunames now contains the text labels. Other changes include generalizing the number of plots to match number of SAU in all cases (I hope), There is also, now, a ConstantCatchHS.R file in the data-raw sub-directory that allows for a constant catch to be applied to each SAU. This demonstrates that the issue when using the TasmanianHS.R file is not the aMSE package but how the Tasmanian HS functions work with the available data.

* 2021-08-19 aMSE 0.0.0.400 Can now have text as SAU names in the Control file, this will use those names to label all plots (still working on tbles and extra polishing).

* 2021-08-13   aMSE 0.0.0.500 removed NumNe (midyear numbers-at-size) and introduced cutcatchN as an argument to do_MSE, which removes all size class data from teh catchN array with size-classes less than cutcatchN, whose default = 56 (2 - 110mm removed). This saves a good deal of space.

* 2021-08-12 aMSE 0.0.0.500 Modified do_MSE and added plothsstats, streamlined the conditioning. No EWN.

* 2021-08-02 aMSE 0.0.0.700 Revised the 'ctrlfiletemplate' and 'datafiletemplate' functions to reflect the new usage. It should now be possible to run an example scenario related to M15h75.

* 2021_07-28 aMSE 0.0.0.800 Added functions (changecolumn, getrecdevcolumn, gettasdevssq) to assist in automated operating model conditioning on the fishery. Also modified the biology_plots and oneyear and oneyearcat to assist with speeding the processes. Lots of additional minor changes but some important ones to do_MSE (read the help, ?do_MSE).

* 2021-07-20 aMSE 0.0.0.900 Added do_condition and compareCPUE. Both used to speed the conditioning of the operating model, though a number of developments are still under development. Currently can automatically search for the AvRec value that will optimize the sum-of-squared residual fit between the
observed CPUE in the historical period, and those predicted by the conditioned model. The undeveloped bit related to the ad hoc recruitment deviates.

* 2021_07-12 aMSE 0.0.0.1000 Announcing some rather large changes. Early on in development I made a strategic decision which has proved to make the documentation of the code and implementation of extraction of results overly complex. It is now clear that it was a strategic mistake to separate the historical conditioning from the future projections. Combining them, for purposes of applying a Harvest control rule, or tabulating or plotting results, was adding complexity to every scenario, new plot and table. For example, previously, plots where there were traces of the historical data added to the projections were greatly flawed by having the final year of historical data repeated in the projection data. Obviously, that plotting error could have been solved by selecting which years to plot, instead I decided to simplify everything (for both me and other maintainers) by combining the conditioned dynamic component of the zone and the projection years of the dynamics. Now, in the 'zoneDP' object (see 'out$zoneDP' after a run), for each replicate, we have a continuous set of years from the first year of historical catches out to the last year of the projections. This simplifies everything from the documentation of the code base, the application of the harvest strategy, the plotting of results, basically everything I can think of. The next steps are to expand the results sections and re-write the documentation. The current Tasmania example has now been expanded to 56 populations among the 8 SAU. The fit of the operating model the historical CPUE is now rather improved because now it is possible to include recruitment deviates in the historical conditioning. Currently in the included example it is very ad hoc. One can use ad hoc recruitment deviates or ones derived from some form of localized model fitting. The ad hoc ones currently in use were mainly used to test the new utility functions provided to assist with such things. 

* 2021_07-08 aMSE 0.0.0.1500 Another big change. Functions added to aid conditioning, and aspects of the projections streamlined.  No errors, warnings, or notes

* 2021-07-04 aMSE 0.0.0.1900 More big changes. I have removed the numbers-at-size large matrices from the zoneDP object into a NAS object. zoneDP is now only 12Mb instead of 450Mb and so can now be stored conveniently for each run. In addition, fixed recruitment deviates have now been implemented (although the default is to not use them - all are set to -1). But, expreimentally, in Tasmania, setting the deviate in 1991 = 2.0 meaning that the base level of recruitment off the recruitment curve is doubled, improves the relationship between predicted cpue and observed cpue so this holds great promise for improved conditioning.

* 2021-06-28 aMSE 0.0.0.2000 The AMSE can now conduct what are effectively retrospective analyses. This is affected by changing the 'CATCHES, 47, if>1 this is the number of historical catch years' line in the control.csv file, which can now be done in a loop using the new 'changevar' function. This has necessitated also changing the 'CEYRS, 20, The first year of catches for which there are CPUE records, ie 1992 is 20th year of catches', which used to take a direct value of how any years of CPUE data there are. We could have just changed two numbers but along that path errors lay.

* 2021-06-24 aMSE 0.0.0.2100 Big jump in number as some larger changes. I have removed a bunch of deprecated functions (which are now in 'deprecated.R' in the data-raw directory). I have tidied many other functions, and have modified listall files to aid in the auto-documentation of the package.

* 2021-05-27 aMSE 0.0.0.2400 Have encapsulated the bulk of running the MSE for a particular scenario into the function 'do_MSE'. Check out its help page but also see it in action in the example within the readme file or in the new version of first_use_saudata.R in ~Dropbox\National abalone MSE\aMSE_files\scenarios\HS652510. I have also updated the aMSE_0.0.0.2400.tar.gz file. Any issues get in touch. 

* 2021-05-26 0.0.0.2500 Some relatively large changes this time. I have pulled most of the code used to run the MSE into a separate source file 'MSE_source.R'. You can find all required files in  \Dropbox\National abalone MSE\aMSE_files\scenarios. I have moved TasmanianHS.R into tasdata as it is currently common to different uses. I will endeavour to generalize the contents of MSE_source.R to turn it into a function that can be included in the package. Then any changes made there will flow naturally to other instances of each scenario file. Also, now there is a minimu requirement of 2 SAU but each SAU can have a minimum of 1 population.

* 2021-05-11 0.0.0.2650 Modifications to allow for greater internal consistency and to prepare for more generalizations of the code.

* 2021-04-25 0.0.0.2700 Some large changes relating to generalization. All plots now have correctly labelled year axes. The doprojection function requires functions to sampleCE, sampleFIS and sampleNaS data. It also needs a getdata function which uses the three sample functions to make a data object of the hcrfun, then we need a calcpopC function that uses the output of the hcrfun to calculate the actual catches to be taken from each population. This all needs better documentation but that will take time to catch up with all these changes.

* 2021-04-20 0.0.0.2900 Continuing process of modifying the core functions to allow for the numbers-at-size prior to fishing and the midyear-exploitB (midyexpB) to be included in teh zoneD object. doprojection still requires fursther modification to account for applying to chosen HCR and HS to the completed conditioned zoneDD to provide the first year of expected catches in the projections (which are to be defined by teh conditoned data prior to projections beginning).

* 2021-04-19 0.0.0.3000 Started development required to generalize the use of the doprojection function (it was previously called doTASprojection). This generalization requires the generation of separate functions to process the data for cpue, any FIS, and the numbers-at-size. Numerous diganostic plots have been added, also providing a template for adding more  (see News.md for history of development).

* 2021-04-14 0.0.0.3200 Added NumNe to zoneDP output ready for its use in estimating FIS results by population and other numbers-at-size related indices. NumNe are the numbers-at-size following growth and the application of half of natural mortality but before any fishing moretality has occurred. 

* 11/04/2021 0.0.0.3300 modified the projection function and HCR function to generate and use the aspirational catches and include the TAC. Added a function to plot the zonescale TAC.

* 30/03/2021 0.0.0.3500 Working on working with length-composition data both within the projections and during the conditioning. Added plotNt, plotCnt, prepareDDNt, and getLFdata.

* 24/03/2021 0.0.0.4000 Fixed the summary of populations to SAU and the plotting of those results.

* 22/03/2021 0.0.0.4500 Added the option of using an SAU based data file for reading in the biological properties of the zone; this simplifies the conditioning, especially of the recruitment levels by population.

* 17/03/2021 0.0.0.5000 Well behind in this documentation. The MSE now operates successfully, with the README page on GitHub providing a worked example. The latest versions of aMSE, makehtml, and rutilsMH, as tar.gz source files along with a control2.csv and westzone.csv data files, and a first-run.R file are to be found in National abalone MSE/aMSE_files. Currently, the development of additional functions to generate further summary plots and tables are still being worked on. And work towards correctly conditioning the Operating Model onto the Tasmanian west coast is underway.

* 26/02/2021 0.0.0.5600 Produced SAU summary using alltosau. 

* 19/02/2021 0.0.0.5750 Starting to improve the HCR to more closely match current usage

* 16/02/2021 0.0.0.6000 Now all components are currently operational.

* 12/01/2021 0.0.0.6200 Increasing generality now have calibrateMCDA and dohistoricC, plus extra HS related functions. More to come.

* 05/01/21 0.0.0.6300 Added the preferred option of using a time-series of historical catches to deplete the unfished zone. Also added 'qmult' so as to adjust the cpue output, but may revert to just using the 'MaxCEpar'

* 24/12/2020 0.0.0.6400 Strangely readctrlfile seemed to have reverted to readctrlzone!!? Fixed that and included applyHS as an explicit pointer to the harvest strategy function in doprojection.

* 10/12/2020 0.0.0.6500 Generalized the selection of the applicable HS but now need to define applyHS before running doprojection.

* 06/12/2020 0.0.0.6600 Started developing a method to calibrate the MCDA using fishery data, without needing to condition the model completely. Fixed up examples so they now work with available data-sets.

* 05/12/2020 0.0.0.6700 Exploring the initiation of projections.

* 03/12/2020 0.0.0.6800 Added separate randomseeds for population definition and projections. Added wtedmean and worked on vignettes. Modified zone.RData.

* 25/11/2020 0.0.0.6900 Added first zone wide summary of outputs, still need to include length-composition of catch

