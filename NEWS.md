## aMSE Recent Activity.

See the ReadMe for the latest two events

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

