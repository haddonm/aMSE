## aMSE Recent Activity.

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

