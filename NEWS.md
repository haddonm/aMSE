## aMSE Recent Activity

* Added a `NEWS.md` file to track/log the development of the package.

* Currently can generate a full region, characterize its productivity, and if run with no catch or variability maintains an equilibrium, even in the presence of larval drift.

* 16-04-2020 worked on setuphtml and addfilename to simplify the use of make-html when adding a plot or table to the results for a simulation.

* 17-04-2020 0.0.0.8920: Cleaned up the make_html code to produce valid HTML5 asnd CSS3, separated the CSS code into a separate file that is linked into the different html files. Cleaned up the imports and fixed a bunch of the examples, so that there are now no errors, no warnings, and no notes. Although the examples for the read[ctrl][data][region1] functions all currently have a 'dontrun' status. An issue to fix. 

* 01-05-2020 0.0.0.8600: Continued with the development of Running_aMSE.Rnd. This has the great advantage of forcing me to clarify the text and the relationships between the R functions and the tasks they do. It even suggests to me ways to improve on the names of the functions. Naming things continues to be one of the harder things (if I want to keep these things meaningful, nothing wrong with a1, a2, x1, x2, etc). Have discovered how to include some CSS code to customize the vignette format. But, more importantly, I have improved the text, changed addfilename to logfilename and automated recording the file type (which is used by make_html to process the file correctly on the local webpage).

* 10-05-2020: 0.0.0.8500 Continued the development of Running_aMSE.Rmd. Clarified the intent of this document. Cleaned up the makehtml code and added the dodepletion section. Now need to implement at least one HCR funciton and associated file.

* 11-05-2020 0.0.0.8450 Put the results generation into internal functions so their production can be automated as new ones are developed. Ensured all internal functions were operating consistently with regard the variolus name changes I have made. e.g. fixed the emergence calculations, which was really messing with the results! Set up a better worked example in the readme in the github version so it now produces plots and tables and can generate the local website shoudl the user wish.

* 12-05-2020 0.0.0.8400 Simplified the intended directory structure to a single results directory _resdir_. Cleaned out old, now unused code, and tidied the help files.

* 15-05-2020 0.0.0.8350 Started including estimates relating to SMUs as well as populations and whole of region. This entailed making changes to the design of regionD, to ensure internal consistency and ease of use. Began adding details of the fleet dynamics to Running-aMSE.Rmd.
