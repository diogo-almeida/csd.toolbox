# To-do list for the csd.toolbox package

## Importing data

For the time being, there is no dedicated import data function. The only solution for the time being is to export the EEG/ERP data as ASCII text, and then use the regular R functions to read the data in. It is unlikely that I will invest time in creating import functions to directly read binary data from proprietary formats, but I should probably collect ASCII-exported data from a couple of different EEG/ERP systems so I can at least create some useful import functions for people who might be discouraged to use the package otherwise.

## Plotting functions

There are no dedicated plotting functions for the EEG and CSD data. I need to add at the very least:

* Single electrode plotting function
* Multiple electrode plotting function

### Topographic maps

* topographic maps based on any arbitrary eeg montage
* port other code from original MATLAB toolbox that might help dealing with the above
* email Dan Wu for her code with Dezhong Yao (paper: Wu & Yao (2007) The azimuth projection for the display of 3-D EEG data). I had this at some point, but I lost it somehow (though look for backups from my old computer).


