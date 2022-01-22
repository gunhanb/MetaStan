
In this submission, I have re-written some code, extend some functionalities, and add unit tests.

## Test environments
* win-builder (devel and release)



## Re-submission

This is a resubmission.

In this submission, I fixed the error in the "convert_data_arm" function.

## Test environments
* win-builder (devel and release)
* Ubuntu (devel and release)




## New submission

This is a new submission, since I missed the deadline.

In this submission, I have replaced "donttest" with "dontrun".

## Test environments
* win-builder (devel and release)
* local OS X install, R 3.6.3



## This is a new submission

In this submission, I have corrected the Stan parameterization used in the "meta_stan" function, and updated the reference paper.

## Test environments
* win-builder (devel and release)
* ubuntu 18.04, R 3.6.1
* local OS X install, R 3.6.1

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs. 


## This is a new submission

In this submission, I have
* extended the package by including functions to fit model-based meta-analysis models
* added two functions: create_MBMA_dat and MBMA_stan
* added an example dataset for model-based meta-anylsis

## Test environments
* win-builder (devel and release)
* ubuntu 18.04, R 3.6.0
* local OS X install, R 3.6.0

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs. 


## This is a resubmission

In this submission, I have
* updated meta_stan function by adding "data" argument to make it more user-friendly.


## This is a new submission, since I missed the deadline.

In this submission, I have
* updated Makevars and Makevars.win files following Prof. Ripley comments.

## This is a resubmission

In this submission, I have
* added a reference describing methods of the package to DESCRIPTION
* added more examples to Rd-files
* replaced \dontrun with \donttest
* fixed the typo in man/meta_stan.Rd

  ## R CMD check results
  There is 1 NOTE:
  
  Possibly mis-spelled words in DESCRIPTION:
    Friede (7:35)
    Guenhan (7:14)
    Roever (7:23)
  
  These are author names.

## Test environments
* ubuntu 18.04, R 3.5.1
* local OS X install, R 3.5.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:

* This is a new submission.

* checking for GNU extensions in Makefiles ... NOTE

  GNU make is a SystemRequirements.

GNU make is added as SystemRequirements in Description 
of the package
