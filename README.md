# Quadrophenia
## A Tetraploid Dosage Calling R Package

__Note: This was edited to remove confidential information. It may have become unstable in the process.__ 

After you have installed R, and/or RStudio, Quadrophenia can be installed like any R Package. The one caveat being that Quadrophenia relies on a Perl interpreter for some of its functionality. So, you will also need to install Perl on your machine. 

### Installation
The first step is to ensure that you have Perl installed on your computer. This can be downloaded at [https://www.perl.org/get.html](https://www.perl.org/get.html).

Once Perl is installed, install RStudio from the Cran homepage at [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/). 

Then install the Devtools package in RStudio by entering the following into the RStudio command line.

    install.packages("devtools")

Load this into the local environment with 

    library(devtools)

Since this Github repository is not public, you will not be able to install this package with the Devtools R library. However, it is simple enough to download this package by clicking ___Clone or download___ to the right, or use this [link](https://github.com/dsherma7/Quadrophenia/archive/master.zip). 

Once downloaded, you can then use the devtools install() function as

    install([path to Quadrophenia]/Quadrophenia/)

This will also install the other dependencies including plyr, gdata, and xlsx. 


