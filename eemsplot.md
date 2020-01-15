# Estimating effective migration surfaces [(EEMS)](https://github.com/dipetkov/eems)

Here is the code used to run eems_sats in command-line and plot the results in R. Detailed steps for running this program can be accessed [here](https://github.com/dipetkov/eems/tree/master/runeems_sats).

The program runeems_sats implements the EEMS method for analyzing spatial population structure. This version uses raw microsatellite data.

### Input data format (microsatellites)
`runeems_sats` requires three data input files that have the same file name but different extension. For this dataset, the file name is `microsat` + extension.

#### microsat.sites
`microsat.sites` is the matrix of allele copies; missing alleles are specified by any negative number. For n individuals at L loci, the sites matrix is n-by-L for L haploid markers and n-by-2L for L diploid markers. The matrix of allele copies was generated from `Genepop.txt`.

#### microsat.coord
`microsat.coord` are the sample coordinates, two coordinates per sample, one sample per line. The sampling locations should be given in the same order as the rows and columns of the dissimilarity matrix. Coordinates were taken from `Calcar_TESS2 - randGeo.txt`.

#### microsat.outer
`microsat.outer` are the habitat coordinates, as a sequence of vertices that form a closed polygon. The habitat vertices should be listed counterclockwise and the first vertex should also be the last vertex, so that the outline is a closed ring. Otherwise, EEMS attempts to "correct" the polygon and prints a warning message. Habitat coordinates were made by A. Zyck based on site coordinates.

All files were created in Excel, copied into terminal using `nano`, and saved as the file name with the appropriate extension.

#### Create separate directory for EEMS files
```bash
mkdir microsat_JP

cd microsat_JP
```

All input files should be saved in this directory.

#### Create a parameter file to set EEMS program arguments
There are a number of program parameters that can be set by the user. Here is the parameter file I created, named `params-microsat-chain1.ini`. Read the [instruction manual](https://github.com/dipetkov/eems/blob/master/Documentation/EEMS-doc.pdf) for descriptions of each parameter.
```R
datapath = ./microsat_JP/microsat
mcmcpath = ./microsat_JP/microsat-D200-chain1
nIndiv = 600
nSites = 6
nDemes = 200
diploid = true
numMCMCIter = 5000000
numBurnIter = 1000000
numThinIter = 9999
qEffctProposalS2 = 0.03
mEffctProposalS2 = 3.5
mrateMuProposalS2 = 0.1
mSeedsProposalS2 = 0.1
qSeedsProposalS2 = 0.1
```
- The first 9 parameter have to be specified (they have no default parameters).
- The last 5 parameters are *variances for the proposal distribution*. (I had to play around with these values, and changed most from the default values)
- As `runeems_sats` is running, it outputs information about the frequency of accepting proposals of different types.
  - The goal is to choose the parameters so that proposals are accepted about 20% − 30% of the time.
  - In practice, it seems to be sufficient to choose the variances so that proposals are not accepted too rarely (less than 10% of the time) or
too often (more than 40% of the time).

I created two additional parameter files, one with 300 demes `params-microsat-chain2.ini`, the other with 600 demes `params-microsat-chain3.ini`. All other parameters were kept the same.

### Run the EEMS program
```bash
runeems_sats --params params-microsat-chain1.ini
```
Repeat with the other two parameter files.

*NOTE: As you increase the number of demes, the program takes longer to run.*

### Plotting
In RStudio

#### Install rEEMSplots package
This R package provides the function `eems.plots` to visualize the estimated effective migration and diversity surfaces. It is not on CRAN, so install it from source instead. In the R console:
```R
## Check that the current directory contains the rEEMSplots source directory (from GitHub)
if (file.exists("./rEEMSplots/")) {
  install.packages("rEEMSplots", repos=NULL, type="source")
} else {
  stop("Move to the directory that contains the rEEMSplots source to install the package.")
}
```
```R
library(rEEMSplots)
```` 	
#### Plotting EEMS after running runeems_snps

```R
path = "/home/azyck/Fiddler_Crab/ddocent_env/microsat_JP/"
dirs = c(paste0(path,"microsat-D200-chain1"), paste0(path,"microsat-D300-chain1"), paste0(path,"microsat-D600-chain1"))
```
```R
eems.plots(mcmcpath = dirs, plotpath = "/home/azyck/Fiddler_Crab/ddocent_env/microsat_JP/microsat-All-plots",
           longlat = FALSE,add.grid=F,add.outline = T,add.demes = T, projection.in = "+proj=longlat +datum=WGS84", projection.out = "+proj=merc +datum=WGS84",
           add.map = T,add.abline = T, add.r.squared = T)
```

This will produce several different .png outputs that are saved to the directory specified in the `plotpath`.  
