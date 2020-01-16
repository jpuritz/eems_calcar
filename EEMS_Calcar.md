EEMS\_Calcar
================
Jon Puritz
1/15/2020

Estimating effective migration surfaces [(EEMS)](https://github.com/dipetkov/eems)
==================================================================================

Here is the code used to run eems\_sats in command-line and plot the results in R. Detailed steps for running this program can be accessed [here](https://github.com/dipetkov/eems/tree/master/runeems_sats).

The program runeems\_sats implements the EEMS method for analyzing spatial population structure. This version uses raw microsatellite data.

### Input data format (microsatellites)

`runeems_sats` requires three data input files that have the same file name but different extension. For this dataset, the file name is `microsat` + extension.

#### microsat.sites

`microsat.sites` is the matrix of allele copies; missing alleles are specified by any negative number. For n individuals at L loci, the sites matrix is n-by-L for L haploid markers and n-by-2L for L diploid markers. The matrix of allele copies was generated from `Genepop.txt`.

#### microsat.coord

`microsat.coord` are the sample coordinates, two coordinates per sample, one sample per line. The sampling locations should be given in the same order as the rows and columns of the dissimilarity matrix. Coordinates were taken from `Calcar_TESS2 - randGeo.txt`.

#### microsat.outer

`microsat.outer` are the habitat coordinates, as a sequence of vertices that form a closed polygon. The habitat vertices should be listed counterclockwise and the first vertex should also be the last vertex, so that the outline is a closed ring. Otherwise, EEMS attempts to "correct" the polygon and prints a warning message. Habitat coordinates were made by A. Zyck based on site coordinates.

### All files are included in the repository

Run Parameters
--------------

``` bash
cat params-microsat-D48-chain1.ini
```

    ## datapath = /home/jpuritz/eems_calcar/microsat
    ## mcmcpath = /home/jpuritz/eems_calcar/microsat-D48-chain1
    ## nIndiv = 600
    ## nSites = 6
    ## nDemes = 48
    ## diploid = true
    ## numMCMCIter = 1000000
    ## numBurnIter = 100000
    ## numThinIter = 999
    ## qEffctProposalS2 = 0.04
    ## mEffctProposalS2 = 3.75  
    ## mrateMuProposalS2 = 0.15 
    ## mSeedsProposalS2 = 0.15  
    ## qSeedsProposalS2 = 0.15

``` bash
cp params-microsat-D48-chain1.ini params-microsat-D48-chain2.ini
cp params-microsat-D48-chain1.ini params-microsat-D48-chain3.ini
sed -i 's/chain1/chain2/g' params-microsat-D48-chain2.ini
sed -i 's/chain1/chain3/g' params-microsat-D48-chain3.ini
```

``` bash
for i in {96,192,384,768}
do
  for j in {1,2,3}
  do
  cp params-microsat-D48-chain1.ini params-microsat-D$i-chain$j.ini
  sed -i "s/chain1/chain$j/g" params-microsat-D$i-chain$j.ini
  sed -i "s/D48/D$i/g" params-microsat-D$i-chain$j.ini
  sed -i "s/nDemes = 48/nDemes = $i/g" params-microsat-D$i-chain$j.ini
  done
done
```

Run eems in parallel 3 chains per deme level
--------------------------------------------

``` bash
ls params-microsat-D* | parallel runeems_sats --params {} >eems.out 2>&1
```

``` r
library("rEEMSplots")
library("ggplot2")
library("dplyr")
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library("rworldmap")
```

    ## Loading required package: sp

    ## ### Welcome to rworldmap ###

    ## For a short introduction type :   vignette('rworldmap')

``` r
library("rworldxtra")
```

``` r
path = "/home/jpuritz/eems_calcar/"
plotpath = "/home/jpuritz/eems_calcar/eems"
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"
```

``` r
dirs = c(paste0(path,"microsat-D48-chain1"),paste0(path,"microsat-D48-chain2"),paste0(path,"microsat-D48-chain3"),paste0(path,"microsat-D96-chain1"),paste0(path,"microsat-D96-chain2"),paste0(path,"microsat-D96-chain3"),paste0(path,"microsat-D192-chain1"),paste0(path,"microsat-D192-chain2"),paste0(path,"microsat-D192-chain3"),paste0(path,"microsat-D384-chain1"),paste0(path,"microsat-D384-chain2"),paste0(path,"microsat-D384-chain3"),paste0(path,"microsat-D768-chain1"),paste0(path,"microsat-D768-chain2"),paste0(path,"microsat-D768-chain3"))

coord <- read.table(paste0(path, "microsat.coord"))
coord.dec <- SpatialPoints(cbind(coord$V2, coord$V1), proj4string=CRS(projection_none))
coord.merc <- spTransform(coord.dec, CRS(projection_mercator))

eems.plots(mcmcpath = dirs, plotpath,longlat = FALSE,add.grid=F,add.outline = T,add.demes = F, projection.in = projection_none, projection.out = projection_mercator, add.map = T,add.abline = T, add.r.squared = T,add.title = FALSE,m.plot.xy = { points(coord.merc, col = "black", cex =1, pch=19) },q.plot.xy = { points(coord.merc, col = "black", cex =1, pch=19) })
```

    ## Input projection: +proj=longlat +datum=WGS84
    ## Output projection: +proj=merc +datum=WGS84

    ## Loading rgdal (required by projection.in)

    ## Loading rworldmap (required by add.map)

    ## Loading rworldxtra (required by add.map)

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Using 'euclidean' distance to assign interpolation points to Voronoi tiles.

    ## Processing the following EEMS output directories :

    ## /home/jpuritz/eems_calcar/microsat-D48-chain1/home/jpuritz/eems_calcar/microsat-D48-chain2/home/jpuritz/eems_calcar/microsat-D48-chain3/home/jpuritz/eems_calcar/microsat-D96-chain1/home/jpuritz/eems_calcar/microsat-D96-chain2/home/jpuritz/eems_calcar/microsat-D96-chain3/home/jpuritz/eems_calcar/microsat-D192-chain1/home/jpuritz/eems_calcar/microsat-D192-chain2/home/jpuritz/eems_calcar/microsat-D192-chain3/home/jpuritz/eems_calcar/microsat-D384-chain1/home/jpuritz/eems_calcar/microsat-D384-chain2/home/jpuritz/eems_calcar/microsat-D384-chain3/home/jpuritz/eems_calcar/microsat-D768-chain1/home/jpuritz/eems_calcar/microsat-D768-chain2/home/jpuritz/eems_calcar/microsat-D768-chain3

    ## Plotting effective migration surface (posterior mean of m rates)

    ## /home/jpuritz/eems_calcar/microsat-D48-chain1

    ## /home/jpuritz/eems_calcar/microsat-D48-chain2

    ## /home/jpuritz/eems_calcar/microsat-D48-chain3

    ## /home/jpuritz/eems_calcar/microsat-D96-chain1

    ## /home/jpuritz/eems_calcar/microsat-D96-chain2

    ## /home/jpuritz/eems_calcar/microsat-D96-chain3

    ## /home/jpuritz/eems_calcar/microsat-D192-chain1

    ## /home/jpuritz/eems_calcar/microsat-D192-chain2

    ## /home/jpuritz/eems_calcar/microsat-D192-chain3

    ## /home/jpuritz/eems_calcar/microsat-D384-chain1

    ## /home/jpuritz/eems_calcar/microsat-D384-chain2

    ## /home/jpuritz/eems_calcar/microsat-D384-chain3

    ## /home/jpuritz/eems_calcar/microsat-D768-chain1

    ## /home/jpuritz/eems_calcar/microsat-D768-chain2

    ## /home/jpuritz/eems_calcar/microsat-D768-chain3

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Plotting effective diversity surface (posterior mean of q rates)

    ## /home/jpuritz/eems_calcar/microsat-D48-chain1

    ## /home/jpuritz/eems_calcar/microsat-D48-chain2

    ## /home/jpuritz/eems_calcar/microsat-D48-chain3

    ## /home/jpuritz/eems_calcar/microsat-D96-chain1

    ## /home/jpuritz/eems_calcar/microsat-D96-chain2

    ## /home/jpuritz/eems_calcar/microsat-D96-chain3

    ## /home/jpuritz/eems_calcar/microsat-D192-chain1

    ## /home/jpuritz/eems_calcar/microsat-D192-chain2

    ## /home/jpuritz/eems_calcar/microsat-D192-chain3

    ## /home/jpuritz/eems_calcar/microsat-D384-chain1

    ## /home/jpuritz/eems_calcar/microsat-D384-chain2

    ## /home/jpuritz/eems_calcar/microsat-D384-chain3

    ## /home/jpuritz/eems_calcar/microsat-D768-chain1

    ## /home/jpuritz/eems_calcar/microsat-D768-chain2

    ## /home/jpuritz/eems_calcar/microsat-D768-chain3

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Plotting posterior probability trace

    ## /home/jpuritz/eems_calcar/microsat-D48-chain1

    ## /home/jpuritz/eems_calcar/microsat-D48-chain2

    ## /home/jpuritz/eems_calcar/microsat-D48-chain3

    ## /home/jpuritz/eems_calcar/microsat-D96-chain1

    ## /home/jpuritz/eems_calcar/microsat-D96-chain2

    ## /home/jpuritz/eems_calcar/microsat-D96-chain3

    ## /home/jpuritz/eems_calcar/microsat-D192-chain1

    ## /home/jpuritz/eems_calcar/microsat-D192-chain2

    ## /home/jpuritz/eems_calcar/microsat-D192-chain3

    ## /home/jpuritz/eems_calcar/microsat-D384-chain1

    ## /home/jpuritz/eems_calcar/microsat-D384-chain2

    ## /home/jpuritz/eems_calcar/microsat-D384-chain3

    ## /home/jpuritz/eems_calcar/microsat-D768-chain1

    ## /home/jpuritz/eems_calcar/microsat-D768-chain2

    ## /home/jpuritz/eems_calcar/microsat-D768-chain3

    ## Plotting average dissimilarities within and between demes

    ## /home/jpuritz/eems_calcar/microsat-D48-chain1

    ## /home/jpuritz/eems_calcar/microsat-D48-chain2

    ## /home/jpuritz/eems_calcar/microsat-D48-chain3

    ## /home/jpuritz/eems_calcar/microsat-D96-chain1

    ## Using 'euclidean' distance to assign interpolation points to Voronoi tiles.
    ## 
    ## 
    ## 
    ## Using 'euclidean' distance to assign interpolation points to Voronoi tiles.

    ## EEMS results for at least two different population grids

### EEMS runs per deme setting

#### 48

``` r
dirs = c(paste0(path,"microsat-D48-chain1"),paste0(path,"microsat-D48-chain2"),paste0(path,"microsat-D48-chain3"))
plotpath = "/home/jpuritz/eems_calcar/eems-D48"

eems.plots(mcmcpath = dirs, plotpath,longlat = FALSE, add.grid=F,add.outline = T,add.demes = F, projection.in = projection_none, projection.out = projection_mercator, add.map = T,add.abline = T, add.r.squared = T,add.title = FALSE,m.plot.xy = { points(coord.merc, col = "black", cex =1, pch=19) },q.plot.xy = { points(coord.merc, col = "black", cex =1, pch=19) })
```

    ## Input projection: +proj=longlat +datum=WGS84
    ## Output projection: +proj=merc +datum=WGS84

    ## Loading rgdal (required by projection.in)

    ## Loading rworldmap (required by add.map)

    ## Loading rworldxtra (required by add.map)

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Using 'euclidean' distance to assign interpolation points to Voronoi tiles.

    ## Processing the following EEMS output directories :

    ## /home/jpuritz/eems_calcar/microsat-D48-chain1/home/jpuritz/eems_calcar/microsat-D48-chain2/home/jpuritz/eems_calcar/microsat-D48-chain3

    ## Plotting effective migration surface (posterior mean of m rates)

    ## /home/jpuritz/eems_calcar/microsat-D48-chain1

    ## /home/jpuritz/eems_calcar/microsat-D48-chain2

    ## /home/jpuritz/eems_calcar/microsat-D48-chain3

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Plotting effective diversity surface (posterior mean of q rates)

    ## /home/jpuritz/eems_calcar/microsat-D48-chain1

    ## /home/jpuritz/eems_calcar/microsat-D48-chain2

    ## /home/jpuritz/eems_calcar/microsat-D48-chain3

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Plotting posterior probability trace

    ## /home/jpuritz/eems_calcar/microsat-D48-chain1

    ## /home/jpuritz/eems_calcar/microsat-D48-chain2

    ## /home/jpuritz/eems_calcar/microsat-D48-chain3

    ## Plotting average dissimilarities within and between demes

    ## /home/jpuritz/eems_calcar/microsat-D48-chain1

    ## /home/jpuritz/eems_calcar/microsat-D48-chain2

    ## /home/jpuritz/eems_calcar/microsat-D48-chain3

#### 96

``` r
dirs = c(paste0(path,"microsat-D96-chain1"),paste0(path,"microsat-D96-chain2"),paste0(path,"microsat-D96-chain3"))
plotpath = "/home/jpuritz/eems_calcar/eems-D96"

eems.plots(mcmcpath = dirs, plotpath,longlat = FALSE, add.grid=F,add.outline = T,add.demes = F, projection.in = projection_none, projection.out = projection_mercator, add.map = T,add.abline = T, add.r.squared = T,add.title = FALSE,m.plot.xy = { points(coord.merc, col = "black", cex =1, pch=19) },q.plot.xy = { points(coord.merc, col = "black", cex =1, pch=19) })
```

    ## Input projection: +proj=longlat +datum=WGS84
    ## Output projection: +proj=merc +datum=WGS84

    ## Loading rgdal (required by projection.in)

    ## Loading rworldmap (required by add.map)

    ## Loading rworldxtra (required by add.map)

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Using 'euclidean' distance to assign interpolation points to Voronoi tiles.

    ## Processing the following EEMS output directories :

    ## /home/jpuritz/eems_calcar/microsat-D96-chain1/home/jpuritz/eems_calcar/microsat-D96-chain2/home/jpuritz/eems_calcar/microsat-D96-chain3

    ## Plotting effective migration surface (posterior mean of m rates)

    ## /home/jpuritz/eems_calcar/microsat-D96-chain1

    ## /home/jpuritz/eems_calcar/microsat-D96-chain2

    ## /home/jpuritz/eems_calcar/microsat-D96-chain3

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Plotting effective diversity surface (posterior mean of q rates)

    ## /home/jpuritz/eems_calcar/microsat-D96-chain1

    ## /home/jpuritz/eems_calcar/microsat-D96-chain2

    ## /home/jpuritz/eems_calcar/microsat-D96-chain3

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Plotting posterior probability trace

    ## /home/jpuritz/eems_calcar/microsat-D96-chain1

    ## /home/jpuritz/eems_calcar/microsat-D96-chain2

    ## /home/jpuritz/eems_calcar/microsat-D96-chain3

    ## Plotting average dissimilarities within and between demes

    ## /home/jpuritz/eems_calcar/microsat-D96-chain1

    ## /home/jpuritz/eems_calcar/microsat-D96-chain2

    ## /home/jpuritz/eems_calcar/microsat-D96-chain3

#### 192

``` r
dirs = c(paste0(path,"microsat-D192-chain1"),paste0(path,"microsat-D192-chain2"),paste0(path,"microsat-D192-chain3"))
plotpath = "/home/jpuritz/eems_calcar/eems-D192"

eems.plots(mcmcpath = dirs, plotpath,longlat = FALSE, add.grid=F,add.outline = T,add.demes = F, projection.in = projection_none, projection.out = projection_mercator, add.map = T,add.abline = T, add.r.squared = T,add.title = FALSE,m.plot.xy = { points(coord.merc, col = "black", cex =1, pch=19) },q.plot.xy = { points(coord.merc, col = "black", cex =1, pch=19) })
```

    ## Input projection: +proj=longlat +datum=WGS84
    ## Output projection: +proj=merc +datum=WGS84

    ## Loading rgdal (required by projection.in)

    ## Loading rworldmap (required by add.map)

    ## Loading rworldxtra (required by add.map)

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Using 'euclidean' distance to assign interpolation points to Voronoi tiles.

    ## Processing the following EEMS output directories :

    ## /home/jpuritz/eems_calcar/microsat-D192-chain1/home/jpuritz/eems_calcar/microsat-D192-chain2/home/jpuritz/eems_calcar/microsat-D192-chain3

    ## Plotting effective migration surface (posterior mean of m rates)

    ## /home/jpuritz/eems_calcar/microsat-D192-chain1

    ## /home/jpuritz/eems_calcar/microsat-D192-chain2

    ## /home/jpuritz/eems_calcar/microsat-D192-chain3

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Plotting effective diversity surface (posterior mean of q rates)

    ## /home/jpuritz/eems_calcar/microsat-D192-chain1

    ## /home/jpuritz/eems_calcar/microsat-D192-chain2

    ## /home/jpuritz/eems_calcar/microsat-D192-chain3

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Plotting posterior probability trace

    ## /home/jpuritz/eems_calcar/microsat-D192-chain1

    ## /home/jpuritz/eems_calcar/microsat-D192-chain2

    ## /home/jpuritz/eems_calcar/microsat-D192-chain3

    ## Plotting average dissimilarities within and between demes

    ## /home/jpuritz/eems_calcar/microsat-D192-chain1

    ## /home/jpuritz/eems_calcar/microsat-D192-chain2

    ## /home/jpuritz/eems_calcar/microsat-D192-chain3

#### 384

``` r
dirs = c(paste0(path,"microsat-D384-chain1"),paste0(path,"microsat-D384-chain2"),paste0(path,"microsat-D384-chain3"))
plotpath = "/home/jpuritz/eems_calcar/eems-D384"

eems.plots(mcmcpath = dirs, plotpath,longlat = FALSE, add.grid=F,add.outline = T,add.demes = F, projection.in = projection_none, projection.out = projection_mercator, add.map = T,add.abline = T, add.r.squared = T,add.title = FALSE,m.plot.xy = { points(coord.merc, col = "black", cex =1, pch=19) },q.plot.xy = { points(coord.merc, col = "black", cex =1, pch=19) })
```

    ## Input projection: +proj=longlat +datum=WGS84
    ## Output projection: +proj=merc +datum=WGS84

    ## Loading rgdal (required by projection.in)

    ## Loading rworldmap (required by add.map)

    ## Loading rworldxtra (required by add.map)

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Using 'euclidean' distance to assign interpolation points to Voronoi tiles.

    ## Processing the following EEMS output directories :

    ## /home/jpuritz/eems_calcar/microsat-D384-chain1/home/jpuritz/eems_calcar/microsat-D384-chain2/home/jpuritz/eems_calcar/microsat-D384-chain3

    ## Plotting effective migration surface (posterior mean of m rates)

    ## /home/jpuritz/eems_calcar/microsat-D384-chain1

    ## /home/jpuritz/eems_calcar/microsat-D384-chain2

    ## /home/jpuritz/eems_calcar/microsat-D384-chain3

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Plotting effective diversity surface (posterior mean of q rates)

    ## /home/jpuritz/eems_calcar/microsat-D384-chain1

    ## /home/jpuritz/eems_calcar/microsat-D384-chain2

    ## /home/jpuritz/eems_calcar/microsat-D384-chain3

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Plotting posterior probability trace

    ## /home/jpuritz/eems_calcar/microsat-D384-chain1

    ## /home/jpuritz/eems_calcar/microsat-D384-chain2

    ## /home/jpuritz/eems_calcar/microsat-D384-chain3

    ## Plotting average dissimilarities within and between demes

    ## /home/jpuritz/eems_calcar/microsat-D384-chain1

    ## /home/jpuritz/eems_calcar/microsat-D384-chain2

    ## /home/jpuritz/eems_calcar/microsat-D384-chain3

#### 768

``` r
dirs = c(paste0(path,"microsat-D768-chain1"),paste0(path,"microsat-D768-chain2"),paste0(path,"microsat-D768-chain3"))
plotpath = "/home/jpuritz/eems_calcar/eems-D768"

eems.plots(mcmcpath = dirs, plotpath,longlat = FALSE, add.grid=F,add.outline = T,add.demes = F, projection.in = projection_none, projection.out = projection_mercator, add.map = T,add.abline = T, add.r.squared = T,add.title = FALSE,m.plot.xy = { points(coord.merc, col = "black", cex =1, pch=19) },q.plot.xy = { points(coord.merc, col = "black", cex =1, pch=19) })
```

    ## Input projection: +proj=longlat +datum=WGS84
    ## Output projection: +proj=merc +datum=WGS84

    ## Loading rgdal (required by projection.in)

    ## Loading rworldmap (required by add.map)

    ## Loading rworldxtra (required by add.map)

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Using 'euclidean' distance to assign interpolation points to Voronoi tiles.

    ## Processing the following EEMS output directories :

    ## /home/jpuritz/eems_calcar/microsat-D768-chain1/home/jpuritz/eems_calcar/microsat-D768-chain2/home/jpuritz/eems_calcar/microsat-D768-chain3

    ## Plotting effective migration surface (posterior mean of m rates)

    ## /home/jpuritz/eems_calcar/microsat-D768-chain1

    ## /home/jpuritz/eems_calcar/microsat-D768-chain2

    ## /home/jpuritz/eems_calcar/microsat-D768-chain3

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Plotting effective diversity surface (posterior mean of q rates)

    ## /home/jpuritz/eems_calcar/microsat-D768-chain1

    ## /home/jpuritz/eems_calcar/microsat-D768-chain2

    ## /home/jpuritz/eems_calcar/microsat-D768-chain3

    ## Using the default DarkOrange to Blue color scheme, with 'white' as the midpoint color.
    ## It combines two color schemes from the 'dichromat' package, which itself is based on
    ## a collection of color schemes for scientific data graphics:
    ##  Light A and Bartlein PJ (2004). The End of the Rainbow? Color Schemes for Improved Data
    ##  Graphics. EOS Transactions of the American Geophysical Union, 85(40), 385.
    ## See also http://geog.uoregon.edu/datagraphics/color_scales.htm

    ## Plotting posterior probability trace

    ## /home/jpuritz/eems_calcar/microsat-D768-chain1

    ## /home/jpuritz/eems_calcar/microsat-D768-chain2

    ## /home/jpuritz/eems_calcar/microsat-D768-chain3

    ## Plotting average dissimilarities within and between demes

    ## /home/jpuritz/eems_calcar/microsat-D768-chain1

    ## /home/jpuritz/eems_calcar/microsat-D768-chain2

    ## /home/jpuritz/eems_calcar/microsat-D768-chain3
