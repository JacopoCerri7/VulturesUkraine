# VulturesUkraine
A repository with the reproducible software code for the analysis of the movement of a cinereous vulture over the active war zone in Ukraine, during March 2022

## Information

The following folder contains the reproducible data and software code for the study entitled "GPS telemetry reveals changes in the movement patterns of a cinereous vulture (Aegypius monachus) flying over an active war zone". Here, there are two files.

### Software code
"ReproducibleCode_20261029.R" contains the software code for replicating the analyses in R. This script has to be executed in R (https://www.r-project.org/). Relevant packages are highlighted in the initial section of the script.


### Dataset

ReproducibleDataset_20261029.Rdata" contains the GPS locations of the tracked individual, with the following attributes:

- "timestamp", or the timestamp of the GPS location (UTC time zone)
- "date", or the date of each GPS location
- "lon", or the longitude of each GPS location (epsg: 4326)
- "lat", or the longitude of each GPS location (epsg: 4326)
- "lon.utm", or the longitude of each GPS location (epsg: 3035)
- "lat.utm", or the longitude of each GPS location (epsg: 3035)
- "ta", or turning angles (in degrees)
- "sl", or the length of the step (in metres) since the previous  location
- "difftime", or the time elapsed (in minutes) since the previous location
- "speed", or the ratio between step length and time
- "elev.msl", or the elevation according to the altimeter (in metres)
- "dem", or the elevation according to the DEM (Amazon AWS Global DEM, in metres)
- "h", or the flying height, calculated as the difference between "elev.msl" and "dem"
- "cross", or the flight over Ukraine (first, second, third, see the manuscript)
- "dist.front.abs", or the distance from the active front line (in metres)
- "dist.front.abs.std", or the distance from the active front line (converted to Z-scores)
- "sl.std", or the length of the step since the previous location (converted to Z-scores)
- "h.std", or the flying height (converted to Z-scores)
- "ta.rad", or turning angles (in radiants)

For further information please contact Jacopo Cerri (Mammal Research Institute, Polish Academy of Sciences) at j.cerri@ibs.bialowieza.pl
