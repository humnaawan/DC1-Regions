# DC1-Regions

This repo hosts the code/output for finding the DC1 region to run through imSim/PhoSim.

The idea is to find a region in the wide-fast-deep (WFD; main) survey region that is 'typical' -- defined in terms of having typical 5\sigma coadded depth (in comparison with the median depth across the survey); defined in order to avoid simulating say a Deep Drilling Field, which wouldnt be representative of the majority of the survey and the LSS studies we can/will carry out.

The region considered is based on 4 contiguous undithered fields of view configured so that the region contains (at least) two triple points; these fields of view are indexed with respective fieldIDs. We consider two geometries that could circumscribe the 4 contiguous fields; see `DC1 Region Types.pdf` for the two geometries. Since such a configuration is really only fixed in an undithered survey and we want to base our DC1 region on a dithered survey, we use the regions defined by the undithered configuration and check for a typical coadded depth in the region resulting from a dithered survey.

There are two main functions in this repository:
## findDC1Regions
This functions finds possible candidate regions for DC1; depends on intermediates.py and plotFunctions.py. The basic flow is:

1. Find RA, Dec that each HEALPix pixels corresponds to.

2. Find out which HEALPix pixels correspond to which fieldID.

3. Run over all the undithered fieldIDs and find the HEALPix pixels that either lie within a circular area (using healpy.query_disc; first geometry in `DC1 Region Types.pdf`). 

4. Take the coadded depth map (for user-specified observing strategy, undithered or a dithered one), consider all the 'regions', and calculate the difference in the depth and the range of depth in the region.

Here, difference= abs(median depth in the region - median depth across the survey) and range= abs(min depth - max depth in the region). 
#### A region is considered 'good' if difference<tolerance, where tolerance is user-defined and scatter is small.

The function returns the fiducial fieldID (i.e. fID on which the region is based on), the pixel numbers in each good region, the fieldIDs that these pixels belong to, etc. Also we get back the scatter value to see what regions are better.

#### Note: fieldRA, fieldDec (corresponding to undithered positions) remain constant and are needed for this region-finding step in order to connect various HEALPix pixels with fieldID; dithered telescope pointings don't play any role here. The user-inputted coadded depth (for either undithered or dithered survey) comes into play in Step 4 above.

## findDC1Chips
The purpose of this function is to find all the chips that are observed in a region through the ten year survey. It takes a fieldID, find the pixels in the region(second figure in `DC1 Region Types.pdf`) and a large enough Nside (=512) such that HEALPix pixel size is less than the chip size.

#### Note: This is the function whose output *would change* if new afterbuner output is available since the dithered columns added by the new afterburner are different than the ones from before. The dithered OpSim columns come into play (assuming we want to find the chips for a dithered survey) alongside rotSkyPos, expMJD, etc characterizing the telescope position at any given visit. 

-------
See the notebooks folder for iPython notebooks runnings these two functions. Output folder contains the output from `findDC1Chips`.

-------
## To Do?
1. Find a way to check the list of chips per visit lies in the intended region.