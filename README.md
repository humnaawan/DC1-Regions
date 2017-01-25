# DC1-Regions

This repo hosts the code/output for finding the DC1 region to run through imSim/PhoSim.

The idea is to find a region in the wide-fast-deep (WFD; main) survey region that is 'typical' -- defined in terms of having typical 5\sigma coadded depth (in comparison with the median depth across the survey); defined in order to avoid simulating say a Deep Drilling Field, which wouldnt be representative of the majority of the survey and the LSS studies we can/will carry out.

The region considered is based on 4 contiguous undithered FOVs configured so that the region contains two triple points. We consider two geometries that could circumscribe the 4 contiguous fields; see `DC1 Region Types.pdf` for the two geometries. Since such a configuration is really only fixed in an undithered survey and we want to base our DC1 region on a dithered survey, we use the regions defined by the undithered configuration and check for a typical coadded depth in the region resulting from a dithered survey.

There are two main functions in this repository:
## findDC1Regions
This functions finds possible candidate regions for DC1; depends on intermediates.py and plotFunctions.py. The basic flow is:

1. Find RA, Dec that each HEALPix pixels corresponds to.

2. Find out which HEALPix pixels correspond to which FOV.

3. Run over all the undithered FOVs and find the HEALPix pixels that either lie within a rectangular region or a circular one (using healpy.query_disc). 

4. Take the coadded depth map (for user-specified observing strategy, undithered or a dithered one), consider all the 'regions', and calculate the difference in the depth and the scatter.

Here, difference= abs(np.mean depth in the region - median depth across the survey) and scatter= abs(min depth - max depth in the region). 
#### A region is considered 'good' if difference<tolerance, where tolerance is user-defined.

The function returns the fiducial fieldID (i.e. fID on which the region is based on), the pixel numbers in each good region, the FOVs that these pixels belong to, etc. Also we get back the scatter value to see what regions are better.

#### Note: Only undithered data columns are used from the OpSim data for this region-finding step, alongside the user-inputted coadded depth (for either undithered or dithered survey) which comes into play in Step 4 above.

## findDC1Chips
The purpose of this function is to find all the chips that are observed in a region through the ten year survey. It takes a fieldID, find the pixels in the region and a large enough Nside (=512) such that HEALPix pixel size is less than the chip size.

#### Note: This is the function whose output would change if new afterbuner output is available since the dithered columns added by the new afterburner are different than the ones from before. The dithered OpSim columns come into play (assuming we want to find the chips for a dithered survey) alongside rotSkyPos, expMJD, etc characterizing the telescope position at any given visit. 

-------
See the notebooks folder for iPython notebooks runnings these two functions. Output folder contains the output from `findDC1Chips`.

-------
## To Do
1. Get all the chips when get the warning with a pixel landing in more than one pixel. Current version of chipNameFromRaDec wont return all so need a way to know if there are exceptions raised.