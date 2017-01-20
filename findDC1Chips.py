import numpy as np
from lsst.obs.lsstSim import LsstSimMapper
from lsst.sims.utils import ObservationMetaData
from lsst.sims.coordUtils import chipNameFromRaDec
    
def findDC1Chips(dither, regionPixels, simdataIndex_for_pixel, pixelNum, pixRA, pixDec, simdata):
    """

    Find the chips that are used in the region.

    Required Parameters
    -------------------
      * dither: str: dither strategy to focus on.
      * regionPixels: list: list of pixel numbers in the region to plot.
      * simdataIndex_for_pixel: dict: dictionary with keys= dither strategy. Each key points to a dictionary
                                      with keys= pixel number, pointing to the list of indices corresponding
                                      to that pixel in simdata array.
      * pixelNum, pixelRA, pixelDec: output of getSurveyHEALPixRADec
      * simdata: output of getSimData (must have fieldRA, fieldDec, ditheredRA, ditheredDec, rotSkyPos, expMJD)

    """
    camera = LsstSimMapper().camera
    
    chipNames= []
    for p, pixel in enumerate(regionPixels):
        simdataInds= simdataIndex_for_pixel[dither][pixel]['idxs']
        
        for index in simdataInds:
            if (dither=='NoDither'):
                pointingRA= np.degrees(simdata[index]['fieldRA'])
                pointingDec= np.degrees(simdata[index]['fieldDec'])
            else:
                pointingRA= np.degrees(simdata[index]['ditheredRA'])
                pointingDec= np.degrees(simdata[index]['ditheredDec'])
            rotSkyPos= np.degrees(simdata[index]['rotSkyPos'])
            expMJD= simdata[index]['expMJD']
            
            obs = ObservationMetaData(pointingRA= pointingRA, pointingDec= pointingDec,
                                      rotSkyPos= rotSkyPos, mjd= expMJD)
            i= np.where(pixel==pixelNum[dither])[0][0]
            chipNames.append(chipNameFromRaDec(np.degrees(pixRA[dither][i]),
                                               np.degrees(pixDec[dither][i]),
                                               camera=camera, obs_metadata=obs))
    return np.unique(chipNames)
