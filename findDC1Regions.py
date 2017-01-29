import numpy as np
from plotFunctions import *
from intermediates import *

def findDC1Regions(coaddBundle, dbpath, plotTestPlots= True,
                   filterBand= 'i', threshold= 0.0001, nside= 256,
                   returnAll= False):
    """

    Find candidate regions for DC1 (based on how wildly the average depth in the region 
    differs from the survey median depth).

    Returns a bunch of stuff ...

    Required Parameters
    -------------------
      * coaddBundle: dict: dictionary with keys= observing strategy names, pointing to corresponding
                           to a metricBundle object.
           NOTE: coaddBundle should have always have NoDither. If want to find regions based on a dithered
                 survey, the bundle should have the dithered data ALONG WITH the undithered one.
      * dbpath: str: path to the OpSim database.

    Optional Parameters
    -------------------
      * plotTestPlots: bool: set to False if dont want to plot things for debugging/testing code.
                             Default: True
      * filterBand: str: filter to consider. Default: 'i'
      * threshold: float: region will be considered good if average depth in the region is within the
                          threshold of survey median depth. Default: 0.0001
      * nside: int: HEALPix resolution parameter. Defaut: 256
      * returnAll: bool: set to True to get (things needed to find chips.. ):
                         [focusDither, output_rect, output_disc, simdata, pixels_in_FOV,
                                                   simdataIndex_for_pixel, pixelNum, pixRA, pixDec]
                         Default: False returns:
                         [focusDither, output_rect, output_disc]
    """
    FOV_radius= 0.0305
    printProgress('Getting RA, Dec for HEALPix pixels ...', highlight= True)
    pixelNum= getSurveyHEALPixRADec(coaddBundle)   # each output is a dicitonary.
    
    printProgress('Getting simdata ...', highlight= True)
    simdata= getSimData(dbpath, filterBand)    # contains fieldID, fieldRA, fieldDec, rotSkyPos, expMJD, ditheredRA, ditheredDec
    
    printProgress('Getting pixels_in_FOV ...', highlight= True)
    pixels_in_FOV= getFOVsHEALPixReln(pixelNum, simdata, nside= nside) # each output is a dicitonary.

    #########################################################################################################################
    if plotTestPlots:
        # plot sample FOV for each dither strategy (NoDither and anohter one)
        printProgress('Code test: Plotting using plotFOV ...', highlight= True)
        IDsToTestWith= [1421, 1365, 1447]
        plotFOV(coaddBundle, pixels_in_FOV, IDsToTestWith, filterBand)
    
        printProgress('Code test: Plotting regions using buildAndPlotRegion ...', highlight= True)
        fID= 1421
        printProgress('Code test: buildAndPlotRegion: Plots with a NON-DISC region of interest ...')
        buildAndPlotRegion(fID, simdata, coaddBundle, FOV_radius, nside= nside,
                           pixels_in_FOV= pixels_in_FOV,
                           disc= False)
        printProgress('Code test: buildAndPlotRegion: Plots with a CIRCULAR region of interest ...')
        buildAndPlotRegion(fID, simdata, coaddBundle, FOV_radius, nside= nside,
                            pixels_in_FOV= pixels_in_FOV,
                           disc= True)
    ########################################################################################################################
    printProgress('Finding good regions ...', highlight= True)
    
    if (len(coaddBundle.keys())==1): # only NoDither provided.
        focusDither= 'NoDither'
    else:
        for dither in coaddBundle.keys():
            if (dither != 'NoDither'): focusDither= dither
            
    printProgress('Finding good regions with threshold= %f using %s' % (threshold, focusDither), highlight= True)
    output_nonDisc= findGoodRegions(focusDither, simdata, coaddBundle, FOV_radius, pixels_in_FOV,
                                 allIDs= True, nside= nside,
                                 disc= False, threshold= threshold)
    output_disc= findGoodRegions(focusDither, simdata, coaddBundle, FOV_radius, pixels_in_FOV,
                                 allIDs= True, nside= nside,
                                 disc= True, threshold= threshold)

    printProgress('Plotting good regions with threshold= %f using %s' % (threshold, focusDither), highlight= True)
    printProgress('Rectangular regions (using plotRegion):')
    plotRegion(coaddBundle,focusDither, pixels_in_FOV,
               output_nonDisc['goodFiducialIDs'], output_nonDisc['regionPixels'], filterBand= filterBand)
    printProgress('Cicular regions (using plotRegion):')
    plotRegion(coaddBundle, focusDither, pixels_in_FOV,
               output_disc['goodFiducialIDs'], output_disc['regionPixels'], filterBand= filterBand)
    if returnAll:
        return [focusDither, output_nonDisc, output_disc, simdata, pixels_in_FOV, pixelNum, pixRA, pixDec]
    else:
        return [focusDither, output_nonDisc, output_disc]
