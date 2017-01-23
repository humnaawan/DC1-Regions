import healpy as hp
import matplotlib.pyplot as plt
import copy
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import sys
import time
import numpy as np

__all__= ['printProgress', 'getSurveyHEALPixRADec', 'getSimData',
          'getFOVsHEALPixReln', 'enclosingPolygon', 'findRegionPixels',
          'findGoodRegions', 'findRegionFOVs' ]

def printProgress(whatToPrint, highlight= False, newLine= True):
    """
    Print statements while making sure statements are printed before the next thing starts.

    Required Parameter
    ------------------
    whatToPrint: str: whatever to print

    Optional Parameters
    -------------------
    * highlight: bool: set to True to add an extra line of ## to highlight the print statement.
                       Default: False
    * newline: bool: set to False to not add a new line.
                     Default: True

    """
    append= ''
    if highlight: append= '\n############################################'
    if newLine: print append + '\n## ' + whatToPrint
    else: print '## ' + whatToPrint
    sys.stdout.flush()
    time.sleep(1.0)
        
def getSurveyHEALPixRADec(coaddBundle):
    """

    Get the RA, Dec (in radians) corresponding to each HEALPix pixel.
    Method returns three dictionaries with keys= keys in coaddBundle: pixelNumber, pixelRA, pixelDec

    Required Parameter:  coaddBundle: dict: dictionary with keys= observing strategy names, pointing to
                                            corresponding to a metricBundle object.

    """
    # create dictionaries giving pixelNumbers and their correspondong RA, Dec for all dither strategies.
    # need to worry about each strategy separately since the mask is generally different.
    pixelNum, pixRA, pixDec= {}, {}, {}
    for dither in coaddBundle:
        pixelNum[dither], pixRA[dither], pixDec[dither]= [], [], []
        for pix in range(len(coaddBundle[dither].slicer)):
            if not coaddBundle[dither].metricValues.mask[pix]:   # only consider the unmasked pixels
                temp= coaddBundle[dither].slicer._pix2radec(pix)    # radians returned
                pixelNum[dither].append(pix)
                pixRA[dither].append(temp[0])
                pixDec[dither].append(temp[1])
    return [pixelNum, pixRA, pixDec]

def getSimData(dbpath, filterBand, extraCols= []):
    """

    Get OpSim data columns (for WFD).

    Required Parameters
    -------------------
      * dbpath: str: path to the OpSim database.
      * filterBand: str: filter to consider, e.g. 'r'

    Optional Parameters
    -------------------
      * extraCols: list of str: list of additional columns needed from the database.
                                e.g. ['expDate', 'obsHistID']

    """
    # get the columns we care about in simdata.
    import lsst.sims.maf.db as db
    import lsst.sims.maf.utils as mafUtils
    #dbfile = path+'minion_1016_sqlite.db'
    opsdb = db.OpsimDatabase(dbpath)
    propIds, propTags = opsdb.fetchPropInfo()
    wfdWhere = mafUtils.createSQLWhere('WFD', propTags)
    sqlconstraint= wfdWhere + ' and filter=="' + filterBand + '"'
    colnames = ['fieldID', 'fieldRA', 'fieldDec', 'rotSkyPos', 'expMJD', 'ditheredRA', 'ditheredDec'] + extraCols
    simdata = opsdb.fetchMetricData(colnames, sqlconstraint)
    
    return simdata
    
def getFOVsHEALPixReln(pixelNum, pixRA, pixDec, simdata):
    """

    Finds the correspondences between different FOVs and HEALPix pixels.
    Returns two dictionaries (keys= keys in pixelNum): For each keys:
    1. pixels_in_FOV is a dictionary of keys= fieldID, pointing to the list of HEALPix pixels
                      corresponding to the FOV.
    2. simdataIndex_for_pixel is  a dictionary with keys= pixel number, pointing to the list of
                              indices corresponding to that pixel in simdata array.

    Required Parameters
    -------------------
      * pixelNum, pixelRA, pixelDec: output of getSurveyHEALPixRADec
      * simdata: output of getSimData

    """
    # each of pixelNum, pixRA, pixDec is a dicitonary.
    import lsst.sims.maf.slicers as slicers
    pixels_in_FOV= {}
    simdataIndex_for_pixel= {}

    for dither in pixelNum:
        pixels_in_FOV[dither]= {}
        simdataIndex_for_pixel[dither]= {}
        slicer = slicers.UserPointsSlicer(ra= pixRA[dither], dec= pixDec[dither])   # inputting radians (<=> HEALPix survey pixels)
        slicer.setupSlicer(simdata)
    
        for i in range(len(slicer)):  # running over only the sky pixels
            ind = slicer._sliceSimData(i)
            simdataIndex_for_pixel[dither][pixelNum[dither][i]]= ind
            ids = simdata[ind['idxs']]['fieldID']   # fieldIDs corresponding to pixelNum[i]
            for uniqID in np.unique(ids):
                key= uniqID
                if key not in pixels_in_FOV[dither].keys():
                    pixels_in_FOV[dither][key]= []
                pixels_in_FOV[dither][key].append(pixelNum[dither][i])
    print 'Number of fieldIDs in pixel_in_FOV for %s: %f' %(dither, len(pixels_in_FOV[dither].keys()))
    return [pixels_in_FOV, simdataIndex_for_pixel]


def enclosingPolygon(radius, fieldRA, fieldDec):
    """

    Returns the corners of the rectangular region to input into query_polygon.

    Required Parameters
    -------------------
      * radius: float: radius of the FOV in radians.
      * fieldRA: float: RA (in radians) of the FOV center on which to base the rectangle.
      * fieldDec: float: Dec (in radians) of the FOV center on which to base the rectangle.

    """
    corners= np.zeros(shape=(4,3))
    
    x_pt= fieldRA+radius
    y_pt= fieldDec-np.sqrt(3)*radius
    c = SkyCoord(ra=x_pt*u.radian, dec=y_pt*u.radian)
    corners[0,]= c.cartesian.xyz
    
    x_pt= fieldRA+radius
    y_pt= fieldDec+np.sqrt(3)*radius
    c = SkyCoord(ra=x_pt*u.radian, dec=y_pt*u.radian)
    corners[1,]= c.cartesian.xyz
    
    x_pt= fieldRA-4*radius
    y_pt= fieldDec+np.sqrt(3)*radius
    c = SkyCoord(ra=x_pt*u.radian, dec=y_pt*u.radian)
    corners[2,]= c.cartesian.xyz
    
    x_pt= fieldRA-4*radius
    y_pt= fieldDec-np.sqrt(3)*radius
    c = SkyCoord(ra=x_pt*u.radian, dec=y_pt*u.radian)
    corners[3,]= c.cartesian.xyz

    return corners

def findRegionPixels(fID, simdata, nside, disc, FOV_radius):
    """

    Find the region (disc or rectangular) based on the specified field ID.

    Required Parameters
    -------------------
      * fID: int: fieldID for the FOV on which to base the region.
      * simdata: np.array: array containing OpSim columns (must have fieldID, fieldRA, fieldDec).
      * nside: int: HEALPix resolution parameter.
      * disc: bool: set to true if want disc-like region; False for rectangular.
      * FOV_radius: float: radius of the FOV in radians.

    """
    ind= np.where(simdata[:]['fieldID']== fID)[0]
    # fieldRA, fieldDec remain fixed for NoDither; dont change with expMJD.
    # use as the 'center' of the enclosing region (disc or rectangle).
    fixedRA= simdata[ind[0]]['fieldRA']
    fixedDec= simdata[ind[0]]['fieldDec']
    if not disc:
        centralRA, centralDec= fixedRA, fixedDec
        corners= enclosingPolygon(FOV_radius, centralRA, centralDec)
        diskPixels= hp.query_polygon(nside, corners)    # HEALpixel numbers
    else:
        centralRA, centralDec= fixedRA, fixedDec-FOV_radius*np.sqrt(3)/2.
        c = SkyCoord(ra=centralRA*u.radian, dec= centralDec*u.radian)
        diskPixels= hp.query_disc(nside= nside, vec=c.cartesian.xyz, radius= 2.5*FOV_radius)

    return [centralRA, centralDec, diskPixels]

def findRegionFOVs(regionPixels, dither, simdataIndex_for_pixel, simdata):
    """

    Find the FOVs that corresponds to any HEALPix pixels in the region.

    Required Parameters
    -------------------
      * regionPixels: array: array containing the HEALPIx pixel numbers in the region of interest.
      * dither: str: dither strategy to focus on.
      * simdataIndex_for_pixel: dict: dictionary with keys= dither strategy. Each key points to a dictionary
                                      with keys= pixel number, pointing to the list of indices corresponding
                                      to that pixel in simdata array.
      * simdata: np.array: array containing OpSim columns (must have fieldID for here).

    """
    idList= []
    for p in regionPixels:
        ind= simdataIndex_for_pixel[dither][p]
        ids = simdata[ind['idxs']]['fieldID']   # fieldIDs corresponding to pixelNum[i]
        uniqID= np.unique(ids)
        idList+= list(uniqID)
    return np.unique(idList)

def findGoodRegions(simdata, coaddBundle, surveyMedianDepth, FOV_radius, pixels_in_FOV,
                    nside= 256, threshold= 0.01,
                    allIDs= True, IDsToTestWith= [],
                    disc= False,  raRange= [-180,180], decRange= [-70,10]):
    """

    Find good regions (i.e. regions with average depth within specified threshold of survey
    median depth).

    Returns: arrays: good pixel numbers, good IDs, difference between mean depth in the region and
    median survey depth, abs(max-min) region depth, center RA, center Dec

    Required Parameters
    -------------------
      * simdata: np.array: array containing OpSim columns (must have fieldID, fieldRA, fieldDec).
      * coaddBundle: dict: dictionary with keys= observing strategy names, pointing to corresponding
                           to a metricBundle object.
      * surveyMedianDepth: float: median depth to base the 'goodness' on.
      * FOV_radius: float: radius of the FOV in radians.
      * pixels_in_FOV: dict: dictionary with keys= field ID, pointing to the list of HEALPix pixels
                             that fall in the FOV.

    Optional Parameters
    -------------------
      * nside: int: HEALPix resolution parameter. Defaut: 256
      * threshold: float: region will be considered good if average depth in the region is within the
                          threshold of survey median depth. Default: 0.01
      * allIDs: bool: set to False to consider only a subset of FOVs; will plot things out. Default: True
      * IDsToTestWith: list: list of fieldIDs to consider. Default: []
      * disc: bool: set to True if want disc-like region; False for rectangular. Default: False
      * raRange: list: constraint on the RA in cartview. Default: [-180,180]
      * decRange: list: constraint on the Dec in cartview. Default: [-70,10]

    """
    # a region is 'good' if abs(typicalDepth in the region -surveyMedianDepth)<threshold
    goodIDs= []
    goodPixelNums= []
    scatterInDepth= []
    diffMeanMedian= []
    centerRA= []
    centerDec= []
    
    focusDither= surveyMedianDepth.keys()[0]
    considerIDs= pixels_in_FOV[focusDither].keys()
    if not allIDs: considerIDs= IDsToTestWith
        
    for ID in considerIDs:
        centralRA, centralDec, diskPixels= findRegionPixels(ID, simdata, nside, disc, FOV_radius)

        typicalDepth= np.mean(coaddBundle[focusDither].metricValues.data[diskPixels])
        diff= abs(typicalDepth-surveyMedianDepth[focusDither])
        if (diff<threshold): 
            goodIDs.append(ID)
            goodPixelNums.append(diskPixels)
            diffMeanMedian.append(diff)
            scatterInDepth.append(abs(max(coaddBundle[focusDither].metricValues.data[diskPixels])-min(coaddBundle[focusDither].metricValues.data[diskPixels])))
            centerRA.append(centralRA)
            centerDec.append(centralDec)
            
        if not allIDs:
            check= copy.deepcopy(coaddBundle[[focusDither]])
            check.metricValues.data[:]= 0
            check.metricValues.data[diskPixels]= 1000.
            check.metricValues.data[pixels_in_FOV[focusDither][ID]]= 26.4
            
            plt.clf()
            hp.cartview(check.metricValues.filled(check.slicer.badval), 
                        flip='astro', rot=(0,0,0) ,
                        lonra= raRange, latra= decRange,
                        min= 26.3, max= 26.5, title= '', cbar=False)
            hp.graticule(dpar=20, dmer=20, verbose=False)
            plt.title('fID ' + str(ID) , size= 22)
            ax = plt.gca()
            im = ax.get_images()[0]
            fig= plt.gcf()
            cbaxes = fig.add_axes([0.1, 0.25, 0.8, 0.04]) # [left, bottom, width, height]
            cb = plt.colorbar(im,  orientation='horizontal',
                              format= '%.1f', cax = cbaxes) 
            #cb.set_label(str('i-Band Coadded Depth'), fontsize=18)
            #cb.ax.tick_params(labelsize= 18)
            plt.show()
            
    return [np.array(goodPixelNums), np.array(goodIDs), np.array(diffMeanMedian), np.array(scatterInDepth), np.array(centerRA), np.array(centerDec)]


