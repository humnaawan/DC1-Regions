import healpy as hp
import matplotlib.pyplot as plt
import copy
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic
import pandas as pd
import sys
import time
import numpy as np

__all__= ['printProgress', 'getSurveyHEALPixRADec', 'getSimData',
          'getFOVsHEALPixReln', 'enclosingPolygon', 'findRegionPixels',
          'findContigFOVs', 'findRegionFOVs', 'findGoodRegions']

def printProgress(whatToPrint, highlight= False, newLine= True):
    """
    Print statements while making sure statements are printed before the next thing starts.

    Required Parameter
    ------------------
    * whatToPrint: str: whatever to print

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

    Required Parameter
    ------------------
    * coaddBundle: dict: dictionary with keys= observing strategy names, pointing to
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

def getSimData(dbpath, filterBand, extraCols= [], newAfterburner= False):
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
    * newAfterburner: bool: set to True if the opsim database is using the
                            new afterbuner. Default: False

    """
    # get the columns we care about in simdata.
    import lsst.sims.maf.db as db
    import lsst.sims.maf.utils as mafUtils
    #dbfile = path+'minion_1016_sqlite.db'
    opsdb = db.OpsimDatabase(dbpath)
    propIds, propTags = opsdb.fetchPropInfo()
    wfdWhere = mafUtils.createSQLWhere('WFD', propTags)
    sqlconstraint= wfdWhere + ' and filter=="' + filterBand + '"'
    if newAfterburner:
        colnames = ['fieldID', 'fieldRA', 'fieldDec', 'rotSkyPos', 'expMJD',
                    'hexDitherPerNightRA', 'hexDitherPerNightDec',
                    'randomDitherFieldPerVisitRA', 'randomDitherFieldPerVisitDec'] + extraCols
    else:
        colnames = ['fieldID', 'fieldRA', 'fieldDec', 'rotSkyPos', 'expMJD', 'ditheredRA', 'ditheredDec'] + extraCols
    simdata = opsdb.fetchMetricData(colnames, sqlconstraint)
    
    return simdata
    
def getFOVsHEALPixReln(pixelNum, simdata, nside= 256):
    """

    Finds the correspondences between different FOVs and HEALPix pixels.
    Returns a dictionary pixels_in_FOV (keys= keys in pixelNum; observing stragies):
    each key points to a dictionary with keys= fieldID, pointing to the list of HEALPix pixels
    corresponding to the FOV.
 

    Required Parameters
    -------------------
    * pixelNum: dict: keys= observing strategies; pointing to pixelNumbers in the survey region.
    * simdata: np.array: array containing OpSim columns (must have fieldID for this function).

    Optional Parameters
    --------------------
    * nside: int: HEALPix resolution parameters. Default: 256

    """
    import lsst.sims.maf.slicers as slicers
    pixels_in_FOV= {}

    for dither in pixelNum:
        pixels_in_FOV[dither]= {}
        slicer= slicers.HealpixSlicer(nside= nside)
        slicer.setupSlicer(simdata)

        for pixel in pixelNum[dither]:
            indObsInPixel = slicer._sliceSimData(pixel)
            ids = simdata[indObsInPixel['idxs']]['fieldID']   # fieldIDs corresponding to pixel
            for uniqID in np.unique(ids):
                key= uniqID
                if key not in pixels_in_FOV[dither].keys():
                    pixels_in_FOV[dither][key]= []
                pixels_in_FOV[dither][key].append(pixel)
    print 'Number of fieldIDs in pixel_in_FOV for %s: %f' %(dither, len(pixels_in_FOV[dither].keys()))
    return pixels_in_FOV

def enclosingPolygon(radius, fieldRA, fieldDec):
    """

    Returns the corners of the rectangular region to input into query_polygon.

    Required Parameters
    -------------------
    * radius: float: radius of the FOV in radians.
    * fieldRA: float: RA (in radians) of the FOV center on which to base the rectangle.
    * fieldDec: float: Dec (in radians) of the FOV center on which to base the rectangle.

    """
    def returnXYZ(ra, dec):
        c = SkyCoord(ra=ra*u.radian, dec=dec*u.radian)
        return c.cartesian.xyz
 
    corners= np.zeros(shape=(4,3))
    corners[0,]= returnXYZ(fieldRA+radius, fieldDec-np.sqrt(3)*radius)
    corners[1,]= returnXYZ(fieldRA+radius, fieldDec+np.sqrt(3)*radius)
    corners[2,]= returnXYZ(fieldRA-4*radius, fieldDec+np.sqrt(3)*radius)
    corners[3,]= returnXYZ(fieldRA-4*radius, fieldDec-np.sqrt(3)*radius)

    return corners

def findRegionPixels(fiducialID, simdata, nside, disc, FOV_radius):
    """

    Find the region (disc or rectangular) based on the specified field ID.

    Required Parameters
    -------------------
    * fiducialID: int: fieldID for the FOV on which to base the region.
    * simdata: np.array: array containing OpSim columns (must have fieldID, fieldRA, fieldDec).
    * nside: int: HEALPix resolution parameter.
    * disc: bool: set to true if want disc-like region; False for rectangular.
    * FOV_radius: float: radius of the FOV in radians.

    """
    ind= np.where(simdata[:]['fieldID']== fiducialID)[0]
    # fieldRA, fieldDec remain fixed for NoDither; dont change with expMJD.
    # use as the 'center' of the enclosing region (disc or rectangle).
    fixedRA= simdata[ind[0]]['fieldRA']
    fixedDec= simdata[ind[0]]['fieldDec']
    if not disc:
        fiducialRA, fiducialDec= fixedRA, fixedDec
        corners= enclosingPolygon(FOV_radius, fiducialRA, fiducialDec)
        diskPixels= hp.query_polygon(nside, corners)    # HEALpixel numbers
    else:
        fiducialRA, fiducialDec= fixedRA, fixedDec-FOV_radius*np.sqrt(3)/2.
        c = SkyCoord(ra=fiducialRA*u.radian, dec= fiducialDec*u.radian)
        diskPixels= hp.query_disc(nside= nside, vec=c.cartesian.xyz, radius= 2.5*FOV_radius)

    return [fiducialRA, fiducialDec, diskPixels]

def findContigFOVs(fiducialRA, fiducialDec, fiducialID, FOV_radius, simdata,
                   disc= True, nside= 256):
    """

    Find contiguous fields in an undithered survey based on the specified ID;
    basically will return the IDs of the three fields that surrounding the fiducialID in the
    arrangement set up for the DC1 region.

    Returns the IDs of the three neighbors alongwith the fiducialID.

    Required Parameters
    -------------------
    * fiducialRA: float: RA of fiducial center of the region (radians).
    * fiducialDec: float: Dec of fiducial center of the region (radians).
    * fiducialID: int: fieldID on which the region is based on.
    * FOV_radius: float: radius of the FOV (radians)
    * simdata: np.array: array containing OpSim columns (must have fieldID for this function to work).
    
    Optional Parameters
    -------------------
    * disc: bool: set to True if want disc-like region; False for rectangular. Default: True
    * nside: int: HEALPix resolution parameters. Default: 256

    """
    import lsst.sims.maf.slicers as slicers
    slicer= slicers.HealpixSlicer(nside= nside)
    slicer.setupSlicer(simdata)
    
    dither= 'NoDither'
    contList= []
    contList.append(fiducialID)
    
    if disc:
        # FOV below
        ra= fiducialRA
        dec= fiducialDec-FOV_radius*np.sqrt(3)/2.
        p= hp.ang2pix(nside, np.pi/2.0-dec, ra)
        ind = slicer._sliceSimData(p)
        ids = simdata[ind['idxs']]['fieldID']   # fieldIDs corresponding to pixelNum[i]
        if len(np.unique(ids))>1:
            print 'Something is wrong. Should have only one FOV corresponding to the central pixel.'
            return
        contList.append(ids[0])
        
        ra= fiducialRA-FOV_radius*3/2.
        dec= fiducialDec
        p= hp.ang2pix(nside, np.pi/2.0-dec, ra)
        ind = slicer._sliceSimData(p)
        ids = simdata[ind['idxs']]['fieldID']   # fieldIDs corresponding to pixelNum[i]
        if len(np.unique(ids))>1:
            print 'Something is wrong. Should have only one FOV corresponding to the central pixel.'
            return
        contList.append(ids[0])
        
        ra= fiducialRA+FOV_radius*3/2.
        dec= fiducialDec
        p= hp.ang2pix(nside, np.pi/2.0-dec, ra)
        ind = slicer._sliceSimData(p)
        ids = simdata[ind['idxs']]['fieldID']   # fieldIDs corresponding to pixelNum[i]
        if len(np.unique(ids))>1:
            print 'Something is wrong. Should have only one FOV corresponding to the central pixel.'
            return
        contList.append(ids[0])

        return contList
    else:
        print 'Not developed the function for anything but disc-like region.'
        return

def findRegionFOVs(regionPixels, dither, simdata, nside= 256):
    """

    Find the FOVs that corresponds to any HEALPix pixels in the region.

    Required Parameters
    -------------------
    * regionPixels: array: array containing the HEALPIx pixel numbers in the region of interest.
    * dither: str: dither strategy to focus on.
    * simdata: np.array: array containing OpSim columns (must have fieldID for here).
    
    Optional Parameters
    -------------------
    * nside: int: HEALPix resolution parameters. Default: 256

    """
    import lsst.sims.maf.slicers as slicers
    slicer= slicers.HealpixSlicer(nside= nside)
    slicer.setupSlicer(simdata)
    
    idList= []
    for p in regionPixels:
        ind = slicer._sliceSimData(p)
        ids = simdata[ind['idxs']]['fieldID']   # fieldIDs corresponding to pixelNum[i]
        uniqID= np.unique(ids)
        idList+= list(uniqID)
    return np.unique(idList)

def findGoodRegions(focusDither, simdata, coaddBundle, FOV_radius, pixels_in_FOV,
                    nside= 256, threshold= 0.01,
                    allIDs= True, IDsToTestWith= [],
                    disc= False,  raRange= [-180,180], decRange= [-70,10]):
    """

    Find good regions (i.e. regions with average depth within specified threshold of survey
    median depth).

    Returns: arrays: good pixel numbers, good IDs, difference between mean depth in the region and
    median survey depth, abs(max-min) region depth, fiducial RA, fidicual Dec

    Required Parameters
    -------------------
    * focusDither: str: the observing strategy to focus on.
                        Options: 'NoDither', 'SequentialHexDitherPerNight'.
    * simdata: np.array: array containing OpSim columns (must have fieldID, fieldRA, fieldDec).
    * coaddBundle: dict: dictionary with keys= observing strategy names, pointing to corresponding
                         to a metricBundle object.
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
    goodCenterIDs, goodPixelNums, scatterInDepth, diffMeanMedian, fiducialRAs, fiducialDecs, contigIDs, fiducialGalacticLat= [], [], [], [], [], [], [], []
    
    inSurvey= np.where(coaddBundle[focusDither].metricValues.mask==False)[0]
    surveyMedianDepth= np.median(coaddBundle[focusDither].metricValues.data[inSurvey])
    printProgress('Mean survey depth for %s: %f'% (focusDither, surveyMedianDepth))
    
    considerIDs= pixels_in_FOV[focusDither].keys()
    if not allIDs: considerIDs= IDsToTestWith
        
    for ID in considerIDs:
        fiducialRA, fiducialDec, diskPixels= findRegionPixels(ID, simdata, nside, disc, FOV_radius)

        typicalDepth= np.mean(coaddBundle[focusDither].metricValues.data[diskPixels])
        diff= abs(typicalDepth-surveyMedianDepth)
        
        if (diff<threshold):
            goodCenterIDs.append(ID)
            goodPixelNums.append(diskPixels)
            diffMeanMedian.append(diff)
            scatterInDepth.append(abs(max(coaddBundle[focusDither].metricValues.data[diskPixels])-min(coaddBundle[focusDither].metricValues.data[diskPixels])))
            fiducialRAs.append(fiducialRA)
            fiducialDecs.append(fiducialDec)
            c= SkyCoord(ra= fiducialRA*u.radian, dec= fiducialDec*u.radian)
            fiducialGalacticLat.append(c.transform_to(Galactic).b.radian)
            if disc:
                contigIDs.append(findContigFOVs(fiducialRA, fiducialDec, ID, FOV_radius,
                                                simdata, disc= True, nside= nside))
        if not allIDs:
            check= copy.deepcopy(coaddBundle[focusDither])
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
            plt.show()

    output= {}
    output['regionPixels']= np.array(goodPixelNums)
    output['goodFiducialIDs']= np.array(goodCenterIDs)
    output['diffMeanMedian']= np.array(diffMeanMedian)
    output['scatterInDepth']= np.array(scatterInDepth)
    output['fiducialRA']= np.array(fiducialRAs)
    output['fiducialDec']= np.array(fiducialDecs)
    output['fiducialGalacticLat']= np.array(fiducialGalacticLat)
    if disc: output['contigIDs']= contigIDs
    
    return output

