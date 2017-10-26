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
          'findContigFOVs', 'findRegionFOVs', 'findGoodRegions', 'findFOVPixels']

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
    if newLine: print(append + '\n## ' + whatToPrint)
    else: print('## ' + whatToPrint)
    sys.stdout.flush()
    time.sleep(1.0)

def getSurveyHEALPixRADec(coaddBundle):
    """

    Get the RA, Dec (in radians) corresponding to each HEALPix pixel.
    Method returns a dictionary with keys= keys in coaddBundle: pixelNumber

    Required Parameter
    ------------------
    * coaddBundle: dict: dictionary with keys= observing strategy names, pointing to
                         corresponding to a metricBundle object.

    """
    # create dictionaries giving pixelNumbers and their correspondong RA, Dec for all dither strategies.
    # need to worry about each strategy separately since the mask is generally different.
    pixelNum= {}
    for dither in coaddBundle:
        pixelNum[dither]= []
        for pix in range(len(coaddBundle[dither].slicer)):
            if not coaddBundle[dither].metricValues.mask[pix]:   # only consider the unmasked pixels
                pixelNum[dither].append(pix)

    return pixelNum

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
    slicer= slicers.HealpixSlicer(nside= nside)
    slicer.setupSlicer(simdata)

    pixels_in_FOV= {}
    for dither in pixelNum:
        pixels_in_FOV[dither]= {}
        for pixel in pixelNum[dither]:
            indObsInPixel = slicer._sliceSimData(pixel)
            ids = simdata[indObsInPixel['idxs']]['fieldID']   # fieldIDs corresponding to pixel
            for uniqID in np.unique(ids):
                key= uniqID
                if key not in pixels_in_FOV[dither].keys():
                    pixels_in_FOV[dither][key]= []
                pixels_in_FOV[dither][key].append(pixel)
    print('Number of fieldIDs in pixel_in_FOV for %s: %d' %(dither, len(pixels_in_FOV[dither].keys())))
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
        #corners= enclosingPolygon(FOV_radius, fiducialRA, fiducialDec)
        #regionPixels= hp.query_polygon(nside, corners)    # HEALpixel numbers

        contigIDList= findContigFOVs(fiducialRA, fiducialDec, fiducialID, FOV_radius, simdata,
                                     disc= False, nside= nside)

        # flatten the list:
        flatList= []
        for ID in contigIDList:
            if ((type(ID)==np.ndarray) or (type(ID)==list)):
                for i in range(len(ID)):
                    flatList.append(ID[i])
            else:
                flatList.append(ID)
        allPixels= []
        for ID in flatList:
            ind= np.where(simdata[:]['fieldID']== ID)[0]
            fixedRA= simdata[ind[0]]['fieldRA']
            fixedDec= simdata[ind[0]]['fieldDec']
            c = SkyCoord(ra=fixedRA*u.radian, dec= fixedDec*u.radian)
            regionPixels= hp.query_disc(nside= nside, vec=c.cartesian.xyz, radius= FOV_radius, inclusive= True)

            allPixels+= list(regionPixels)
        regionPixels= np.unique(allPixels)
    else:
        fiducialRA, fiducialDec= fixedRA, fixedDec-FOV_radius*np.sqrt(3)/2.
        c = SkyCoord(ra=fiducialRA*u.radian, dec= fiducialDec*u.radian)
        regionPixels= hp.query_disc(nside= nside, vec=c.cartesian.xyz, radius= 2.5*FOV_radius)

    return [fiducialRA, fiducialDec, regionPixels]

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
    slicer= slicers.HealpixSlicer(nside= nside, verbose= False)
    slicer.setupSlicer(simdata)

    def returnIDs(ra, dec):
        p= hp.ang2pix(nside, np.pi/2.0-dec, ra)
        ind = slicer._sliceSimData(p)
        ids = simdata[ind['idxs']]['fieldID']   # fieldIDs corresponding to pixelNum[i]
        uniqID= np.unique(ids)
        #if len(uniqID)>1:
        #    print 'Should really have only one FOV corresponding to the central pixel (fID: %d) but have %s'%(fiducialID, uniqID)
        return list(uniqID)

    contList= []
    contList.append(fiducialID)

    if disc:
        # assume given fiduvicalRA, Dec are the center of the big disc.
        # FOV below
        ra= fiducialRA
        dec= fiducialDec-FOV_radius*np.sqrt(3)/2.
        contList+= returnIDs(ra, dec)

        ra= fiducialRA-FOV_radius*3/2.
        dec= fiducialDec
        contList+= returnIDs(ra, dec)

        ra= fiducialRA+FOV_radius*3/2.
        dec= fiducialDec
        contList+= returnIDs(ra, dec)
    else:
        # FOV below
        ra= fiducialRA
        dec= fiducialDec-2*FOV_radius*np.sqrt(3)/2.
        contList+= returnIDs(ra, dec)

        ra= fiducialRA-FOV_radius*3/2.
        dec= fiducialDec-FOV_radius*np.sqrt(3)/2.
        contList+= returnIDs(ra, dec)

        ra= fiducialRA+FOV_radius*3/2.
        dec= fiducialDec-FOV_radius*np.sqrt(3)/2.
        contList+= returnIDs(ra, dec)

    return np.unique(contList)

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
                    nside= 256, depthDiffThreshold= 0.01, rangeThreshold= 0.3,
                    allIDs= True, IDsToTestWith= [],
                    disc= False,  raRange= [-180,180], decRange= [-70,10]):
    """

    Find good regions (i.e. regions with median depth within specified threshold of survey
    median depth AND with max depth-min depth within a specific threshold).

    Returns: arrays: good pixel numbers, good IDs, difference between median depth in the region and
    median survey depth, abs(max-min) region depth, fiducial RA, fiduciual Dec

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
    * depthDiffThreshold: float: region will be considered good if median depth in the region is within
                                 this threshold of survey median depth. Default: 0.01
    * rangeThreshold: flat: region will be considered if good if range of depth in the region (i.e.
                            difference max-min depth) is within the specified rangeThreshold.
                            Default: 0.3
    * allIDs: bool: set to False to consider only a subset of FOVs; will plot things out. Default: True
    * IDsToTestWith: list: list of fieldIDs to consider. Default: []
    * disc: bool: set to True if want disc-like region; False for rectangular. Default: False
    * raRange: list: constraint on the RA in cartview. Default: [-180,180]
    * decRange: list: constraint on the Dec in cartview. Default: [-70,10]

    """
    # a region is 'good' if abs(median depth in the region -surveyMedianDepth)<depthDiffThreshold and rangeDepth in the region < rangeThreshold
    goodCenterIDs, goodPixelNums, rangeInDepth, diffMedianRegionSurvey= [], [], [], []
    fiducialRAs, fiducialDecs, contigIDs, fiducialGalacticLat=  [], [], [], []

    inSurvey= np.where(coaddBundle[focusDither].metricValues.mask==False)[0]
    surveyMedianDepth= np.median(coaddBundle[focusDither].metricValues.data[inSurvey])
    printProgress('Mean survey depth for %s: %f'% (focusDither, surveyMedianDepth))

    considerIDs= pixels_in_FOV[focusDither].keys()
    if not allIDs: considerIDs= IDsToTestWith

    printProgress('Checking for regions with median depth within %f of the median survey depth AND depth range < %f'
                  % (depthDiffThreshold, rangeThreshold))
    for ID in considerIDs:
        fiducialRA, fiducialDec, regionPixels= findRegionPixels(ID, simdata, nside, disc, FOV_radius)


        notInSurvey= np.where(coaddBundle[focusDither].metricValues.mask[regionPixels]==True)[0]
        if not (len(notInSurvey)>0): # the whole region is in the survey
            typicalDepth= np.median(coaddBundle[focusDither].metricValues.data[regionPixels])
            diffDepth= abs(typicalDepth-surveyMedianDepth)
            depthRange= abs(max(coaddBundle[focusDither].metricValues.data[regionPixels])-min(coaddBundle[focusDither].metricValues.data[regionPixels]))

            if ((diffDepth<depthDiffThreshold) and (depthRange<rangeThreshold)):
                goodCenterIDs.append(ID)
                goodPixelNums.append(regionPixels)
                diffMedianRegionSurvey.append(diffDepth)
                rangeInDepth.append(depthRange)
                fiducialRAs.append(fiducialRA)
                fiducialDecs.append(fiducialDec)
                c= SkyCoord(ra= fiducialRA*u.radian, dec= fiducialDec*u.radian)
                fiducialGalacticLat.append(c.transform_to(Galactic).b.radian)
                contigIDs.append(findContigFOVs(fiducialRA, fiducialDec, ID, FOV_radius,
                                                simdata, disc= disc, nside= nside))
        if not allIDs:
            check= copy.deepcopy(coaddBundle[focusDither])
            check.metricValues.data[:]= 0
            check.metricValues.data[regionPixels]= 1000.
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

    printProgress('Total number of good regions founds: %d'%(len(goodCenterIDs)), highlight= True)
    output= {}
    output['regionPixels']= np.array(goodPixelNums)
    output['goodFiducialIDs']= np.array(goodCenterIDs)
    output['diffMedianRegionSurvey']= np.array(diffMedianRegionSurvey)
    output['rangeInDepth']= np.array(rangeInDepth)
    output['fiducialRA']= np.array(fiducialRAs)
    output['fiducialDec']= np.array(fiducialDecs)
    output['fiducialGalacticLat']= np.array(fiducialGalacticLat)
    output['contigIDs']= contigIDs

    return output


def findFOVPixels(fiducialID, simdata, nside, FOV_radius):
    """

    Find the pixels inside specified (undithered) FOV (identified by a field ID).
    ** Uses the undithered pointing info; looks at circular region
       centered on FOV's pointingRA, pointingDec.

    Required Parameters
    -------------------
    * fiducialID: int: fieldID for the FOV.
    * simdata: np.array: array containing OpSim columns (must have fieldID, fieldRA, fieldDec).
    * nside: int: HEALPix resolution parameter.
    * FOV_radius: float: radius of the FOV in radians.

    """
    ind= np.where(simdata[:]['fieldID']== fiducialID)[0]
    # fieldRA, fieldDec remain fixed for NoDither; dont change with expMJD.
    # use as the 'center' of the enclosing region (disc or rectangle).
    fixedRA= simdata[ind[0]]['fieldRA']
    fixedDec= simdata[ind[0]]['fieldDec']

    fiducialRA, fiducialDec= fixedRA, fixedDec
    c = SkyCoord(ra=fiducialRA*u.radian, dec= fiducialDec*u.radian)
    regionPixels= hp.query_disc(nside= nside, vec=c.cartesian.xyz, radius= FOV_radius)

    return [fiducialRA, fiducialDec, regionPixels]
