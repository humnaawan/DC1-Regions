import numpy as np
import lsst.sims.maf.slicers as slicers
from lsst.obs.lsstSim import LsstSimMapper
from lsst.sims.utils import ObservationMetaData
from lsst.sims.coordUtils import chipNameFromRaDec
from DC1code.intermediates import printProgress, getSimData, findRegionPixels, findRegionFOVs
import pickle
import time
import os
import datetime
from numpy.lib.recfunctions import append_fields

__all__= ['findDC1Chips']

def findDC1Chips(dbpath, newAfterburner, fiducialDither, fiducialID,
                 addRotDither= False, nside= 512,
                 filterBand= 'r', disc= True, FOV_radius= 0.0305,
                 saveData= True, outputPath= None):
    """

    Find the chips in each visit that come into play in the region based on
    the given fieldID.

    Required Parameters
    -------------------
    * dbpath: str: path to the OpSim database.
    * newAfterburner: bool: set to True if the opsim database is using the 
                           new afterbuner. If not, False.
    * fiducialDither: str: observing strategy to base the chips on.
                           Possibilities: 'NoDither', 'SequentialHexDitherPerNight'
                           If newAfterburner= True, could also be 'RandomDitherFieldPerVisit'
    * fiducialID: int: fieldID for the FOV on which to base the region.

    Optional Parameters
    -------------------
    * addRotDither: bool: set to True to add (random) rotational dithers; will only work
                          if newAfterburner= True. Default: False
      ***** NOTE: assume here that the parallactic angle does not change with rotational
                  dithers; not entirely accurate.
    * filterBand: str: filter to consider. Default: 'r'
    * disc: bool: set to False if don't want disc-like region; will implement the one
                  in intermediates.enclosingRegion. Default: True
    * FOV_radius: float: radius of the FOV (radians). Default: 0.0305
    * saveData: bool: set to False if don't want to save the output data (includes
                      obsIDsList, expDatesList, fIDsList, chipNamesList; saved as a
                      pickle). Default: True
    * outputPath: str: path to the output directory where the data should be saved.
                       If None given, will save in the current working directory.
                       Default: None

    Outputs
    -------
    obsIDsList, expDatesList, fIDsList, chipNamesList, pixNumsList
    Default parameters save the output data as a pickle. See above.

    """
    if newAfterburner:
        if (fiducialDither!= 'NoDither') and (fiducialDither!= 'SequentialHexDitherPerNight') and (fiducialDither!= 'RandomDitherFieldPerVisit'):
            printProgress('Problem: fiducialDither invalid. Must be one of NoDither, SequentialHexDitherPerNight, RandomDitherFieldPerVisit.', highlight= True)
            return
    else:
        if (fiducialDither!= 'NoDither') and (fiducialDither!= 'SequentialHexDitherPerNight'):
            printProgress('Problem: fiducialDither invalid. Must be one of NoDither or SequentialHexDitherPerNight.', highlight= True)
            return

    if addRotDither and not newAfterburner:
        printProgress('Problem: need the new afterburner output to add rotational dithers.', highlight= True)
        return
        
    #nside= 512   # needed for HEALPix pixels to be smaller than the chips themselves
    extraCols= ['expDate', 'obsHistID']
    if addRotDither: extraCols+=['ditheredRotTelPos', 'rotTelPos']
    
    printProgress(datetime.datetime.isoformat(datetime.datetime.now())+'\nGetting simData ... ', highlight= True)
    simdata= getSimData(dbpath, filterBand,
                        extraCols= extraCols,
                        newAfterburner= newAfterburner)  # need the two extra columns to tag different visits.

    # set up
    pointingRACol, pointingDecCol, rotSkyCol= [], [], []
    if (fiducialDither=='NoDither'):
        pointingRACol= 'fieldRA'
        pointingDecCol= 'fieldDec'
    else:
        if newAfterburner:
            if (fiducialDither=='SequentialHexDitherPerNight'):
                pointingRACol= 'hexDitherPerNightRA'
                pointingDecCol= 'hexDitherPerNightDec'
            else:
                # above check should've ensured that fiducialDither= 'RandomDitherFieldPerVisit'
                pointingRACol= 'randomDitherFieldPerVisitRA'
                pointingDecCol= 'randomDitherFieldPerVisitDec'
        else:
            # fiducialDither should be 'SequentialHexDitherPerNight'
            pointingRACol= 'ditheredRA'
            pointingDecCol= 'ditheredDec'

    if addRotDither:
        ### NOTE: assuming that the parallactic angle does not change with rotational dithers,
        ###       which is not entirely accurate.
        parallacticAngle= simdata['rotTelPos']-simdata['rotSkyPos'][:]
        ditheredRotSkyPos= simdata['ditheredRotTelPos'][:]-parallacticAngle
        simdata = append_fields(simdata, 'ditheredRotSkyPos', ditheredRotSkyPos,
                                dtypes= ditheredRotSkyPos.dtype, usemask=False, asrecarray=False)
        rotSkyCol= 'ditheredRotSkyPos'
    else:
        rotSkyCol= 'rotSkyPos'
        
    printProgress('Finding region pixels ... ', highlight= True)
    centralRA, centralDec, regionPixels= findRegionPixels(fiducialID, simdata,
                                                          nside,
                                                          disc, FOV_radius)
    totPixels= len(regionPixels)
    printProgress('Total number of pixels in the region: %d'%(totPixels))
    
    # set up the slicer
    hpSlicer= slicers.HealpixSlicer(nside= nside, 
                                    lonCol= pointingRACol, latCol= pointingDecCol,
                                    rotSkyPosColName= rotSkyCol)
    hpSlicer.setupSlicer(simdata)    # slice data: know which pixels are observed in which visit
    
    camera = LsstSimMapper().camera
    chipNames, obsIDs, expDates, fIDs, pixNums= [], [], [], [], []

    # find the possible fIDs that could come into the region.
    #valid_fields= findRegionFOVs(regionPixels, 'NoDither', simdata, nside= nside)
    #indObsInPixel= []   # all possibilites
    #for ID in valid_fields:   # run over all IDs that could come into the region
    #    ind= np.where(simdata['fieldID']==ID)[0]
    #    indObsInPixel+= list(ind)
        
    prevPercent= 0.
    startTime= time.time()
    for p, pixel in enumerate(regionPixels):  # run over all the pixels in the region
        pixRA, pixDec= hpSlicer._pix2radec(pixel)    # radians returned
        indObsInPixel = hpSlicer._sliceSimData(pixel)   # indices in simData for when an observation
                                                        # happened in this pixel
        
        for index in indObsInPixel['idxs']:
            # get data from simdata for identifying each visit
            expDate= simdata[index]['expDate']
            obsID= simdata[index]['obsHistID']
            fID= simdata[index]['fieldID']
            
            # for chip finding
            pointingRA= simdata[index][pointingRACol] # radians
            pointingDec= simdata[index][pointingDecCol] # radians
            rotSkyPos= simdata[index][rotSkyCol] # radians
            expMJD= simdata[index]['expMJD']
            
            # set up for the finding the chips
            obs = ObservationMetaData(pointingRA= np.degrees(pointingRA), pointingDec= np.degrees(pointingDec),
                                      rotSkyPos= np.degrees(rotSkyPos), mjd= expMJD)
            chipsInVisit= chipNameFromRaDec(np.degrees(pixRA), np.degrees(pixDec),
                                            camera=camera, obs_metadata=obs)
            if chipsInVisit is not None:
                 if not (chipsInVisit.__contains__('A')) and not (chipsInVisit.__contains__('B')): # no wavefront sensors
                     obsIDs.append(obsID)
                     expDates.append(expDate)
                     chipNames.append(chipsInVisit)
                     fIDs.append(fID)
                     pixNums.append(pixel)

        percentDone= 100.*(p+1)/totPixels
        delPercent= percentDone-prevPercent
        if (delPercent>5):
            printProgress('%f%% pixels done\nTime passed (min): %f'%(percentDone, (time.time()-startTime)/60.), newLine= False)
            prevPercent= percentDone
        #if (percentDone>1.):
        #    break
    obsIDs, expDates, fIDs, chipNames, pixNums= np.array(obsIDs), np.array(expDates), np.array(fIDs), np.array(chipNames), np.array(pixNums)
    
    printProgress('Unique obsHistIDs: %d \n## Unique expDates: %d \n## Unique chipNames: %d \n'%(len(np.unique(obsIDs)),
                                                                                                 len(np.unique(expDates)),
                                                                                                 len(np.unique(chipNames))), newLine= True)

    #  get rid of repeated entries; consolidate the data from unique observations.
    printProgress('Consolidating the data ... ', highlight= True)
    obsIDsList, expDatesList, fIDsList, chipNamesList, pixNumsList= [], [], [], [], []
    for obs in np.unique(obsIDs):
        obsIDsList.append(obs)
        ind= np.where(obsIDs==obs)[0]
        expDatesList.append(np.unique(expDates[ind]))
        fIDsList.append(np.unique(fIDs[ind]))
        chipNamesList.append(np.unique(chipNames[ind]))
        pixNumsList.append(np.unique(pixNums[ind]))
        
    # see how many chips are added by any given visit
    numChips= []
    for i in range(len(obsIDsList)):
        numChips.append(len(chipNamesList[i]))
    printProgress('Max number of chips added by any given visit: %d \n## Min number: %d \n## Total number of chips (across all visits to be simulated): %d'
                  %(max(numChips), min(numChips), sum(numChips)), highlight= True)
    # save the data?
    if saveData:
        printProgress('Saving the data ... ', highlight= True)
        dataToSave = {'obsHistID': obsIDsList, 'expDate': expDatesList,
                      'fIDs': fIDsList, 'chipNames': chipNamesList,
                      'pixNum': pixNumsList}
        if outputPath is not None:
            currentDir= os.getcwd()
            os.chdir(outputPath)

        shapeTag= ''
        if disc: shapeTag= 'disc'
        else: shapeTag= 'nonDisc'

        burnerTag= ''
        if newAfterburner: burnerTag = 'newAfterburnerOutput_'
        rotDitherTag= ''
        if addRotDither: rotDitherTag= 'randomRotDithered_'
        filename= '%s_chipPerVisitData_%sfID%d_%s_%s%sRegion_nside%s_%dNonWFChipsToSimulate'%(str(datetime.date.isoformat(datetime.date.today())),
                                                                                            burnerTag, fiducialID, fiducialDither,
                                                                                            rotDitherTag, shapeTag, nside, sum(numChips))
        with open(filename+'.pickle', 'wb') as handle:
            pickle.dump(dataToSave, handle, protocol=pickle.HIGHEST_PROTOCOL)

        printProgress('Saved the data as %s.pickle'%(filename), highlight= True)
        if outputPath is not None:
            os.chdir(currentDir)

    return [obsIDsList, expDatesList, fIDsList, chipNamesList, pixNumsList]
