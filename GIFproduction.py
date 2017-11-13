from intermediates import *
import numpy as np
import matplotlib.pyplot as plt
from lsst.sims.coordUtils import raDecFromPixelCoords
from numpy.lib.recfunctions import append_fields
from lsst.obs.lsstSim import LsstSimMapper
from lsst.sims.utils import ObservationMetaData, altAzPaFromRaDec
from lsst.sims.coordUtils import getCornerRaDec
from readPython2Pickle import pickleRead # needed since working in Python 3.
import healpy as hp
import os
import imageio
import matplotlib.patches as mpatches

def obsMetaDataDict(obsHistIDs, simdata, transDithered, rotDithered):
    """
    Function to get the metadata for each observation corresponding to input obsHistIDs.

    Required Inputs
    ---------------
    * obsHistIDs: list/array: visit IDs.
    * simdata: OpSim database file.
    * transDithered: bool: True if want metadata for pointings with translational dithers (randomFieldPerVisit rn).
                           Else will have undithered pointings.
    * rotDithered: bool: True if want metadata for pointings with rotational dithers (random dithers).
                         Else will have metadata for undithered pointings.

    Returns
    -------
    * Dictionary of ObservationMetaData objects for each objID; keys= objIDs.

    """
    if transDithered:
        pointingRACol= 'randomDitherFieldPerVisitRA'
        pointingDecCol= 'randomDitherFieldPerVisitDec'
    else:
        pointingRACol= 'fieldRA'
        pointingDecCol= 'fieldDec'
    if rotDithered:
        rotSkyCol= 'ditheredRotSkyPos'
    else:
        rotSkyCol= 'rotSkyPos'

    obsMetaData= {}
    for hid in np.unique(obsHistIDs):
        inds= np.where(simdata['obsHistID']==hid)[0] # all the entries corresponding to the objID.

        if (len(inds)>1):
            raise TypeError("More than one visit found for an objID.")

        for index in inds:
            pointingRA= simdata[index][pointingRACol] # radians
            pointingDec= simdata[index][pointingDecCol] # radians
            rotSkyPos= simdata[index][rotSkyCol] # radians
            expMJD= simdata[index]['expMJD']
            obsMetaData[hid] = ObservationMetaData(pointingRA= np.degrees(pointingRA), pointingDec= np.degrees(pointingDec),
                                              rotSkyPos= np.degrees(rotSkyPos), mjd= expMJD)
    return obsMetaData

def movie(filepath, nside, simdata, fIDmovie, histIDmovie, dataTag,
          transDithered, rotDithered, outputDir, testRun= False,
          obsHistIndMin= 50, obsHistIndMax= 150):
    """
    Function to create a timeseries out of pngs. Doesnt return anything; saves GIFs (temporarily the pngs).

    Required Inputs
    ---------------
    * filepath: str: path to the file containing the pickled dictionary with chipNames, obsHistIDs, etc.
    * nside: int: HEALpix resolution parameter used to produce the file in filepath.
    * simdata: OpSim database.
    * fIDmovie: bool: True if want the GIF to be w.r.t fID, i.e. all the visits corresponding to a given FOV
                      will be plotted sequentially.
    * histIDmovie: bool: True if want the GIF to be w.r.t objsID, i.e. visits for a set of histIDs will be plotted.
    * dataTag: str: tag to identify the saved GIF; will also appear in the title.
    * transDithered: bool: True if want output with translational dithers (randomFieldPerVisit rn).
                           Else will have undithered pointings.
    * rotDithered: bool: True if want output for pointings with rotational dithers (random dithers).
                         Else will have undithered pointings.
    * outputDir: str: path where the GIF will be saved.

    Optional Inputs:
    ---------------
     * testRun: bool: set to run if just want to check the method runs. If histIDmovie= False, will only
                      only consider 2-3 visits; if fIDmovie= True, will only run for one fID.
                     Default: False
    * obsHistIndMin, obsHistIndMax: int, int: index of the obsHistID to consider. If histIDmovie= True,
                                        obsHistIndMin:obsHistIndMax will be plotted.
                                        Default: 50, 150 (arbitrarily chosen).

    """
    # non-science chips
    wavefront= ['R:4,0 S:0,2,A', 'R:4,0 S:0,2,B','R:0,4 S:2,0,A', 'R:0,4 S:2,0,B',
                'R:0,0 S:2,2,A', 'R:0,0 S:2,2,B','R:4,4 S:0,0,A', 'R:4,4 S:0,0,B']
    guiders= ['R:0,0 S:1,2', 'R:0,0 S:2,1', 'R:0,4 S:1,0',
                'R:0,4 S:2,1', 'R:4,0 S:0,1', 'R:4,0 S:1,2', 'R:4,4 S:0,1', 'R:4,4 S:1,0']

    def pixNum2RaDec(nside, pix):
        # convert from HEALpix pixels to ra, dec.
        theta, phi = hp.pix2ang(pix, nside)
        return phi, np.pi/2.0-theta

    # for DC1, fiducial fieldID= 1447. Get the region pixels.
    eh, eh, DC1Pixels= findRegionPixels(1447, simdata, nside= nside, disc= False, FOV_radius= 0.0305)
    DC1pixRA, DC1pixDec= pixNum2RaDec(DC1Pixels, nside)

    # read in the data
    with open(filepath, 'rb') as handle:
        inputData= pickleRead(handle)
    # check on max index
    if (obsHistIndMax>len(inputData['obsHistID'])):
        raise ValueError("obsHistIndMax must be less than total number of histIDs.")

    # get the observation metadata for the various visits.
    obsMetaData= obsMetaDataDict(inputData['obsHistID'], simdata, transDithered, rotDithered)
    camera= LsstSimMapper().camera

    # plot helpers
    factor= .01
    xmin, xmax= 1.54, 1.72
    ymin, ymax= -0.56, -0.42
    kargs= {'fps': 1}
    fontsize= 14

    # set up output directory.
    os.chdir(outputDir)
    outDir= 'DC1_gifs'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    os.chdir(outDir)

    # set up for circular FOV footprints.
    stepsize = np.pi/50.
    theta = np.arange(0, np.pi*2.+stepsize, stepsize)
    radius= 0.0305
    delx= radius*np.cos(theta)
    dely= radius*np.sin(theta)

    # find pixels related to each FOV for plotting purposes.
    FOVpixs= {}
    for fid in np.unique(inputData['fIDs']):
        centralRA, centralDec, regionPixels= findFOVPixels(fid, simdata, nside= nside,
                                                            FOV_radius= radius)
        FOVpixs[fid]= regionPixels

   # function to get chips boundaries.
    def getXYchips(chip, camera, obsData):
        x, y= np.zeros(5), np.zeros(5)
        [(x[0], y[0]), (x[1], y[1]), (x[3], y[3]), (x[2], y[2])]= getCornerRaDec(chip, camera, obsData)
        # returns [(xmin, ymin), (xmin, ymax), (xmax, ymin), (xmax, ymax)] so need to order things
        # for one direction.
        x[4], y[4]= x[0], y[0]  # loop back for closed chip boundaries.
        x, y= np.deg2rad(x), np.deg2rad(y)  # working in radians.
        return x,y

    if histIDmovie:
        filenames= []
        ###########################################################################
        # loop over histIDs
        if testRun: obsHistIndMax= obsHistIndMin+2
        for index, hid in enumerate(inputData['obsHistID'][obsHistIndMin:obsHistIndMax]):
            fig, axes= plt.subplots(1,1)
            # plot the full DC1 region
            axes.plot(DC1pixRA, DC1pixDec, 'o', color= 'g', alpha= 0.03)
            i= np.where(inputData['obsHistID']==hid)[0][0]

            # plot the undithered FOV pixels
            fid= inputData['fIDs'][i][0]
            pixRA, pixDec= pixNum2RaDec(FOVpixs[fid], nside)
            axes.plot(pixRA, pixDec, '.', color= 'k', alpha= 0.05)

            # plot the FOV boundary/center.
            pointingRA, pointingDec= np.radians(obsMetaData[hid].pointingRA), np.radians(obsMetaData[hid].pointingDec)
            axes.plot(pointingRA, pointingDec, 'x', color= 'k')
            axes.plot(pointingRA+delx/np.cos(pointingDec), pointingDec+dely, color= 'k', label= 'FOV')

            # get/plot the chips
            for chip in inputData['chipNames'][i]:
                x, y= getXYchips(chip, camera, obsMetaData[hid])  # chip corners
                axes.fill(x,y, 'b', alpha=0.2, edgecolor='r')

            # plot the wavefront sensors
            wavefrontColors= ['r', 'c', 'b', 'm', 'r', 'b', 'c', 'm']
            x0, y0, x1, y1= 0., 0., 0., 0
            for ith, chip in enumerate(wavefront):
                x, y= getXYchips(chip, camera,obsMetaData[hid])
                if (ith==0):
                    x0, y0= x[3], y[3]
                if (ith==2):
                    x1, y1= x[3], y[3]
                axes.fill(x,y, wavefrontColors[ith], alpha=1., edgecolor=wavefrontColors[ith])
            delx_arrow= pointingRA-x0
            dely_arrow= pointingDec-y0
            axes.axes.arrow(x0, y0, delx_arrow, dely_arrow, color= 'k',
                            width= 0.001/2., head_width= 0.003)
            axes.axes.arrow(x0+delx_arrow, y0+dely_arrow, delx_arrow, dely_arrow,
                            color= 'k', width= 0.001/2., head_width= 0.0001/2.)
            #axes.axes.arrow(x0, y0, delx_arrow1, dely_arrow1, color= 'k', width= 0.001/2.)
            #axes.axes.arrow(x0+delx_arrow1, y0+dely_arrow1, delx_arrow1*2, dely_arrow1*2,
            #                color= 'k', width= 0.001/2., head_width= 0.0002)
            #axes.axes.arrow(x0+3*delx_arrow1, y0+3*dely_arrow1, delx_arrow1, dely_arrow1,
            #                color= 'k', width= 0.001/2., head_width= 0.001/2.)

            # plot the guidors
            for chip in guiders:
                x, y= getXYchips(chip, camera,obsMetaData[hid])
                axes.fill(x,y, 'y', alpha=1.0, edgecolor='y')

            # set up for title: include histIS, filter, rot angle
            inds= np.where(simdata['obsHistID']==hid)[0]
            filterObs= simdata['filter'][inds]
            rotTelPos= simdata['rotTelPos'][inds]
            dithrotTelPos= simdata['ditheredRotTelPos'][inds]
            rotSkyPos= simdata['rotSkyPos'][inds]
            dithrotSkyPos= simdata['ditheredRotSkyPos'][inds]

            title= '%s\nhID: %s ; fid: %s ; expDate: %s ; filter: %s \n\n'%(dataTag,
                                                                        inputData['obsHistID'][i] ,inputData['fIDs'][i],
                                                                        inputData['expDate'][i],
                                                                        filterObs)
            title+= 'rotSkyPos: %s, dithRotSkyPos: %s \n'%(rotSkyPos, dithrotSkyPos)
            title+= 'rotTelPos: %s, dithRotTelPos: %s \n'%(rotTelPos, dithrotTelPos)
            title+= 'parallacticAngle= rotTelPos-rotSkyPos : %s \n'%(rotTelPos-rotSkyPos)
            #title+= 'parallacticAngle from utils: %s \n'%np.radians(altAzPaFromRaDec(obsMetaData[hid].pointingRA,
            #                                                    obsMetaData[hid].pointingDec,
            #                                                    obsMetaData[hid])[2])
            title+= '$\Delta$rotSkyPos= dithrotSkyPos-rotSkyPos : %s \n'%(dithrotSkyPos-rotSkyPos)
            title+= '$\Delta$rotTelPos= dithRotTelPos-rotTelPos : %s \n'%(dithrotTelPos-rotTelPos)

            axes.set_title(title,  fontsize= fontsize, )

            # axes labels
            axes.set_xlabel('RA (rad)', fontsize= fontsize)
            axes.set_ylabel('Dec (rad)', fontsize= fontsize)
            #axes.axis('equal')

            # axes limits
            axes.set_xlim([xmin-factor, xmax+factor])
            axes.set_ylim([ymin-factor, ymax+factor])
            axes.tick_params(axis='x', labelsize=fontsize)
            axes.tick_params(axis='y', labelsize=fontsize)

            # set up legend
            blue = mpatches.Patch(facecolor='b', edgecolor= 'red', label= 'DC1 Chips', alpha= 0.2)
            yellow = mpatches.Patch(facecolor='y' , label= 'guiders', alpha= 1.0)
            handles=  [blue, yellow]
            for c in wavefrontColors[0:4]:
                handles.append(mpatches.Patch(color=c, label= 'wavefront'))
            axes.legend(handles= handles, loc= 'upper left')

            fig.set_size_inches(10,10)
            # save each png.
            filename= '%sth.png'%index
            plt.savefig(filename, format= 'png', bbox_inches = 'tight')
            filenames.append(filename)
            plt.close()

        # compile saved pngs; delete them after.
        images = []
        for filename in filenames:
            images.append(imageio.imread(filename))
            os.remove(filename)  # delete png.
        if testRun:
            movName= 'hIDs_%s-%sIndex_%s_testRun.gif'%(obsHistIndMin, obsHistIndMax, dataTag)
        else:
            movName= 'hIDs_%s-%sIndex_%s.gif'%(obsHistIndMin, obsHistIndMax, dataTag)
        imageio.mimsave(movName, images, **kargs)
        print('Saved ', movName)

    if fIDmovie:
        ###########################################################################
        fIDsToConsider= np.unique(inputData['fIDs'])
        if testRun: fIDsToConsider= [fIDsToConsider[0]]
        for fid in fIDsToConsider:
            filenames= []
            inds= np.where(inputData['fIDs']==fid)[0]

            for index, fid_visitInd in enumerate(inds):  # run over all the visits related to the fID.
                fig, axes= plt.subplots(1,1)
                # plot the full DC1 region
                axes.plot(DC1pixRA, DC1pixDec, 'o', color= 'g', alpha= 0.03)
                # plot the undithered FOV pixels
                pixRA, pixDec= pixNum2RaDec(FOVpixs[fid], nside)
                axes.plot(pixRA, pixDec, '.', color= 'k', alpha= 0.05)

                hid= inputData['obsHistID'][fid_visitInd]
                # plot the FOV center/boundary.
                pointingRA, pointingDec= np.radians(obsMetaData[hid].pointingRA), np.radians(obsMetaData[hid].pointingDec)
                axes.plot(pointingRA, pointingDec, 'x', color= 'k')
                axes.plot(pointingRA+delx/np.cos(pointingDec), pointingDec+dely, color= 'k', label= 'FOV')

                # get/plot the chips
                for chip in inputData['chipNames'][fid_visitInd]:
                    x, y= getXYchips(chip, camera, obsMetaData[hid])
                    axes.fill(x,y, 'b', alpha=0.2, edgecolor='r')

                # plot the wavefront sensors
                wavefrontColors= ['r', 'c', 'b', 'm', 'r', 'b', 'c', 'm']
                x0, y0, x1, y1= None, None, None, None
                for ith, chip in enumerate(wavefront):
                    x, y= getXYchips(chip, camera,obsMetaData[hid])
                    if (ith==0):
                        x0, y0= x[0], y[0]
                    if (ith==len(wavefront)):
                        x1, y1= x[0], y[0]
                    axes.fill(x,y, wavefrontColors[ith], alpha=1., edgecolor=wavefrontColors[ith])
                delx_arrow= pointingRA-x0
                dely_arrow= pointingDec-y0
                axes.axes.arrow(x0, y0, delx_arrow, dely_arrow, color= 'k',
                                width= 0.001/2., head_width= 0.003)
                axes.axes.arrow(x0+delx_arrow, y0+dely_arrow, delx_arrow, dely_arrow,
                                color= 'k', width= 0.001/2., head_width= 0.0001/2.)

                # plot the guiders
                for chip in guiders:
                    x, y= getXYchips(chip, camera,obsMetaData[hid])
                    axes.fill(x,y, 'y', alpha=1., edgecolor='y')

                # set up for title
                inds= np.where(simdata['obsHistID']==hid)[0]
                filterObs= simdata['filter'][inds]
                rotTelPos= simdata['rotTelPos'][inds]
                dithrotTelPos= simdata['ditheredRotTelPos'][inds]
                rotSkyPos= simdata['rotSkyPos'][inds]
                dithrotSkyPos= simdata['ditheredRotSkyPos'][inds]
                title= '%s\nhID: %s ; fid: %s ; expDate: %s ; filter: %s \n\n'%(dataTag,
                                                                inputData['obsHistID'][fid_visitInd] ,
                                                                inputData['fIDs'][i],
                                                                inputData['expDate'][fid_visitInd],
                                                                filterObs)
                title+= 'rotSkyPos: %s, dithRotSkyPos: %s \n'%(rotSkyPos, dithrotSkyPos)
                title+= 'rotTelPos: %s, dithRotTelPos: %s \n'%(rotTelPos, dithrotTelPos)
                title+= 'parallacticAngle= rotTelPos-rotSkyPos : %s \n'%(rotTelPos-rotSkyPos)
                title+= '$\Delta$rotSkyPos= dithrotSkyPos-rotSkyPos : %s \n'%(dithrotSkyPos-rotSkyPos)
                title+= '$\Delta$rotTelPos= dithRotTelPos-rotTelPos : %s \n'%(dithrotTelPos-rotTelPos)
                axes.set_title(title,  fontsize= fontsize, )

                # axes labels
                axes.set_xlabel('RA (rad)', fontsize= fontsize)
                axes.set_ylabel('Dec (rad)', fontsize= fontsize)
                axes.tick_params(axis='x', labelsize=fontsize)
                axes.tick_params(axis='y', labelsize=fontsize)

                # legend
                blue = mpatches.Patch(facecolor='b', edgecolor= 'red', label= 'DC1 Chips', alpha= 0.2)
                yellow = mpatches.Patch(facecolor='y' , label= 'guiders', alpha= 1.0)
                handles=  [blue, yellow]
                for c in wavefrontColors[0:4]:
                    handles.append(mpatches.Patch(color=c, label= 'wavefront'))
                axes.legend(handles= handles, loc= 'upper left')

                # axes lims.
                axes.set_xlim([xmin-factor, xmax+factor])
                axes.set_ylim([ymin-factor, ymax+factor])
                fig.set_size_inches(10,10)
                # save pngs
                filename= '%sth.png'%index
                plt.savefig(filename, format= 'png', bbox_inches = 'tight')
                filenames.append(filename)
                plt.close()
            # compile saved pngs; delete them after.
            images = []
            for filename in filenames:
                images.append(imageio.imread(filename))
                os.remove(filename) # remove.
            if testRun:
                movName= 'fID%s_%s_testRun.gif'%(fid, dataTag)
            else:
                movName= 'fID%s_%s.gif'%(fid, dataTag)
            imageio.mimsave(movName, images, **kargs)

            print('Saved ', movName)

def GIFproduction(dbpath, filepath, nside, fIDmovie, histIDmovie, dataTag,
                    transDithered, rotDithered, outputDir, testRun= False,
                    obsHistIndMin= 50, obsHistIndMax= 150):
    """
    Method to produce GIFs for either the timeseries wrt objHistIDs or the visit for each fID.

    Required Inputs:
    ---------------
    * dbpath: str: path to the OpSim database. MUST contain dithered rotational and translational
                   dither columns.
   * filepath: str: path to the file containing the pickled dictionary with chipNames, obsHistIDs, etc.
   * nside: int: HEALpix resolution parameter used to produce the file in filepath.
   * fIDmovie: bool: True if want the GIF to be w.r.t fID, i.e. all the visits corresponding to a given FOV
                     will be plotted sequentially.
   * histIDmovie: bool: True if want the GIF to be w.r.t objsID, i.e. visits for a set of histIDs will be plotted.
   * dataTag: str: tag to identify the saved GIF; will also appear in the title.
   * transDithered: bool: True if want output with translational dithers (randomFieldPerVisit rn).
                          Else will have undithered pointings.
   * rotDithered: bool: True if want output for pointings with rotational dithers (random dithers).
                        Else will have undithered pointings.
   * outputDir: str: path where the GIF will be saved.

   Optional Inputs:
   ---------------
    * testRun: bool: set to run if just want to check the method runs. If histIDmovie= False, will only
                     only consider 2-3 visits; if fIDmovie= True, will only run for one fID.
                    Default: False
    * obsHistIndMin, obsHistIndMax: int, int: index of the obsHistID to consider. If histIDmovie= True,
                                                obsHistIndMin:obsHistIndMax will be plotted.
                                                Default: 50, 150 (arbitrarily chosen).
   """
    if ((obsHistIndMin<0) or (obsHistIndMax<0)):
        raise ValueError("Must have obsHistIndMin, obsHistIndMax > 0. They are indices to the histID array.")
    # read in the OpSim data.
    extraCols= ['expDate', 'obsHistID', 'ditheredRotTelPos', 'rotTelPos', 'filter']
    simdata= getSimData(dbpath, 'r', newAfterburner= True, extraCols= extraCols)

    parallacticAngle= simdata['rotTelPos']-simdata['rotSkyPos'][:]
    ditheredRotSkyPos= simdata['ditheredRotTelPos'][:]-parallacticAngle
    simdata = append_fields(simdata, 'ditheredRotSkyPos', ditheredRotSkyPos,
                    dtypes= ditheredRotSkyPos.dtype, usemask=False, asrecarray=False)

    # call the function to produce the gifs.
    movie(filepath, nside, simdata, fIDmovie, histIDmovie, dataTag,
            transDithered, rotDithered, outputDir, testRun= testRun,
            obsHistIndMin= obsHistIndMin, obsHistIndMax= obsHistIndMax)

    print('All done.')
    return
