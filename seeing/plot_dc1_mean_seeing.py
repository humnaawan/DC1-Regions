# Goal here is to produce maps for specified quantities baseline2018a and
# pontus_2002 (wider footprint) after 1yr or 10-yrs.
#
# Humna Awan: humna.awan@rutgers.edu
#
########################################################################
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.stackers as mafStackers
import time
import healpy as hp
import numpy as np
########################################################################
outdir = '/global/cscratch1/sd/awan/seeing_minion_2016_new_dithers'
db_path = '/global/cscratch1/sd/awan/dbs_old_unzipped/minion_1016_sqlite_new_dithers.db'
nside = 256 #1024
quantity = 'FWHMeff'
db_tag = 'minion1016_new_diths'
########################################################################
# set up
resultsDb = db.ResultsDb(outDir=outdir)

avg_seeing_bundle = {}
bands = ['r'] #['u', 'g', 'r', 'i', 'z', 'y']

year_tag = 'y10'
########################################################################
print('\n## Running for %s\n'%(year_tag))

for dith in ['nodith', 'random_pernight']:
    startTime = time.time()
    # OpSim database
    opsdb = db.OpsimDatabase(db_path)
    
    # WFD only
    propIds, propTags = opsdb.fetchPropInfo()
    wfdWhere = opsdb.createSQLWhere('WFD', propTags)

    # implement per night dithers
    stackerList = [mafStackers.RandomDitherPerNightStacker(degrees=opsdb.raDecInDeg, randomSeed=1000)]
    if dith == 'nodith':
        lonCol, latCol = 'fieldRA', 'fieldDec'
    elif dith == 'random_pernight':
        lonCol, latCol = 'randomDitherPerNightRa', 'randomDitherPerNightDec'
    else:
        raise ValueError('Invalid dith: %s'%dith)
    
    slicer = slicers.HealpixSlicer(lonCol=lonCol, latCol=latCol,
                               latLonDeg=opsdb.raDecInDeg, nside=nside, useCache=False)
    meanMetric= metrics.MeanMetric(col=quantity)   # for avg quantity per HEALpix pixel

    avg_seeing_bundle = {}

    for band in bands:
        sqlconstraint = '%s and filter=="%s"'%(wfdWhere, band)

        avg_seeing_bundle[band] = metricBundles.MetricBundle(meanMetric, slicer, sqlconstraint,
                                                                       stackerList=stackerList)
    grp = metricBundles.MetricBundleGroup(avg_seeing_bundle, opsdb, outDir=outdir,
                                                 resultsDb=resultsDb,saveEarly= False)
    grp.runAll()
    print('\n## Time taken: %.2f min\n'%((time.time()-startTime)/60.))

    # plot skymaps
    for i, band in enumerate(bands):
        plt.clf()
        if (i==0): # define the color range
            inSurveyIndex = np.where(avg_seeing_bundle[band].metricValues.mask == False)[0]
            median = np.median(avg_seeing_bundle[band].metricValues.data[inSurveyIndex])
            stddev = np.std(avg_seeing_bundle[band].metricValues.data[inSurveyIndex])

            nTicks = 5        
            colorMin = 0.8
            colorMax = 1.4

            increment = (colorMax-colorMin)/float(nTicks)
            ticks = np.arange(colorMin+increment, colorMax, increment)

        hp.mollview(avg_seeing_bundle[band].metricValues.filled(avg_seeing_bundle[band].slicer.badval),
                    flip='astro', rot=(0,0,0), min=colorMin, max=colorMax,
                    title= 'WFD', cbar=False)
        hp.graticule(dpar=20, dmer=20, verbose=False)

        fig = plt.gcf()
        im = plt.gca().get_images()[0]
        cbaxes = fig.add_axes([0.25, 0.18, 0.5, 0.03]) # [left, bottom, width, height]
        cb = plt.colorbar(im, orientation='horizontal',
                              ticks=ticks, format='%.2f', cax=cbaxes)
        cb.set_label('%s-band: %s after %s'%(band, quantity, year_tag), fontsize=12)
        cb.ax.tick_params(labelsize=10)

        fig.set_size_inches(5, 5)
        filename = '%s_%s_%s_nside%s_%sband_%s.png'%(quantity, db_tag, year_tag, nside, band, dith)
        plt.savefig('%s/%s'%(outdir, filename), format='png',  bbox_inches='tight')
        print('## Saved %s in outdir.\n'%filename)


    for band in avg_seeing_bundle:
        outfile = 'seeing_%s_%s.npz'%(band, dith)
        avg_seeing_bundle[band].slicer.writeData(outfilename='%s/%s'%(outdir, outfile),
                                                 metricValues=avg_seeing_bundle[band].metricValues,
                                                 metricName=avg_seeing_bundle[band].metric.name,
                                                 simDataName=avg_seeing_bundle[band].runName,
                                                 constraint=avg_seeing_bundle[band].constraint,
                                                 metadata=avg_seeing_bundle[band].metadata,
                                                 displayDict=avg_seeing_bundle[band].displayDict,
                                                 plotDict=avg_seeing_bundle[band].plotDict)