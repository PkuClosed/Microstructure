# -*- coding: utf-8 -*-
"""
Created on Wed May 18 16:51:56 2016

@author: cyye
"""
import amico 
import sys
import os

dwinames = sys.argv[1]
masknames = sys.argv[2]
bvecsnames = sys.argv[3]
bvalsnames = sys.argv[4]
sel = int(sys.argv[5])

with open(dwinames) as f:
    allDwiNames = f.readlines()
with open(masknames) as f:
    allMaskNames = f.readlines()
with open(bvecsnames) as f:
    allBvecsNames = f.readlines()
with open(bvalsnames) as f:
    allBvalsNames = f.readlines()

allDwiNames = [x.strip('\n') for x in allDwiNames]
allMaskNames = [x.strip('\n') for x in allMaskNames]
allBvecsNames = [x.strip('\n') for x in allBvecsNames]
allBvalsNames = [x.strip('\n') for x in allBvalsNames]

if sel == 0:
    alg = "NODDI"
if sel == 1:
    alg = "CylinderZeppelinBall"

amico.core.setup()

   
for iMask in range(len(allMaskNames)):
    print "Processing subject:", iMask
    dwiname = allDwiNames[iMask]
    maskname = allMaskNames[iMask]
    
    subjectFolder = os.path.dirname(dwiname)
    studyFolder = os.path.dirname(subjectFolder)
    ae = amico.Evaluation(studyFolder, subjectFolder)
    
    bvalname = bvalsnames[iMask]
    bvecname = bvecsnames[iMask]
    amico.util.fsl2scheme(bvalname,bvecname)
    
    ae.load_data(dwi_filename = os.path.basename(dwiname), 
                 scheme_filename = "NODDI_protocol.scheme", 
                 mask_filename = os.path.basename(maskname), b0_thr = 0)
    # ae.set_model("NODDI")
    ae.set_model(alg)
    ae.generate_kernels(regenerate=True)
    ae.load_kernels()
    ae.fit()
    ae.save_results()