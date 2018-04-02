# -*- coding: utf-8 -*-
"""
Created on Wed May 18 16:51:56 2016

@author: cyye
"""
import amico 
import sys

#if len(sys.argv) < 7:
#    studyFolder = "/DATA/249/cyye/data/NODDI"
#    subjectFolder = "/DATA/249/cyye/data/NODDI/NODDI_example_dataset"
#    dwiName = "NODDI_DWI.img"
#    bvalName = "/DATA/249/cyye/data/NODDI/NODDI_example_dataset/NODDI_protocol.bval"
#    bvecName = "/DATA/249/cyye/data/NODDI/NODDI_example_dataset/NODDI_protocol.bvec"
#    maskName = "brain_mask.img"
if len(sys.argv) == 7:
    studyFolder = sys.argv[1]
    subjectFolder = sys.argv[2]
    dwiName = sys.argv[3]
    bvalName = sys.argv[4]
    bvecName = sys.argv[5]
    maskName = sys.argv[6]
    alg = "NODDI"
if len(sys.argv) == 8:
    studyFolder = sys.argv[1]
    subjectFolder = sys.argv[2]
    dwiName = sys.argv[3]
    bvalName = sys.argv[4]
    bvecName = sys.argv[5]
    maskName = sys.argv[6]
    sel = int(sys.argv[7])
    if sel == 0:
        alg = "NODDI"
    if sel == 1:
        alg = "CylinderZeppelinBall"
amico.core.setup()
ae = amico.Evaluation(studyFolder, subjectFolder)
amico.util.fsl2scheme(bvalName,bvecName)
ae.load_data(dwi_filename = dwiName, scheme_filename = "NODDI_protocol.scheme", mask_filename = maskName, b0_thr = 0)
# ae.set_model("NODDI")
ae.set_model(alg)
ae.generate_kernels(regenerate=True)
ae.load_kernels()
ae.fit()
ae.save_results()