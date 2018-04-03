import sys
import os
import nibabel as nib
import numpy as np
from dipy.core.gradients import gradient_table
from dmipy.core.acquisition_scheme import gtab_dipy2mipy
from dmipy.distributions.distribute_models import SD1WatsonDistributed
from dmipy.signal_models import cylinder_models, gaussian_models
from dmipy.core.modeling_framework import MultiCompartmentModel
from dmipy.core.acquisition_scheme import acquisition_scheme_from_bvalues
from dipy.data import get_sphere
from dmipy.core.acquisition_scheme import gtab_dipy2mipy

# Load data
if len(sys.argv) > 1:
    dwinames = sys.argv[1]
    masknames = sys.argv[2]
    bvecsnames = sys.argv[3]
    bvalsnames = sys.argv[4]
    directory = sys.argv[5]

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

if os.path.exists(directory) == False:
    os.mkdir(directory)
# Processing

for iMask in range(len(allMaskNames)):
    print "Processing subject", iMask

    print "Loading"
    dwi_nii = nib.load(allDwiNames[iMask])
    dwi = dwi_nii.get_data()
    mask = nib.load(allMaskNames[iMask]).get_data()

    gradient_directions = np.loadtxt(allBvecsNames[iMask])  # on the unit sphere
    bvalues = np.loadtxt(allBvalsNames[iMask])  # given in s/mm^2
    gtab_dipy = gradient_table(bvalues, gradient_directions, atol=3e-2)
    acq_scheme = gtab_dipy2mipy(gtab_dipy)

    acq_scheme.print_acquisition_info

    ball = gaussian_models.G1Ball()
    stick = cylinder_models.C1Stick()
    zeppelin = gaussian_models.G2Zeppelin()

    watson_dispersed_bundle1 = SD1WatsonDistributed(models=[stick, zeppelin])

    print watson_dispersed_bundle1.parameter_names

    watson_dispersed_bundle1.set_tortuous_parameter('G2Zeppelin_1_lambda_perp','C1Stick_1_lambda_par','partial_volume_0')
    watson_dispersed_bundle1.set_equal_parameter('G2Zeppelin_1_lambda_par', 'C1Stick_1_lambda_par')
    watson_dispersed_bundle1.set_fixed_parameter('G2Zeppelin_1_lambda_par', 1.7e-9)

    watson_dispersed_bundle2 = watson_dispersed_bundle1.copy()
    NODDIx_mod = MultiCompartmentModel(models=[ball, watson_dispersed_bundle1, watson_dispersed_bundle2])

    print NODDIx_mod.parameter_names

    NODDIx_mod.set_fixed_parameter('G1Ball_1_lambda_iso', 3e-9)

    NODDIx_fit = NODDIx_mod.fit(acq_scheme, dwi, mask = mask, solver='mix', maxiter=300)

    hdr = dwi_nii.header

    # create output folder
    dir_sub = os.path.join(directory, "%03d" % iMask)
    if os.path.exists(dir_sub) == False:
        os.mkdir(dir_sub)
        
    sphere = get_sphere(name='repulsion200').subdivide()
    fods = NODDIx_fit.fod(sphere.vertices, visual_odi_lower_bound=0.05)
    results_nii = nib.Nifti1Image(fods, dwi_nii.get_affine(), hdr)
    results_name = os.path.join(dir_sub, "FODs.nii.gz")
    results_nii.to_filename(results_name)

    fitted_parameters = NODDIx_fit.fitted_parameters
    for name, values in fitted_parameters.items():
        # if values.squeeze().ndim != 2:
        #     continue

        results_nii = nib.Nifti1Image(values, dwi_nii.get_affine(), hdr)
        results_name = os.path.join(dir_sub, name + ".nii.gz")
        results_nii.to_filename(results_name)
