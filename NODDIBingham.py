import sys
import os
import nibabel as nib
import numpy as np
from dipy.core.gradients import gradient_table
from dmipy.core.acquisition_scheme import gtab_dipy2mipy
from dmipy.distributions.distribute_models import SD2BinghamDistributed
from dmipy.signal_models import cylinder_models, gaussian_models
from dmipy.core.modeling_framework import MultiCompartmentModel
from dmipy.core.acquisition_scheme import acquisition_scheme_from_bvalues
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

    bingham_dispersed_bundle = SD2BinghamDistributed(models=[stick, zeppelin])

    print bingham_dispersed_bundle.parameter_names.parameter_names

    bingham_dispersed_bundle.set_tortuous_parameter('G2Zeppelin_1_lambda_perp', 'C1Stick_1_lambda_par', 'partial_volume_0')
    bingham_dispersed_bundle.set_equal_parameter('G2Zeppelin_1_lambda_par', 'C1Stick_1_lambda_par')
    bingham_dispersed_bundle.set_fixed_parameter('G2Zeppelin_1_lambda_par', 1.7e-9)

    NODDI_bingham_mod = MultiCompartmentModel(models=[ball, bingham_dispersed_bundle])

    print NODDI_bingham_mod.parameter_names

    NODDI_bingham_mod.set_fixed_parameter('G1Ball_1_lambda_iso', 3e-9)

    NODDI_fit_hcp = NODDI_bingham_mod.fit(acq_scheme, dwi, mask = mask)

    fitted_parameters = NODDI_fit_hcp.fitted_parameters

    hdr = dwi_nii.header

    # create output folder
    dir_sub = os.path.join(directory, "%03d" % iMask)
    if os.path.exists(dir_sub) == False:
        os.mkdir(dir_sub)    

    for name, values in fitted_parameters.items():
        # if values.squeeze().ndim != 2:
        #     continue
        results_nii = nib.Nifti1Image(values, dwi_nii.get_affine(), hdr)
        results_name = os.path.join(dir_sub, name + ".nii.gz")
        results_nii.to_filename(results_name)
