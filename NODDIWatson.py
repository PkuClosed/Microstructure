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

# Load data
if len(sys.argv) >= 10:
    dwinames = sys.argv[1]
    masknames = sys.argv[2]
    bvecsnames = sys.argv[3]
    bvalsnames = sys.argv[4]
    odnames = sys.argv[5]
    directory = sys.argv[9]

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

# Processing

for iMask in range(len(allMaskNames)):
    print "Processing subject", iMask

    print "Loading"
    dwi_nii = nib.load(allDwiNames[iMask])
    dwi = dwi_nii.get_data()
    mask = nib.load(allMaskNames[iMask]).get_data()

    bvalues = np.loadtxt(allBvalsNames[iMask])  # given in s/mm^2
    bvalues_SI = bvalues * 1e6  # now given in SI units as s/m^2
    gradient_directions = np.loadtxt(allBvecsNames[iMask])  # on the unit sphere
    acq_scheme = acquisition_scheme_from_bvalues(bvalues_SI, gradient_directions)

    # gtab_dipy = gradient_table(allBvalsNames[iMask], allBvecsNames[iMask], atol=3e-2)
    # acq_scheme_mipy = gtab_dipy2mipy(gtab_dipy)

    acq_scheme.print_acquisition_info

    ball = gaussian_models.G1Ball()
    stick = cylinder_models.C1Stick()
    zeppelin = gaussian_models.G2Zeppelin()

    watson_dispersed_bundle = SD1WatsonDistributed(models=[stick, zeppelin])

    print watson_dispersed_bundle.parameter_names

    watson_dispersed_bundle.set_tortuous_parameter('G2Zeppelin_1_lambda_perp','C1Stick_1_lambda_par','partial_volume_0')
    watson_dispersed_bundle.set_equal_parameter('G2Zeppelin_1_lambda_par', 'C1Stick_1_lambda_par')
    watson_dispersed_bundle.set_fixed_parameter('G2Zeppelin_1_lambda_par', 1.7e-9)

    NODDI_mod = MultiCompartmentModel(models=[ball, watson_dispersed_bundle])

    print NODDI_mod.parameter_names

    NODDI_mod.set_fixed_parameter('G1Ball_1_lambda_iso', 3e-9)

    NODDI_fit_hcp = NODDI_mod.fit(acq_scheme, dwi, mask = mask)

    fitted_parameters = NODDI_fit_hcp.fitted_parameters

    hdr = dwi_nii.header


    for name, values in fitted_parameters.items():
        # if values.squeeze().ndim != 2:
        #     continue
        results_nii = nib.Nifti1Image(values, dwi_nii.get_affine(), hdr)
        results_name = os.path.join(directory, name + "%02d" % iMask + ".nii.gz")
        results_nii.to_filename(results_name)
