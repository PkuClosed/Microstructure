import sys
import os
import nibabel as nib
import numpy as np
from dipy.core.gradients import gradient_table
from dmipy.core.acquisition_scheme import gtab_dipy2mipy
from dmipy.distributions import distribute_models
from dmipy.distributions.distribute_models import SD2BinghamDistributed
from dmipy.signal_models import cylinder_models, gaussian_models
from dmipy.core.modeling_framework import MultiCompartmentModel
from dmipy.core.acquisition_scheme import acquisition_scheme_from_bvalues

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

    bvalues = np.loadtxt(allBvalsNames[iMask])  # given in s/mm^2
    bvalues_SI = bvalues * 1e6  # now given in SI units as s/m^2
    gradient_directions = np.loadtxt(allBvecsNames[iMask])  # on the unit sphere
    acq_scheme = acquisition_scheme_from_bvalues(bvalues_SI, gradient_directions)

    # gtab_dipy = gradient_table(allBvalsNames[iMask], allBvecsNames[iMask], atol=3e-2)
    # acq_scheme_mipy = gtab_dipy2mipy(gtab_dipy)

    acq_scheme.print_acquisition_info

    ball = gaussian_models.G1Ball()
    cylinder = cylinder_models.C4CylinderGaussianPhaseApproximation()
    gamma_cylinder = distribute_models.DD1GammaDistributed(models=[cylinder])

    axcaliber_gamma = MultiCompartmentModel(models=[ball, gamma_cylinder])
    print axcaliber_gamma.parameter_cardinality

    axcaliber_gamma.set_fixed_parameter('DD1GammaDistributed_1_C4CylinderGaussianPhaseApproximation_1_lambda_par', 1.7e-9)
    axcaliber_gamma.set_fixed_parameter('DD1GammaDistributed_1_C4CylinderGaussianPhaseApproximation_1_mu', [0, 0])

    axcaliber_gamma_fit = axcaliber_gamma.fit(acq_scheme, dwi, mask = mask, solver='mix', maxiter=100)

    fitted_parameters = axcaliber_gamma_fit.fitted_parameters

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
