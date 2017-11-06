#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8

import logging

logger = logging.getLogger(__name__)

import argparse

teigen_keysdoc = {
    "note": "Note is stored to the data. It can be used for experiment identification",
    "force_rewrite": "rewrite data on the hard-drive",
    "noise_exponent": "exponent that controls the ratio of the individual components to the wavelength",
    "noise_lambda0": "Minimum and maximum noise wavelengths in millimeters",
    "noise_lambda1": "Minimum and maximum noise wavelengths in millimeters",
    "noise_std": "Standard deviation of noise intensity",
    "noise_mean": "Noise mean intensity",
    "noise_rng_seed": "noise random generator seed"
}

teigendoc = "\
Input settings<br></br>\
appearance:<br></br>\
show_surface:<br></br>\
areasampling:<br></br>\
areasize_mm:<br></br>\
areasize_px:<br></br>\
voxelsize_mm:<br></br>\
filepattern:<br></br>\
filepattern_series_number:<br></br>\
generator_id:<br></br>\
generators:<br></br>\
- - - Cylinder generator<br></br>\
      - - [element_number<br></br>\
        - [radius_distribution_uniform<br></br>\
        - [radius_distribution_normal<br></br>\
        - [radius_distribution_fixed<br></br>\
        - [radius_distribution_minimum<br></br>\
        - [radius_distribution_maximum<br></br>\
        - [radius_distribution_mean<br></br>\
        - [radius_distribution_standard_deviation<br></br>\
        - - intensity_profile_radius<br></br>\
          - &id001 [0.4, 0.7, 1.0, 1.3]<br></br>\
        - - intensity_profile_intensity<br></br>\
          - &id002 [195, 190, 200, 30]<br></br>\
        - [random_generator_seed<br></br>\
  - - Gensei generator<br></br>\
    - !!python/object/apply:collections.OrderedDict<br></br>\
      - - [n_objects<br></br>\
  - - Cylinder continues<br></br>\
    - !!python/object/apply:collections.OrderedDict<br></br>\
      - - [element_number<br></br>\
        - [radius_distribution_uniform<br></br>\
        - [radius_distribution_normal<br></br>\
        - [radius_distribution_fixed<br></br>\
        - [radius_distribution_minimum<br></br>\
        - [radius_distribution_maximum<br></br>\
        - [radius_distribution_mean<br></br>\
        - [radius_distribution_standard_deviation<br></br>\
        - - intensity_profile_radius<br></br>\
          - *id001<br></br>\
        - - intensity_profile_intensity<br></br>\
          - *id002<br></br>\
        - [random_generator_seed<br></br>\
  - - Unconnected cylinders: This generator produces objects in shape of cylinders of known length and radius. The cylinders end with hemispheres with the same radius attached to each end of the cylinders. If the cylinder length is set to zero, only the two endpoint hemispheres are generated, which results in generating a sphere. The objects do not interfere with each other.<br></br>\
element_number: number of elements at which the generator stops generating further elements, even before reaching the expected volume fraction<br></br>\
radius_distribution_uniform: all the values of radius within the given limits appear with the same probability<br></br>\
radius_distribution_normal: the values of radius will be generated using the normal (Gaussian) parametric distribution<br></br>\
radius_distribution_fixed: only the mean value of radius will be used<br></br>\
radius_distribution_minimum: the lower limit of the radius of cylinders that are to be generated<br></br>\
radius_distribution_maximum: the upper limit of the radius of cylinders that are to be generated<br></br>\
radius_distribution_mean: the mean or expectation of the normal (Gaussian) parametric distribution of the radius<br></br>\
radius_distribution_standard_deviation: the standard deviation of the normal (Gaussian) parametric distribution of the radius<br></br>\
length_distribution_mean: the mean or expectation of the normal (Gaussian) parametric distribution of the cylinder length; when length_distribution_mean and length_distribution_standard deviations are both set to 0, the generated objects become spheres<br></br>\
length_distribution_standard_deviation: the standard deviation of the normal (Gaussian) parametric distribution of the cylinder length<br></br>\
        - - intensity_profile_radius: the relative distance of pixels from the axis of each generated cylinder, where 1 denotes the radius of the cylinder, values <1 denote pixels inside the cylinder, and values >1 are pixels outside the cylinder<br></br>\
        - - intensity_profile_intensity: the intensity of pixels at the given intensity_profile_radius<br></br>\
orientation_anisotropic: if checked, cylinders will be generated not randomly and isotropically, but with preferential orientations<br></br>\
        - - orientation_main: Cartesian coordinates of the terminal point of an Euclidean vector showing the preferred orientation (active only when orientation_anisotropic is checked)<br></br>\
orientation_variance_rad: statistical variance of the angle between the orientation_main and the vectors representing the axes of generated cylinders (active only when orientation_anisotropic is checked)<br></br>\
volume_fraction: the expected volume fraction of the cylinders within the region of interest (ROI);<br></br>\
maximum_1000_iteration_number: the maximum number of thousands of iterations at which the generator stops<br></br>\
random_generator_seed: a number used to initialize the pseudorandom number generator. Can be any integer between 0 and (2^32 - 1) inclusive, an array (or other sequence) of such integers, or None (the default). If seed is None (the value is set to -1), then RandomState will try to read data from /dev/urandom (or the Windows analogue) if available or seed from the clock otherwise must be convertible to 32bit unsigned integers<br></br>\
last_element_can_be_smaller:<br></br>\
<br></br>\
Postprocessing:<br></br>\
gaussian_blur: Gaussian blur is used if this parameter is set True.<br></br>\
gaussian_filter_sigma_mm: Standard deviation for Gaussian kernel. The standard deviations of the Gaussian filter<br></br>\
add_noise<br></br>\
noise_preview<br></br>\
limit_negative_intensities<br></br>\
noise_random_generator_seed<br></br>\
exponent<br></br>\
lambda_start<br></br>\
lambda_stop<br></br>\
noise_std<br></br>\
noise_mean<br></br>\
surface_measurement<br></br>\
measurement_multiplier<br></br>\
measurement_resolution<br></br>\
output_dtype<br></br>\
negative<br></br>\
required_teigen_version: 0.2.14<br></br>\
<br></br>\
<br></br>\
<br></br>\
Outputs:<br></br>\
<br></br>\
length [mm]: analytically computed total length of actually generated cylinders (i.e., the distance between the endpoints of the cylinder axis) within the ROI<br></br>\
volume [mm^3]: analytically computed total volume of actually generated objects within the ROI; this includes both the volume of the cylinders and the volume of the hemispheres<br></br>\
surface [mm^2]: analytically computed total surface of actually generated objects with the ROI; this includes both the surface of the cylinders and the surface of the hemispheres<br></br>\
area volume [mm^3]: total volume of ROI<br></br>\
count []: total number of objects actually generated within the ROI<br></br>\
numeric volume [mm^2]: numerically computed total volume of actually generated objects within the ROI; this includes both the volume of the cylinders and the volume of the hemispheres; The calculation is based on surface triangulation of the parametrically defined objects. Cylinder and spheres primitives are connected into one object using vtkBooleanOperationPolyDataFilter (http://hdl.handle.net/10380/3262) and volume is measured using the vtkMassProperties.<br></br>\
<br></br>\
.  using the……..algorithm.<br></br>\
numeric surface [mm^2]: numerically computed total surface of actually generated objects within the ROI; this includes both the surface of the cylinders and the surface of the hemispheres; The calculation is based on surface triangulation of the parametrically defined objects using the……..algorithm.<br></br>\
length d. [mm^-2]: the length density of the actually generated cylinders i.e., the analytically computed length per volume of the ROI<br></br>\
volume d. []: the volume density (or volume fraction) of the actually generated objects (including cylinders and hemispheres) i.e., the analytically computed volume per volume of the ROI<br></br>\
surface d. [mm^-1]:  the surface density the actually generated objects (including cylinders and hemispheres) i.e., the analytically computed surface per volume of the ROI<br></br>\
point1: Cartesian coordinates of the initial point of a vector representing the axis of the generated cylinder<br></br>\
point2: Cartesian coordinates of the terminal point of a vector representing the axis of the generated cylinder<br></br>\
radius: the radius of the actually generated cylinder and its corresponding hemisphere<br></br>\
vector: triples of scalar components identifying the vector representing the axis of the generated cylinder<br></br>\
<br></br>\
<br></br>\
      "
def main():
    pass
if __name__ == "__main__":
    main()
