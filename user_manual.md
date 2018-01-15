

    


# Algorithm
The basic concept of the algorithm includes the definition of objects to generate, the generation of the framework of the fiber structure, the surface representation, the quantitative description, the volume representation, and finally the file storage. 

The basic element used in this task is a tube. It is a cylindrical body that ends with hemispheres. The user can set parameters for object length, object radius, and parameters that affect the direction and isotropy of objects. 
Four different generators can be used.

# Graphical User Ingerface

Layout of the Graphical User Interface (GUI) is divided into three columns. 
First column is for setting all parameters that determine application output. Some of the parameters set the behavior of the user interface. 

Second column is used to trigger first and second step of generation procedure. Some of the most importatnt parameters are in this part of GUI.

Third column show the outputs. The 3D visualization, graphs and tabels are displayed here.

<img src="https://raw.githubusercontent.com/mjirik/teigen/master/graphics/paper/screenshots/teigen_screenshot1_400.png" width="800">



# Parameters and input settings

* Appearance
     * `skip_volume_generation` - can be used if slice images are not required
     * `surface_3d_preview` - enable 3D preview
     * `show_aposteriori_surface` - enable 3D aposteriori surface. This surface is based on output voxel-grid information.
       No apriori information about shape is used.
     * `noise_preview` - show 1 slice of noise
     * `force_rewrite` - rewrite output files
* Output
     * `note` - can be used for description of experiment
     * `aposteriori_measurement` - turn on aposteriori measurement. It is based on threshold segmentation
     of output data and marching cubes algorithm
     * `one_row_filenam` - CSV file where each measurement is stored as one row. It can be usefull for experiment series
* Postprocessing
     * `gaussian_blur` - add gaussian blur to output data. Gaussian blur is used if this parameter is set True. 
     * `gaussian_filter_sigma_mm` - control gaussian blur standard deviation parameter
     * `add_noise` - add noise to output data
     * `limit_negative_intensities` - output intensities under the zero value are set to zero
     * `noise_rng_seed` - random generator seed
     * `noise_exponent` - exponent that controls the ratio of the individual components to the wavelength
     * `noise_lambda0` - Minimum and maximum noise wavelengths in millimeters
     * `noise_lambda1` - Minimum and maximum noise wavelengths in millimeters
     * `noise_std` - Standard deviation of noise intensity
     * `noise_mean` - Noise mean intensity
     * `intensity_profile_radius` the relative distance of pixels from the axis of each generated cylinder, where 1 denotes the radius of the cylinder, values <1 denote pixels inside the cylinder, and values >1 are pixels outside the cylinder
     * `intensity_profile_intensity` the intensity of pixels at the given intensity_profile_radius
     * `background_intensity` - Control the background intensity of output image.
* Measurement
     * `tube_shape` - Tube-like objects are generated if set True. Otherwise cylinders are generated.
     * `polygon_radius_selection_method` - specify compensatin method
* Batch processing - allow to run multiple configuration 
* Generators
     * Unconnected Tubes - generates separated tube objects. This generator produces objects in shape of cylinders of known length and radius. The cylinders end with hemispheres with the same radius attached to each end of the cylinders. If the cylinder length is set to zero, only the two endpoint hemispheres are generated, which results in generating a sphere. The objects do not interfere with each other. 
     * Connected Tubes - generate objects with overlap
     * Unconnected Porosity - generates tube shaped cavities with no connection
     * Connected Porosity - generates tube shaped connected cavities
* Generator parameters
     * `element_number`: number of elements at which the generator stops generating further elements, even before reaching the expected volume fraction
     * `uniform_radius_distribution`: all the values of radius within the given limits appear with the same probability 
     * `normal_radius_distribution`: the values of radius will be generated using the normal (Gaussian) parametric distribution
     * `fixed_radius_distribution`: only the mean value of radius will be used
     * `radius_distribution_minimum`: the lower limit of the radius of cylinders that are to be generated
     * `radius_distribution_maximum`: the upper limit of the radius of cylinders that are to be generated
     * `radius_distribution_mean`: the mean or expectation of the normal (Gaussian) parametric distribution of the radius
     * `radius_distribution_standard_deviation`: the standard deviation of the normal (Gaussian) parametric distribution of the radius
     * `length_distribution_mean`: the mean or expectation of the normal (Gaussian) parametric distribution of the cylinder length; when length_distribution_mean and length_distribution_standard deviations are both set to 0, the generated objects become spheres
     * `length_distribution_standard_deviation`: the standard deviation of the normal (Gaussian) parametric distribution of the cylinder length
     * `orientation_anisotropic`: if checked, cylinders will be generated not randomly and isotropically, but with preferential orientations 
     * `orientation_main`: Cartesian coordinates of the terminal point of an Euclidean vector showing the preferred orientation (active only when orientation_anisotropic is checked)
     * `orientation_variance_rad`: statistical variance of the angle between the orientation_main and the vectors representing the axes of generated cylinders (active only when orientation_anisotropic is checked)
     * `volume_fraction`: the expected volume fraction of the cylinders within the region of interest (ROI); 
     * `maximum_1000_iteration_number`: the maximum number of thousands of iterations at which the generator stops
     * `random_generator_seed`: a number used to initialize the pseudorandom number generator. Can be any integer between 0 and (2^32 - 1) inclusive, an array (or other sequence) of such integers, or None (the default). If seed is None (the value is set to -1), then RandomState will try to read data from /dev/urandom (or the Windows analogue) if available or seed from the clock otherwise must be convertible to 32bit unsigned integers
     * `last_element_can_be_smaller`: 


# Outputs:

length [mm]: analytically computed total length of actually generated cylinders (i.e., the distance between the endpoints of the cylinder axis) within the ROI
volume [mm^3]: analytically computed total volume of actually generated objects within the ROI; this includes both the volume of the cylinders and the volume of the hemispheres
surface [mm^2]: analytically computed total surface of actually generated objects with the ROI; this includes both the surface of the cylinders and the surface of the hemispheres
area volume [mm^3]: total volume of ROI
count []: total number of objects actually generated within the ROI
numeric volume [mm^2]: numerically computed total volume of actually generated objects within the ROI; this includes both the volume of the cylinders and the volume of the hemispheres; The calculation is based on surface triangulation of the parametrically defined objects. Cylinder and spheres primitives are connected into one object using vtkBooleanOperationPolyDataFilter (http://hdl.handle.net/10380/3262) and volume is measured using the vtkMassProperties.

.  using the……..algorithm.   
numeric surface [mm^2]: numerically computed total surface of actually generated objects within the ROI; this includes both the surface of the cylinders and the surface of the hemispheres; The calculation is based on surface triangulation of the parametrically defined objects using the……..algorithm. 
length d. [mm^-2]: the length density of the actually generated cylinders i.e., the analytically computed length per volume of the ROI
volume d. []: the volume density (or volume fraction) of the actually generated objects (including cylinders and hemispheres) i.e., the analytically computed volume per volume of the ROI
surface d. [mm^-1]:  the surface density the actually generated objects (including cylinders and hemispheres) i.e., the analytically computed surface per volume of the ROI
point1: Cartesian coordinates of the initial point of a vector representing the axis of the generated cylinder
point2: Cartesian coordinates of the terminal point of a vector representing the axis of the generated cylinder
radius: the radius of the actually generated cylinder and its corresponding hemisphere
vector: triples of scalar components identifying the vector representing the axis of the generated cylinder




# Command line

    python teigen.__main__ -p parameters.yaml
