# Algorithm
The basic concept of the algorithm includes the definition of objects to generate, the generation of the framework of the fiber structure, the surface representation, the quantitative description, the volume representation, and finally the file storage. 

The basic element used in this task is a tube. It is a cylindrical body that ends with hemispheres. The user can set parameters for object length, object radius, and parameters that affect the direction and isotropy of objects. 

## Generators

Generátor Unconnected Tubes (patrně přejmenovat na něco jako Independent Tubes).
V každé iteraci je vygenerován náhodný bod, vektor, délka a poloměr.
Náhodný bod je vždy uvnitř definované oblasti. Vektor je dán směrovými parametry. Délka a poloměr denzitní funkcí se zvolenými parametry. Dále jsou určeny kolize tohoto objektu s objekty, které byly do oblasti umístěny dříve. Pokud kolize není, dojde ke vložení objektu do oblasti a pokračuje se do další iterace. Pokud dochází ke kolizi, generátor se osmkrát pokusí zkrátit délku na 75% a pokaždé testuje kolize. V případě úspěchu dojde ke vložení objektu do oblasti, jinak se pokračuje na další iteraci. 

Postup generování je ukončen, pokud je překročena hodnota volume_fraction, nebo element_number, nebo maximální dovolený počet iterací (v tisících) přesáhne hodnotu maximum_1000_iteration_number.


# Volume and surface measurement

Pro objekty bez dotyku lze povrch a objem stanovit analyticky. 
*vzorec*
V případě objektů s dotykem je využívána numerická metoda výpočtu povrchu a objemu založená na triangulaci tvaru okraje objektu. Na základě zvoleného resolution je válcová část objektu reprezentována pravidelným mnohoúhelníkem s počtem hran rovným resolution. Mnohoúhelník je vepsán do kružnice o zvoleném poloměru. Kuloplochy na konci objektu jsou repezentovány mnohostěnem vytvořeným z trojúhelníků. Vrcholy trojúhelníků se nacházejí na průsečíku pomyslných rovnoběžek a polendíků. Počet rovnoběžek i poledníků je roven resolution.

Vzhledem k tomu, že je triangulovný objekt konvexní a je zcela vepsán do modelovaného tvaru, jeho objem a povrch je menší, než analytická hodnota. Tento nedostatek je kompenzován kompenzováním poloměru válcové a kulové plochy. Volba poloměrů je dána následujícími možnostmi.

inscribed: vepsaný mnohostěn - výchozí možnost
circumscribed: opsaný mnohostěn
average: poloměr mnohostěnu je mezi vepsanou a opsanou varianou

cylinder volume: poloměr je kompenzován na základě analytické ekvivalence objemu útvaru s podstavou pravidelného mnohoúhelníku. Stejnou hodnotou je kompenzován i poloměr kulových ploch. 
cylinder surface: poloměr je kompenzován na základě analytické ekvivalence povrchu útvaru s podstavou pravidelného mnohoúhelníku. Stejnou hodnotou je kompenzován i poloměr kulových ploch. 

cylinder volume + sphere error: poloměr kulové plochy je kompenzován na základě změřené chyby dané rozdílem objemu koule a mnohostěnu. Válcová oblast je kompenzována stejně jako u cylinder volume

cylinder surface + sphere error: poloměr kulové plochy je kompenzován na základě změřené chyby dané rozdílem povrchu koule a mnohostěnu. Válcová oblast je kompenzována stejně jako u cylinder surface

cylinder surface + sphere error: poloměr kulové plochy je kompenzován na základě změřené chyby dané rozdílem povrchu koule a mnohostěnu. Válcová oblast je kompenzována stejně jako u cylinder surface

cylinder surface + sphere error + join error: poloměr válce je kompenzován stejně jako u cylinder surface + sphere error. Stejně tak kulová plocha. Navíc je poloměr kulové plochy kompenzován naměřenou chybou spoje mezi kulovou plochou a válcovou plochou


cylinder volume + sphere error + join error: poloměr válce je kompenzován stejně jako u cylinder volume + sphere error. Stejně tak kulový objem. Navíc je poloměr kulového objemu kompenzován naměřenou chybou spoje mezi kulovým objemem a válcovým objemem




#Input settings
appearance: 
show_surface: 
areasampling:
areasize_mm: 
areasize_px: 
voxelsize_mm: 
filepattern: 
filepattern_series_number: 
generator_id: 
generators:
- - - Cylinder generator
      - - [element_number
        - [uniform_radius_distribution
        - [normal_radius_distribution
        - [fixed_radius_distribution
        - [radius_distribution_minimum
        - [radius_distribution_maximum
        - [radius_distribution_mean
        - [radius_distribution_standard_deviation
        - - intensity_profile_radius
          - &id001 [0.4, 0.7, 1.0, 1.3]
        - - intensity_profile_intensity
          - &id002 [195, 190, 200, 30]
        - [random_generator_seed
  - - Gensei generator
    - !!python/object/apply:collections.OrderedDict
      - - [n_objects
  - - Cylinder continues
    - !!python/object/apply:collections.OrderedDict
      - - [element_number
        - [uniform_radius_distribution
        - [normal_radius_distribution
        - [fixed_radius_distribution
        - [radius_distribution_minimum
        - [radius_distribution_maximum
        - [radius_distribution_mean
        - [radius_distribution_standard_deviation
        - - intensity_profile_radius
          - *id001
        - - intensity_profile_intensity
          - *id002
        - [random_generator_seed
  - - Unconnected cylinders: This generator produces objects in shape of cylinders of known length and radius. The cylinders end with hemispheres with the same radius attached to each end of the cylinders. If the cylinder length is set to zero, only the two endpoint hemispheres are generated, which results in generating a sphere. The objects do not interfere with each other. 
element_number: number of elements at which the generator stops generating further elements, even before reaching the expected volume fraction
uniform_radius_distribution: all the values of radius within the given limits appear with the same probability 
normal_radius_distribution: the values of radius will be generated using the normal (Gaussian) parametric distribution
fixed_radius_distribution: only the mean value of radius will be used
radius_distribution_minimum: the lower limit of the radius of cylinders that are to be generated
radius_distribution_maximum: the upper limit of the radius of cylinders that are to be generated
radius_distribution_mean: the mean or expectation of the normal (Gaussian) parametric distribution of the radius
radius_distribution_standard_deviation: the standard deviation of the normal (Gaussian) parametric distribution of the radius
length_distribution_mean: the mean or expectation of the normal (Gaussian) parametric distribution of the cylinder length; when length_distribution_mean and length_distribution_standard deviations are both set to 0, the generated objects become spheres
length_distribution_standard_deviation: the standard deviation of the normal (Gaussian) parametric distribution of the cylinder length
        - - intensity_profile_radius: the relative distance of pixels from the axis of each generated cylinder, where 1 denotes the radius of the cylinder, values <1 denote pixels inside the cylinder, and values >1 are pixels outside the cylinder
        - - intensity_profile_intensity: the intensity of pixels at the given intensity_profile_radius
orientation_anisotropic: if checked, cylinders will be generated not randomly and isotropically, but with preferential orientations 
        - - orientation_main: Cartesian coordinates of the terminal point of an Euclidean vector showing the preferred orientation (active only when orientation_anisotropic is checked)
orientation_variance_rad: statistical variance of the angle between the orientation_main and the vectors representing the axes of generated cylinders (active only when orientation_anisotropic is checked)
volume_fraction: the expected volume fraction of the cylinders within the region of interest (ROI); 
maximum_1000_iteration_number: the maximum number of thousands of iterations at which the generator stops
random_generator_seed: a number used to initialize the pseudorandom number generator. Can be any integer between 0 and (2^32 - 1) inclusive, an array (or other sequence) of such integers, or None (the default). If seed is None (the value is set to -1), then RandomState will try to read data from /dev/urandom (or the Windows analogue) if available or seed from the clock otherwise must be convertible to 32bit unsigned integers
last_element_can_be_smaller: 

Postprocessing:
gaussian_blur: Gaussian blur is used if this parameter is set True. 
gaussian_filter_sigma_mm: Standard deviation for Gaussian kernel. The standard deviations of the Gaussian filter 
add_noise 
noise_preview
limit_negative_intensities
noise_random_generator_seed
exponent
lambda_start
lambda_stop
noise_amplitude
noise_mean
surface_measurement
measurement_multiplier
measurement_resolution
output_dtype
negative
required_teigen_version: 0.2.14



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
