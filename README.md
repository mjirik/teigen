
[![Build Status](https://travis-ci.org/mjirik/teigen.svg?branch=master)](https://travis-ci.org/mjirik/teigen)
[![Coverage Status](https://coveralls.io/repos/github/mjirik/teigen/badge.svg?branch=master)](https://coveralls.io/github/mjirik/teigen?branch=master)

# teigen
Test Image Generator. The basic concept of the algorithm includes the definition of objects to generate, the generation of the framework of the fiber structure, the surface representation, the quantitative description, the volume representation, and finally the file storage.

[user manual](https://github.com/mjirik/teigen/blob/master/user_manual.md)

# Installation

On Linux, OS X or Windows the conda system can be used:

    conda install -c mjirik -c simpleitk -c conda-forge teigen

Or you can use [Windows installer](http://home.zcu.cz/~mjirik/lisa/install/setup_teigen.exe)

# Parameters


python teigen.__main__ -p parameters.yaml

Detailed description of paremeters can be found in 
[user manual](https://github.com/mjirik/teigen/blob/master/user_manual.md)

# Paper datasets and figures

[All experiment configuration](https://github.com/mjirik/teigen/blob/master/examples/paper_experiments_params.ipynb)

[Run all experimerimens and generate all data](https://github.com/mjirik/teigen/blob/master/examples/paper_run_experiments.ipynb)

[Figures generator script](https://github.com/mjirik/teigen/blob/master/examples/paper_figures.ipynb)


* [Dataset1](https://raw.githubusercontent.com/mjirik/teigen/master/data/Dataset1.csv) table was used for measure on the micro-CT console. 
* [Dataset2](https://raw.githubusercontent.com/mjirik/teigen/master/data/Dataset2.csv) was used for comparing known surfaces and volumes with numerically estimated surfaces and volumes. 
* [Dataset3](https://raw.githubusercontent.com/mjirik/teigen/master/data/Dataset3.csv) was used for estimating the sensitivity of the estimates upon various numbers of objects and resolution of the measurement.

# Sample output

<img src="https://raw.githubusercontent.com/mjirik/teigen/master/graphics/teigen_volume_fraction_22_unconnected_n335_paraview.png" width="200">
<img src="https://raw.githubusercontent.com/mjirik/teigen/master/graphics/teigen_volume_fraction_48_overlap4_n400_paraview.png" width="200">
<img src="https://raw.githubusercontent.com/mjirik/teigen/master/graphics/teigen_volume_fraction_89_porosities_n400_paraview.png" width="200">
On the first image is shown Unconnected tubes generator with volume fraction 22 %. Number of elements is 335.
Middle image contain tubes with overlap. Volume fraction 48 %. Number of elements is 400.
Bottom image visualize porosities with no overlap. Volume fraction 89 %, number of elements is 400.


# Screenshots

<img src="https://raw.githubusercontent.com/mjirik/teigen/master/graphics/paper/screenshots/teigen_screenshot1_400.png" width="800">
<img src="https://raw.githubusercontent.com/mjirik/teigen/master/graphics/paper/screenshots/teigen_screenshot2.png" width="800">
<img src="https://raw.githubusercontent.com/mjirik/teigen/master/graphics/paper/screenshots/teigen_screenshot3.png" width="800">
<img src="https://raw.githubusercontent.com/mjirik/teigen/master/graphics/paper/screenshots/teigen_screenshot4.png" width="800">
