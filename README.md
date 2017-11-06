
[![Build Status](https://travis-ci.org/mjirik/teigen.svg?branch=master)](https://travis-ci.org/mjirik/teigen)
[![Coverage Status](https://coveralls.io/repos/github/mjirik/teigen/badge.svg?branch=master)](https://coveralls.io/github/mjirik/teigen?branch=master)

# teigen
Test Image Generator. The basic concept of the algorithm includes the definition of objects to generate, the generation of the framework of the fiber structure, the surface representation, the quantitative description, the volume representation, and finally the file storage.

[user manual](https://github.com/mjirik/teigen/blob/master/user_manual.md)

# Installation

On Linux, OS X or Windows the conda system can be used:

    conda install -c mjirik -c simpleitk teigen

Or you can use [Windows installer](http://147.228.240.61/queetech/install/setup_teigen.exe)

# Parameters


python teigen.__main__ -p parameters.yaml

Detailed description of paremeters can be found in 
[user manual](https://github.com/mjirik/teigen/blob/master/user_manual.md)

# Paper datasets and figures

[All experiments setup](https://github.com/mjirik/teigen/blob/master/examples/paper_experiments_params.ipynb)

[Run all experimerimens and generate all data](https://github.com/mjirik/teigen/blob/master/examples/paper_run_experiments.ipynb)

[Figures generator script](https://github.com/mjirik/teigen/blob/master/examples/paper_figures.ipynb)


* [Dataset1](https://raw.githubusercontent.com/mjirik/teigen/master/data/Dataset1.csv) was used for generating image stacks that were measured on the micro-CT console. 
* [Dataset2](https://raw.githubusercontent.com/mjirik/teigen/master/data/Dataset2.csv) was used for generating image stacks for comparing known surfaces and volumes with numerically estimated surfaces and volumes. 
* [Dataset3](https://raw.githubusercontent.com/mjirik/teigen/master/data/Dataset3.csv) was used for generating image stacks for estimating the sensitivity of the estimates upon various numbers of objects and resolution of the measurement.

# Sample output

Uncoonected tubes with volume fraction 22 %. Number of elements is 335.
![unconnected](https://raw.githubusercontent.com/mjirik/teigen/master/graphics/teigen_volume_fraction_22_unconnected_n335_paraview.png)

Tubes with overlap. Volume fraction 48 %. Number of elements is 400.
![connected](https://raw.githubusercontent.com/mjirik/teigen/master/graphics/teigen_volume_fraction_48_overlap4_n400_paraview.png)
