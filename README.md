
[![Build Status](https://travis-ci.org/mjirik/teigen.svg?branch=master)](https://travis-ci.org/mjirik/teigen)
[![Coverage Status](https://coveralls.io/repos/github/mjirik/teigen/badge.svg?branch=master)](https://coveralls.io/github/mjirik/teigen?branch=master)

# teigen
Test Image Generator

# Installation

    conda install -c mjirik -c simpleitk teigen

## Windows

[Windows installer](http://147.228.240.61/queetech/install/setup_teigen.exe)

# Parameters

python teigen.__main__ -p parameters.yaml

# Sample output

![unconnected](https://raw.githubusercontent.com/mjirik/teigen/master/graphics/teigen_volume_fraction_22_unconnected_n400_paraview.png)


# Paper datasets and figures

[All experiments setup](https://github.com/mjirik/teigen/blob/master/examples/paper_experiments_params.ipynb)

[Run all experimerimens and generate all data](https://github.com/mjirik/teigen/blob/master/examples/paper_run_experiments.ipynb)

[Figures generator script](https://github.com/mjirik/teigen/blob/master/examples/paper_figures.ipynb)


* [Dataset1](https://raw.githubusercontent.com/mjirik/teigen/master/data/Dataset1.csv) was used for generating image stacks that were measured on the micro-CT console. 
* [Dataset2](https://raw.githubusercontent.com/mjirik/teigen/master/data/Dataset2.csv) was used for generating image stacks for comparing known surfaces and volumes with numerically estimated surfaces and volumes. 
* [Dataset3](https://raw.githubusercontent.com/mjirik/teigen/master/data/Dataset3.csv) was used for generating image stacks for estimating the sensitivity of the estimates upon various numbers of objects and resolution of the measurement.