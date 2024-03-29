# MultiSAR-Soil-Moisture-Retrieval

## Two supporting components are required 

1) The look up tables (LUT) built by the Improved IEM (Fung et al., 2002) can be found at:
https://drive.google.com/file/d/17Hnj6FOr71klS8rMjkRljb-5cVHPae5d/view?usp=sharing

2) The Genetic Algorithm TOOLBOX for inversion at https://github.com/UoS-CODeM/GA-Toolbox#genetic-algorithm-toolbox-for-matlab-v12

## Generate the synthetic data

Use syntheticCaseGenerate_v3.m to create the three synthetic data sets that used in Zhu et. al. (2020)

## Soil moisture ensemble retrieval

See Demo.m for an example

## Reference:

Zhu L., Walker J., & Shen X., (2020). Stochastic ensemble methods for multi-SAR-mission soil moisture retrieval. Remote Sensing of Environment. 251, 112099.

Zhu L., Walker J. Tsang L., Huang H., Ye N., & Rüdiger C., (2019). A multi-frequency framework for soil moisture retrieval from time series radar data. Remote Sensing of Environment , 235, 111433.

Zhu L., Walker J. Tsang L., Huang H., Ye N., & Rüdiger C., (2019). Soil moisture retrieval from time series multi angular radar data using a dry do wn constraint. Remote Sensing of Environment, 231, 111237.
