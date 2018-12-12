# BE562 Automating Nephritis Detection in Kidney Biopsies Using Nuclei Segmentation and GMMs

Written in MATLAB R2018a

## MANA

The MANA code includes three files: MANA.m, PWM.m, and MedianArea.m

Due to lack of clarity in the MANA paper, the algorithm we reproduced can either have a slow runtime or much faster runtime. The code located here is the slower version that was run on the SCC to segment kidney tissue slices. If you would like the faster version, please contact me. 

To generate nuclei centers, run MANA.m 

MANA paper: https://biomedical-engineering-online.biomedcentral.com/articles/10.1186/s12938-018-0518-0

## GMMs

To implement the GMMs and density thresholding, run GMM_DensityThresholding.m

The GMM contains 4 files, a description of each is below. 

GMM_DensityThresholding.m:
This program takes in a set of coordinates (position of nuclei centers in the x and y coordinates) and applies GMM clustering to create contours on top of the image using the estimated parameters of each cluster. Density thresholding is then applied to create an output image with only the contours which meet the threshold density value.

findarea.m:
Helper function which finds the pixel area of a given contour.

gaussianND.m: 
Helper function which outputs the probability density function of a N dimension gaussian distribution. 

weightedAverage.m:
Helper function which calculates the weighted average of values given their weights. 

