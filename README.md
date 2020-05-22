# MISTy pipelines
MISTy pipelines used to generate results for the paper


## Example pipeline

### Running MISTy
To run MISTy on the provided in silico data, install the MISTy package and execute the script synthetic_pipeline.R

### Processing time
The processing of all 10 in silico slidies takes approximately 3 minutes on a standard configuration using 4 threads in parralel.

### Output
The results for each sample are stored in separate folders in the folder results/synthetic. In each folder the file coefficients.txt stores the values of the fusion parameters from the meta-modelfor each view and the corresponding p-value for each parameter. The file performance.txt stores the values of the performance measured with root mean squared error (RMSE) and the variance explained (R2) for a model containing only the basic, intrinsic view and for a model considering all defined views. The p values reported in this file are obtained by a one sided t-test comparing the performance per target using the values for all folds. The raw values of the importances of each marker as predictors are stored in separate files, one for each target marker.