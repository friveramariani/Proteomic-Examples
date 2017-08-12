%% Example of proteomics in cancer biology with Matlab
% This report represents a represents an example Matlab proteomic data
% analysis. The data set analyze in this report can be found <https://home.ccr.cancer.gov/ncifdaproteomics/ppatterns.as here>
% which is the FDA-NCI Clinical Proteomics Program Databank. The samples
% downloaded from the FDA-NIC Proteomics Programa Databank corresponds to
% Mass-Spec of overian cancer samples: *Cancer Group* vs  *Normal Group*. 

%% loading the pre-processed dataset
load OvarianCancerQAQCdataset
whos

%% initializing variables to be used downstream in the data analysis
% workflow

N = numel(grp);                         % vector of number of samples
Cidx = strcmp('Cancer',grp);            % logical index vector for Cancer samples
Nidx = strcmp('Normal',grp);            % logical index vector for Normal samples
Cvec = find(Cidx);                      % index vector for Cancer samples
Nvec = find(Nidx);                      % index vector for Normal samples
xAxisLabel = 'Mass/Charge (M/Z)';       % x-axis label for plots
yAxisLabel = 'Ion Intensity';           % y-axis label for plots

