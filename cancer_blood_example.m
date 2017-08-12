% This code is an example of proteomics in cancer biology

% loading the pre-processed dataset
load OvarianCancerQAQCdataset
whos

% initializing variables to be used downstream in the data analysis
% workflow

N = numel(grp);                         % vector of number of samples
Cidx = strcmp('Cancer',grp);            % logical index vector for Cancer samples
Nidx = strcmp('Normal',grp);            % logical index vector for Normal samples
Cvec = find(Cidx);                      % index vector for Cancer samples
Nvec = find(Nidx);                      % index vector for Normal samples
xAxisLabel = 'Mass/Charge (M/Z)';       % x-axis label for plots
yAxisLabel = 'Ion Intensity';           % y-axis label for plots

