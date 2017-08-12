% preprocessing of mass-spec data
%  

repository = 'C:/Users/Felix/Dropbox/DataScience/Projects/computational biology/OvarianCD_PostQAQC/'
repositoryC = [repository 'Cancer/'];
repositoryN = [repository 'Normal/'];

filesCancer = dir([repositoryC '*.txt']);
NumberCancerDatasets = numel(filesCancer);
fprintf('Found %d Cancer mass-spectrograms.\n',NumberCancerDatasets)
filesNormal = dir([repositoryN '*.txt']);
NumberNormalDatasets = numel(filesNormal);
fprintf('Found %d Control mass-spectrograms.\n',NumberNormalDatasets)

files = [ strcat('Cancer/',{filesCancer.name}) ...
          strcat('Normal/',{filesNormal.name})];
N = numel(files);   % total number of files

fprintf('Total %d mass-spectrograms to process...\n',N)

[MZ,Y] = msbatchprocessing(repository,files);

disp('Finished; normalizing and saving to OvarianCancerQAQCdataset.mat.')
Y = msnorm(MZ,Y,'QUANTILE',0.5,'LIMITS',[3500 11000],'MAX',50);

grp = [repmat({'Cancer'},size(filesCancer));...
       repmat({'Normal'},size(filesNormal))];
   
save OvarianCancerQAQCdataset.mat Y MZ grp