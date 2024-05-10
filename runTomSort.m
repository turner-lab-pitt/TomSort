function runTomSort(datafolder,datalist,savefolder,workfolder)

% Input Argument
%   datafolder: folder where the extracellular data are stored
%     datalist: list of extracellular data (cell format)
%   savefolder: folder to store the results from TomSort
%   workfolder: temporal working folder

if ~exist('workfolder','var')
    workfolder = uigetdir('','Select a work folder');
    if workfolder == 0
        mkdir('tempwork');
        workfolder = [pwd,filesep,tempwork,filesep];
    else
        workfolder = [workfolder,filesep];
    end
end
if ~exist('savefolder','var')
    savefolder = uigetdir('','Select a folder to store the sort results');
    savefolder = [savefolder,filesep];
end
if ~exist('datalist','var')
    [datalist,datafolder] = uigetfile('*.mat','Select one or more mat files to sort','MultiSelect','on');
end
if ~iscell(datalist) 
    datalist = cellstr(datalist); 
end
current_folder = pwd;
cd(workfolder);

ndata = length(datalist);
for nd = 1:ndata
    clear SorterOutput numUnits output
    output = sort_channel([datafolder,datalist{nd}]);
    SorterOutput = output;
    numUnits = output.sort_results.num_units;
    disp(['Writing ',savefolder,datalist{nd}])
    save([savefolder,datalist{nd}],'SorterOutput','numUnits','-v7.3');
end

cd(current_folder);


