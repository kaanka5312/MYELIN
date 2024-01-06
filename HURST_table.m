% Reading each tseries file
mainFolder = 'E:/EIB/100_Subj/';

% Get the directory listing
contents = dir(mainFolder);

% Initialize an empty cell array to store folder names
folderNames = {};

% Iterate through each item in the directory
for i = 1:numel(contents)
    % Check if the item is a directory (folder)
    if contents(i).isdir && ~strcmp(contents(i).name, '.') && ~strcmp(contents(i).name, '..')
        % Add the folder name to the cell array
        folderNames = [folderNames; contents(i).name];
    end
end


HURST_arr = zeros(360,100) ;

% Each subject
for i=1:100
    % Importing Hurst exponent
    HURST_arr(:,i)=csvread([mainFolder,folderNames{i},'/',num2str(folderNames{i}),'_Glasser_cortical_ei_hurst_results.csv']);
end
% Saving each subject as a row parcellation as columns 

save([mainFolder,folderNames{i},'/INT/HURST.mat'], 'HURST_arr');
