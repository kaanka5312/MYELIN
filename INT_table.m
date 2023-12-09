% Reading each tseries file
mainFolder = '/media/kaansocat/Elements/EIB/';

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

% Each subject
for i=1:100
CIFTI = cifti_read([mainFolder,folderNames{i},'/INT/Glasser_cortical_timeseries.ptseries.nii']) ;

ACW50_arr = [] ; 
ACW0_arr = [] ;
% Calculating ACW-50 and ACW-0 for each parcellation 
for t=1:360
    [ACW0, ACW50] = ACW_kaan(CIFTI.cdata(t,:),1/0.72,100,50) ;

    ACW50_arr(t) = ACW50 ;
    ACW0_arr(t) = ACW0 ;
end
% Saving each subject as a row parcellation as columns 

save([mainFolder,folderNames{i},'/INT/INT.mat'], 'ACW0_arr' , 'ACW50_arr');
end

%%%====== AVERAGING FOR SPINNING =========

% ACW data 
ACW50_all = [] ;
ACW0_all = [] ;

for i=1:100

load([mainFolder,folderNames{i},'/INT/INT.mat']) ;
ACW50_all = vertcat(ACW50_all, ACW50_arr) ;
ACW0_all = vertcat(ACW0_all, ACW0_arr) ;

end

% Myelin data 

myelin_all = [] ;

for i=1:100

MY = load([mainFolder,folderNames{i},'/Glasser_myelin_content.txt']) ;
myelin_all = vertcat(myelin_all, MY') ;

end

save([mainFolder,'INT_all.mat'], 'ACW0_all', 'ACW50_all', 'myelin_all')

% Spin permutation testing for two cortical maps
AV_ACW50 = mean(ACW50_all,1) ;
AV_ACW0 = mean(ACW0_all,1) ;
AV_MY = mean(myelin_all,1) ;


% Git URL for the function for spin test. LATER VERIFY WITH ENIGMA TOOLBOX
% https://github.com/frantisekvasa/rotate_parcellation.git
% Coordinates for glasser parcellations (HCP)
CORD = load('/home/kaansocat/Documents/MATLAB/rotate_parcellation-master/sphere_HCP.txt') ;
perm_id = rotate_parcellation(CORD(1:180,:), CORD(181:360,:), 10.000) ;

% Generates p valie from above 
perm_sphere_p(AV_ACW0', AV_MY', perm_id, 'spearman')

% ACW0 content to myelin content
% Example data
x = AV_ACW0' ;
y = AV_MY' ;

% See the r-value
corr(x,y,'type','Spearman')

% Fit a linear regression model
lm = fitlm(x, y);

% Plot the scatter plot
scatter(x, y, 'filled', 'DisplayName', 'Data');
hold on;

% Plot the regression line
xValues = linspace(min(x), max(x), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'Regression Line');

% Plot confidence intervals
[yPred, delta] = predict(lm, xValues', 'Prediction', 'curve', 'Confidence', 'on');
plot(xValues, yPred, 'b--', 'DisplayName', '95% Confidence Interval');
plot(xValues, yPred + delta, 'b--', 'HandleVisibility', 'off');
plot(xValues, yPred - delta, 'b--', 'HandleVisibility', 'off');

% Customize the plot
xlabel('ACW-0 (seconds)','FontSize',20);
ylabel('Intracortical Myelin Content','FontSize',20);
title('Spearman Correlation (p_{spin} < 0.001, r= -0.68)', 'FontSize', 24);
legend('show');
grid on;
hold off;


% ACW50 content to myelin content
% Example data
x = AV_ACW50' ;
y = AV_MY' ;

% Fit a linear regression model
lm = fitlm(x, y);

% See the r-value
corr(x,y,'type','Spearman')

% Plot the scatter plot
scatter(x, y, 'filled', 'DisplayName', 'Data');
hold on;

% Plot the regression line
xValues = linspace(min(x), max(x), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'Regression Line');

% Plot confidence intervals
[yPred, delta] = predict(lm, xValues', 'Prediction', 'curve', 'Confidence', 'on');
plot(xValues, yPred, 'b--', 'DisplayName', '95% Confidence Interval');
plot(xValues, yPred + delta, 'b--', 'HandleVisibility', 'off');
plot(xValues, yPred - delta, 'b--', 'HandleVisibility', 'off');

% Customize the plot
xlabel('ACW-0 (seconds)','FontSize',20);
ylabel('Intracortical Myelin Content','FontSize',20);
title('Spearman Correlation (p_{spin} = 0.20, r= 0.1681)', 'FontSize', 24);
legend('show');
grid on;
hold off;
