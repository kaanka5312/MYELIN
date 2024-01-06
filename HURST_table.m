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

save([mainFolder,'/HURST.mat'], 'HURST_arr');

%% Whole Brain
load('E:/EIB/100_Subj/HURST.mat')
load('E:/EIB/INT_all.mat')
load('E:/EIB/100_Subj/GSCORR.mat')

AV_MY = mean(myelin_all,1) ;
AV_HURST = mean(HURST_arr,2) ;
AV_GSCORR = mean(GSCORR_arr,2) ;

CORD = load('C:/Users/kaan/Documents/MATLAB/rotate_parcellation-master/rotate_parcellation-master/sphere_HCP.txt') ;
perm_id = rotate_parcellation(CORD(1:180,:), CORD(181:360,:), 10.000) ;

% Generates p valie from above 
perm_sphere_p(AV_HURST, AV_MY', perm_id, 'spearman') %Not Signif
perm_sphere_p(AV_GSCORR, AV_MY', perm_id, 'spearman') %Signif

subplot(1,2,1)
% Hurst to Myelin
% Example data
x = AV_GSCORR ;
y = AV_MY' ;

% See the r-value
[RHO,PVAL]=corr(x,y,'type','Spearman')

% Fit a linear regression model
lm = fitlm(x, y);

% Plot the scatter plot
scatter(x, y, 'filled', 'black', 'DisplayName', 'Data');
hold on;

% Plot the regression line
xValues = linspace(min(x), max(x), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'r=0.30, p_{spin} > 0.05');

% Plot confidence intervals
[yPred, delta] = predict(lm, xValues', 'Prediction', 'curve', 'Confidence', 'on');
plot(xValues, yPred, 'b--', 'DisplayName', '95% Confidence Interval');
plot(xValues, yPred + delta, 'b--', 'HandleVisibility', 'off');
plot(xValues, yPred - delta, 'b--', 'HandleVisibility', 'off');

% Customize the plot
xlabel('GSCORR','FontSize',20);
ylabel('İntracranial Myelin Content','FontSize',20);
title('Global (p_{spin} > 0.05, r= 0.30)', 'FontSize', 24);
legend('show', 'FontSize', 20);
grid on;
hold off;

subplot(1,2,2)
% Hurst to GSCORR
% Example data
x = AV_GSCORR ;
y = AV_HURST ;

% See the r-value
[RHO,PVAL]=corr(x,y,'type','Spearman')

% Fit a linear regression model
lm = fitlm(x, y);

% Plot the scatter plot
scatter(x, y, 'filled', 'black', 'DisplayName', 'Data');
hold on;

% Plot the regression line
xValues = linspace(min(x), max(x), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'r=0.729, p_{spin} = 0.05');

% Plot confidence intervals
[yPred, delta] = predict(lm, xValues', 'Prediction', 'curve', 'Confidence', 'on');
plot(xValues, yPred, 'b--', 'DisplayName', '95% Confidence Interval');
plot(xValues, yPred + delta, 'b--', 'HandleVisibility', 'off');
plot(xValues, yPred - delta, 'b--', 'HandleVisibility', 'off');

% Customize the plot
xlabel('GSCORR','FontSize',20);
ylabel('Hurst Exponent','FontSize',20);
title('Global (p_{spin} = 0.05, r= 0.729)', 'FontSize', 24);
legend('show', 'FontSize', 20);
grid on;
hold off;

savefig('E:/EIB/FIGURES/Global_Hurst.fig')
