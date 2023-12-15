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


GSCORR_arr = zeros(360,100) ;
% Each subject
for i=1:100
CIFTI = cifti_read([mainFolder,folderNames{i},'/INT/Glasser_cortical_timeseries.ptseries.nii']) ;
GSCORR_arr(:,i)= GSCORR(CIFTI.cdata) ;

end
% Saving each subject as a row parcellation as columns 

save([mainFolder,folderNames{i},'/INT/GSCORR.mat'], 'GSCORR_arr');

%%
load('E:/EIB/INT_all.mat')

AV_MY = mean(myelin_all,1) ;
AV_GSCORR = mean(GSCORR_arr,2) ;


CORD = load('C:/Users/kaan/Documents/MATLAB/rotate_parcellation-master/rotate_parcellation-master/sphere_HCP.txt') ;
perm_id = rotate_parcellation(CORD(1:180,:), CORD(181:360,:), 10.000) ;

% Generates p valie from above 
perm_sphere_p(AV_GSCORR, AV_MY', perm_id, 'spearman')

% ACW0 content to myelin content
% Example data
x = AV_GSCORR ;
y = AV_MY' ;

% See the r-value
[RHO,PVAL]=corr(x,y,'type','Spearman')

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
%% All self regions with GSCORR 

INT = sort([108, 220, 120, 302, 286, 291, 148, 63, 258, 189]) ;
EXT = sort([138, 82, 78, 109, 318, 145, 116, 48, 292, 20, 297, 230, 57, 249]) ;
MENT = sort([241, 214, 291, 330, 248, 212, 312, 150, 107, 267, 78]) ;

%%

% Self Regions that each group is inbetween
CORD = load('C:/Users/kaan/Documents/MATLAB/rotate_parcellation-master/rotate_parcellation-master/sphere_HCP.txt') ;
SELF = sort([INT,EXT,MENT]) ;

% Droppin one of double 78 
indexToRemove = 5;

% Remove the element at the specified index
SELF = SELF([1:indexToRemove-1, indexToRemove+1:end]);

% Droppin one of double 291 
indexToRemove = 27;

% Remove the element at the specified index
SELF = SELF([1:indexToRemove-1, indexToRemove+1:end]);
% Self regions as a whole, without dividing according to layers 

perm_id = rotate_parcellation(CORD(SELF(1:15),:), CORD(SELF(16:end),:), 100.000) ;
% perm_id = rotate_parcellation(CORD(INT(1:4),:), CORD(INT(5:10),:), 100.000) ;

% Generates p value from above 
perm_sphere_p(AV_GSCORR(SELF), AV_MY(SELF)', perm_id, 'spearman')


% ACW0 content to myelin content
% Example data
% Spin permutation testing for two cortical maps

x = AV_GSCORR(SELF) ;
y = AV_MY(SELF)' ;

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
title('Three Layer of Self (p_{spin} < 0.001, r= -0.601)', 'FontSize', 24);
legend('show');
grid on;
hold off;

%% Each Layer Respectively

% Alternatively, include intersection regions to only one layer
% Including 291 to only interoceptive
INT = sort([108, 220, 120, 302, 286, 291, 148, 63, 258, 189]) ;
EXT = sort([138, 82, 78, 109, 318, 145, 116, 48, 292, 20, 297, 230, 57, 249]) ;
% Including 78 to only exteroceptive 
MENT = sort([241, 214, 330, 248, 212, 312, 150, 107, 267]) ;

% Interoceptive 
perm_id = rotate_parcellation(CORD(INT(1:4),:), CORD(INT(5:10),:), 100.000) ;
% Generates p valie from above 
perm_sphere_p(AV_GSCORR(INT), AV_MY(INT)', perm_id, 'spearman')

% Exteroceptive 
perm_id = rotate_parcellation(CORD(EXT(1:9),:), CORD(EXT(10:14),:), 100.000) ;
% Generates p valie from above 
perm_sphere_p(AV_GSCORR(EXT), AV_MY(EXT)', perm_id, 'spearman')

% Mental 
% perm_id = rotate_parcellation(CORD(MENT(1:3),:), CORD(MENT(4:11),:), 100.000) ;
% Generates p valie from above 
% perm_sphere_p(AV_GSCORR(MENT), AV_MY(MENT)', perm_id, 'spearman')

% Mental Alternative. Intersected regions included to lower levels. 
perm_id = rotate_parcellation(CORD(MENT(1:2),:), CORD(MENT(3:end),:), 100.000) ;
% Generates p valie from above 
perm_sphere_p(AV_GSCORR(MENT), AV_MY(MENT)', perm_id, 'spearman')

%%============= FIGURE OF SELF REGIONS ACW-0 to Myelin Content ============
scatter(AV_GSCORR(INT), AV_MY(INT)', 'filled','red','DisplayName','Interoceptive' )
hold on 

% Fit a linear regression model
lm = fitlm( AV_GSCORR(INT), AV_MY(INT)' );
% Plot the regression line
xValues = linspace(min(AV_GSCORR(INT)), max(AV_GSCORR(INT)), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'r_{INT}= 0.478, p_{spin} > 0.05');
corr(AV_GSCORR(INT), AV_MY(INT)','type','Spearman' );

scatter(AV_GSCORR(EXT), AV_MY(EXT)', 'filled','blue','DisplayName','Interoceptive' )
hold on 

% Fit a linear regression model
lm = fitlm( AV_GSCORR(EXT), AV_MY(EXT)' );
% Plot the regression line
xValues = linspace(min(AV_GSCORR(EXT)), max(AV_GSCORR(EXT)), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'b-', 'LineWidth', 2, 'DisplayName', 'r_{EXT}= 0.63, p_{spin} < 0.05');
corr(AV_GSCORR(EXT), AV_MY(EXT)','type','Spearman' );

scatter(AV_GSCORR(MENT), AV_MY(MENT)', 'filled','greed', 'DisplayName','Mental' )

% Fit a linear regression model
lm = fitlm( AV_GSCORR(MENT), AV_MY(MENT)' );
% Plot the regression line
xValues = linspace(min(AV_GSCORR(MENT)), max(AV_GSCORR(MENT)), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'g-', 'LineWidth', 2, 'DisplayName', 'r_{MENT}= 0.33, p_{spin} > 0.05');
[RHO,PVAL]=corr(AV_GSCORR(MENT), AV_MY(MENT)' ,'type','Spearman' )

% Customize the plot
xlabel('GSCORR','FontSize',20);
ylabel('Intracortical Myelin Content','FontSize',20);
title('Each Layer Respectively', 'FontSize', 24);
legend('show','FontSize', 20);
grid on;
hold off
