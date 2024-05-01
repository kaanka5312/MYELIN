% This script has no practical implication at the project
% Author: kaanka5312
%----%----%----%----%----%----%----%----%----%----%----%----%----%----
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
perm_sphere_p(AV_HURST, AV_GSCORR, perm_id, 'spearman') %Signif

% Hurst to Myelin
% Example data
x = AV_HURST ;
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
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'r=0.15, p_{spin} > 0.05');

% Plot confidence intervals
[yPred, delta] = predict(lm, xValues', 'Prediction', 'curve', 'Confidence', 'on');
plot(xValues, yPred, 'b--', 'DisplayName', '95% Confidence Interval');
plot(xValues, yPred + delta, 'b--', 'HandleVisibility', 'off');
plot(xValues, yPred - delta, 'b--', 'HandleVisibility', 'off');

% Customize the plot
xlabel('Hurst Exponent','FontSize',20);
ylabel('İntracranial Myelin Content','FontSize',20);
title('Global (p_{spin} > 0.05, r= 0.15)', 'FontSize', 24);
legend('show', 'FontSize', 20);
grid on;
hold off;

savefig('E:/EIB/FIGURES/MY2Hurst_Global_Hurst.fig')


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
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'r=0.729, p_{spin} < 0.05');

% Plot confidence intervals
[yPred, delta] = predict(lm, xValues', 'Prediction', 'curve', 'Confidence', 'on');
plot(xValues, yPred, 'b--', 'DisplayName', '95% Confidence Interval');
plot(xValues, yPred + delta, 'b--', 'HandleVisibility', 'off');
plot(xValues, yPred - delta, 'b--', 'HandleVisibility', 'off');

% Customize the plot
xlabel('GSCORR','FontSize',20);
ylabel('Hurst Exponent','FontSize',20);
title('Global (p_{spin} < 0.05, r= 0.729)', 'FontSize', 24);
legend('show', 'FontSize', 20);
grid on;
hold off;

savefig('E:/EIB/FIGURES/Hurst2GSCORR_Global_Hurst.fig')

%% All self regions with GSCORR 

INT = sort([108, 220, 120, 302, 286, 291, 148, 63, 258, 189]) ;
EXT = sort([138, 82, 78, 109, 318, 145, 116, 48, 292, 20, 297, 230, 57, 249]) ;
MENT = sort([241, 214, 291, 330, 248, 212, 312, 150, 107, 267, 78]) ;

%% Myelin to Hurst

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

% Spin permutation testing for two cortical maps
perm_id = rotate_parcellation(CORD(SELF(1:15),:), CORD(SELF(16:end),:), 100.000) ;
% perm_id = rotate_parcellation(CORD(INT(1:4),:), CORD(INT(5:10),:), 100.000) ;

% Generates p value from above 
perm_sphere_p(AV_HURST(SELF), AV_MY(SELF)', perm_id, 'spearman')

% GSCORR to myelin content
x = AV_HURST(SELF) ;
y = AV_MY(SELF)' ;

% See the r-value
[RHO, PVAL]=corr(x,y,'type','Spearman')

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
xlabel('GSCORR','FontSize',20);
ylabel('Intracortical Myelin Content','FontSize',20);
title('Three Layer of Self (p_{spin} > 0.05, r=0.413)', 'FontSize', 24);
legend('show');
grid on;
hold off;

savefig('E:/EIB/FIGURES/MY2Hurst_Self_Hurst.fig')

% Non Self Regions
NONSELF = 1:360 ;
NONSELF = setdiff( NONSELF, SELF) ;

perm_id = rotate_parcellation(CORD(NONSELF(1:165),:), CORD(NONSELF(166:end),:), 100.000) ;

% Generates p value from above 
perm_sphere_p(AV_HURST(NONSELF), AV_MY(NONSELF)', perm_id, 'spearman')

x = AV_HURST(NONSELF) ;
y = AV_MY(NONSELF)' ;

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
xlabel('Hurst Exponent','FontSize',20);
ylabel('Intracortical Myelin Content','FontSize',20);
title('Non-Self (p_{spin} > 0.05, r=0.13)', 'FontSize', 24);
legend('show');
grid on;
hold off;

savefig('E:/EIB/FIGURES/MY2Hurst_Nonself_Hurst.fig')

formattedTitle = sprintf('Three Layer of Self\n(p_{spin} > 0.05(\\dag), r=0.413)');

h1 = openfig('E:/EIB/FIGURES/MY2Hurst_Global_Hurst.fig','reuse'); % open figure
ax1 = gca; % get handle to axes of figure
h2 = openfig('E:/EIB/FIGURES/MY2Hurst_Self_Hurst.fig','reuse'); % open figure
ax2 = gca; % get handle to axes of figure
h3 = openfig('E:/EIB/FIGURES/MY2Hurst_Nonself_Hurst.fig','reuse');
ax3 = gca;
% test1.fig and test2.fig are the names of the figure files which you would % like to copy into multiple subplots

h4 = figure; %create new figure
s1 = subplot(1,3,1); %create and get handle to the subplot axes
title('Global (p_{spin} > 0.05, r= 0.15)', 'FontSize', 16);
ylabel('Intracortical Myelin Content','FontSize', 20)
s2 = subplot(1,3,2); %create and get handle to the subplot axes
title('Three Layer of Self (p_{spin} > 0.05(†), r=0.413)','FontSize', 14);
xlabel('Hurst Exponent', 'FontSize', 20)
s3 = subplot(1,3,3);
title('Non-Self (p_{spin} > 0.05, r=0.13)', 'FontSize', 16);
fig1 = get(ax1,'children'); %get handle to all the children in the figure
fig2 = get(ax2,'children');
fig3 = get(ax3,'children');

copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
copyobj(fig2,s2);
copyobj(fig3,s3); %copy children to new parent axes i.e. the subplot axes

%% Hurst to GSCORR

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

% Spin permutation testing for two cortical maps
perm_id = rotate_parcellation(CORD(SELF(1:15),:), CORD(SELF(16:end),:), 100.000) ;
% perm_id = rotate_parcellation(CORD(INT(1:4),:), CORD(INT(5:10),:), 100.000) ;

% Generates p value from above 
perm_sphere_p(AV_HURST(SELF), AV_GSCORR(SELF), perm_id, 'spearman')

% GSCORR to myelin content
x = AV_GSCORR(SELF) ;
y = AV_HURST(SELF) ;

% See the r-value
[RHO, PVAL]=corr(x,y,'type','Spearman')

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
xlabel('GSCORR','FontSize',20);
ylabel('Hurst Exponent','FontSize',20);
title('Three Layer of Self (p_{spin} < 0.05, r=0.77)', 'FontSize', 24);
legend('show');
grid on;
hold off;

savefig('E:/EIB/FIGURES/Hurst2GSCORR_Self_Hurst.fig')

% Non Self Regions
NONSELF = 1:360 ;
NONSELF = setdiff( NONSELF, SELF) ;

perm_id = rotate_parcellation(CORD(NONSELF(1:165),:), CORD(NONSELF(166:end),:), 100.000) ;

% Generates p value from above 
perm_sphere_p(AV_HURST(NONSELF), AV_GSCORR(NONSELF), perm_id, 'spearman')

x = AV_GSCORR(NONSELF) ;
y = AV_HURST(NONSELF) ;

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
xlabel('GSCORR','FontSize',20);
ylabel('Hurst Exponent','FontSize',20);
title('Non-Self (p_{spin} < 0.05, r=0.72)', 'FontSize', 24);
legend('show');
grid on;
hold off;

savefig('E:/EIB/FIGURES/Hurst2GSCORR_Nonself_Hurst.fig')


h1 = openfig('E:/EIB/FIGURES/Hurst2GSCORR_Global_Hurst.fig','reuse'); % open figure
ax1 = gca; % get handle to axes of figure
h2 = openfig('E:/EIB/FIGURES/Hurst2GSCORR_Self_Hurst.fig','reuse'); % open figure
ax2 = gca; % get handle to axes of figure
h3 = openfig('E:/EIB/FIGURES/Hurst2GSCORR_Nonself_Hurst.fig','reuse');
ax3 = gca;
% test1.fig and test2.fig are the names of the figure files which you would % like to copy into multiple subplots

h4 = figure; %create new figure
s1 = subplot(1,3,1); %create and get handle to the subplot axes
title('Global (p_{spin} < 0.05, r= 0.729)', 'FontSize', 16);
ylabel('Intracortical Myelin Content','FontSize', 20)
s2 = subplot(1,3,2); %create and get handle to the subplot axes
title('Three Layer of Self (p_{spin} < 0.05, r=0.77)', 'FontSize', 16);
xlabel('Hurst Exponent', 'FontSize', 20)
s3 = subplot(1,3,3);
title('Non-Self (p_{spin} < 0.05, r=0.72)', 'FontSize', 16);
fig1 = get(ax1,'children'); %get handle to all the children in the figure
fig2 = get(ax2,'children');
fig3 = get(ax3,'children');

copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
copyobj(fig2,s2);
copyobj(fig3,s3); %copy children to new parent axes i.e. the subplot axes

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
perm_sphere_p(AV_HURST(INT), AV_MY(INT)', perm_id, 'spearman')

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
scatter(AV_HURST(INT), AV_MY(INT)', 'filled','red','DisplayName','Interoceptive' )
hold on 

% Fit a linear regression model
lm = fitlm( AV_HURST(INT), AV_MY(INT)' );
% Plot the regression line
xValues = linspace(min(AV_HURST(INT)), max(AV_HURST(INT)), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'r_{INT}= 0.27, p_{spin} > 0.05');
corr(AV_HURST(INT), AV_MY(INT)','type','Spearman' );

scatter(AV_HURST(EXT), AV_MY(EXT)', 'filled','blue','DisplayName','Exteroceptive' )
hold on 

% Fit a linear regression model
lm = fitlm( AV_HURST(EXT), AV_MY(EXT)' );
% Plot the regression line
xValues = linspace(min(AV_HURST(EXT)), max(AV_HURST(EXT)), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'b-', 'LineWidth', 2, 'DisplayName', 'r_{EXT}= 0.60, p_{spin} < 0.05');
corr(AV_HURST(EXT), AV_MY(EXT)','type','Spearman' );

scatter(AV_HURST(MENT), AV_MY(MENT)', 'filled','greed', 'DisplayName','Mental' )

% Fit a linear regression model
lm = fitlm( AV_HURST(MENT), AV_MY(MENT)' );
% Plot the regression line
xValues = linspace(min(AV_HURST(MENT)), max(AV_HURST(MENT)), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'g-', 'LineWidth', 2, 'DisplayName', 'r_{MENT}= 0.67, p_{spin} > 0.05');
[RHO,PVAL]=corr(AV_HURST(MENT), AV_MY(MENT)' ,'type','Spearman' )

% Customize the plot
xlabel('Hurst Exponent','FontSize',20);
ylabel('Intracortical Myelin Content','FontSize',20);
title('Each Layer Respectively', 'FontSize', 24);
legend('show','FontSize', 20);
grid on;
hold off

