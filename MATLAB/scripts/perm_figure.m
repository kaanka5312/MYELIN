% This script makes spin permutation test, and creates respective figures for paper
% Author: kaanka5312

% Adding toolboxes and data paths
addpath("G:/Drive'ım/Research/Myelin_Self/")
addpath("C:/Users/kaan/Documents/MATLAB/rotate_parcellation-master/rotate_parcellation-master/Matlab/")
% loading the data 
load( "E:/EIB/MATLAB_ClassificationApp/MED_TABLE.mat" )

dataTable = movevars(dataTable, 'GSCORR','After','ACW') ;

dataTable.Class1(dataTable.Class1==1)=0;
dataTable.Class1(dataTable.Class1==2)=1;

CORD = load('C:/Users/kaan/Documents/MATLAB/rotate_parcellation-master/rotate_parcellation-master/sphere_HCP.txt') ;
%% GLOBAL
perm_id_g = rotate_parcellation(CORD(1:180,:), CORD(181:360,:), 10000 ) ;

% Generates p value from above 
%[GS_MY_p, GS_MY_d] = perm_sphere_p(MED(:,1), MED(:,3), perm_id_g, 'spearman') ; %GSCORR and MY
%[ACW_MY_p, ACW_MY_d] = perm_sphere_p(MED(:,2), MED(:,3), perm_id_g, 'spearman') ; %GSCORR and ACW

[GS_MY_p, GS_MY_d] = perm_sphere_p(dataTable.GSCORR, dataTable.Myelin, perm_id_g, 'spearman') ; %GSCORR and MY
[ACW_MY_p, ACW_MY_d] = perm_sphere_p(dataTable.ACW, dataTable.Myelin, perm_id_g, 'spearman') ; %ACW and MY
[ACW_GS_p, ACW_GS_d] = perm_sphere_p(dataTable.ACW, dataTable.GSCORR, perm_id_g, 'spearman') ; %ACW and MY

% R values 
[rvals.GSCORR_MY ,~] = corr(dataTable.GSCORR,dataTable.Myelin,'type','Spearman') ;
[rvals.ACW_MY ,~] = corr(dataTable.ACW,dataTable.Myelin,'type','Spearman') ;
[rvals.ACW_GSCORR ,~] = corr(dataTable.ACW,dataTable.GSCORR,'type','Spearman') ;


p_and_d =  cell2struct({[ACW_MY_p; ACW_MY_d], [GS_MY_p; GS_MY_d], [ACW_GS_p; ACW_GS_d]}, ...
                       { 'ACW_MY', 'GSCORR_MY', 'ACW_GSCORR' ...
                        }, 2);


f = figure,
    set(gcf,'color','w');
    set(gcf,'units','normalized','position',[0 0 0.3 0.3])
    fns = fieldnames(p_and_d);

    for k = 1:(numel(fieldnames(rvals))-1)
        % Define plot colors
        if k <= 3; col = [0.66 0.13 0.11]; else; col = [0.2 0.33 0.49]; end

        % Plot null distributions
        t = sgtitle('Global');
        t.FontWeight = 'bold'   ;


        axs = subplot(2, 2, k); hold on
        h = histogram(p_and_d.(fns{k})(2:end), 50, 'Normalization', 'pdf', 'edgecolor', 'w', ...
                      'facecolor', col, 'facealpha', 1, 'linewidth', 0.5);
        l = line([rvals.(fns{k}) rvals.(fns{k})], get(gca, 'ylim'), 'linestyle', '--', ...
                 'color', 'k', 'linewidth', 1.5);
        xlabel(['Null correlations' newline '(' strrep(fns{k}, '_', ' ') ')'])
        ylabel('Density')
        legend(l,['{\it r}=' num2str(round(rvals.(fns{k}), 2)) newline ...
                  '{\it p}=' num2str(round(p_and_d.(fns{k})(1), 3))])
        legend boxoff

    end

            subplot(2,2,3)

        % Fit a linear regression model
lm = fitlm(dataTable.Myelin, dataTable.ACW);

% Plot the scatter plot
scatter(dataTable.Myelin, dataTable.ACW, 'filled', 'black', 'DisplayName', 'Data');
hold on;

% Plot the regression line
xValues = linspace(min(dataTable.Myelin), max(dataTable.Myelin), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'r= -0.68, p_{spin} < 0.05');

% Plot confidence intervals
[yPred, delta] = predict(lm, xValues', 'Prediction', 'curve', 'Confidence', 'on');
plot(xValues, yPred, 'b--', 'DisplayName', '95% Confidence Interval');
plot(xValues, yPred + delta, 'b--', 'HandleVisibility', 'off');
plot(xValues, yPred - delta, 'b--', 'HandleVisibility', 'off');

% Customize the plot
xlabel('ACW-0 (seconds)','FontSize',10);
ylabel('Intracortical Myelin Content','FontSize',10);
title('Spearman Correlation (p_{spin} < 0.05, r= -0.68)', 'FontSize', 14);
legend('show', 'FontSize', 10);
grid on;
hold off;

        
    subplot(2,2,4)

        % Fit a linear regression model
lm = fitlm(dataTable.Myelin, dataTable.GSCORR);

% Plot the scatter plot
scatter(dataTable.Myelin, dataTable.GSCORR, 'filled', 'black', 'DisplayName', 'Data');
hold on;

% Plot the regression line
xValues = linspace(min(dataTable.Myelin), max(dataTable.Myelin), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'r=0.309, p_{spin} = 0.059');

% Plot confidence intervals
[yPred, delta] = predict(lm, xValues', 'Prediction', 'curve', 'Confidence', 'on');
plot(xValues, yPred, 'b--', 'DisplayName', '95% Confidence Interval');
plot(xValues, yPred + delta, 'b--', 'HandleVisibility', 'off');
plot(xValues, yPred - delta, 'b--', 'HandleVisibility', 'off');

% Customize the plot
xlabel('GSCORR','FontSize',10);
ylabel('Intracortical Myelin Content','FontSize',10);
title('Spearman Correlation (p_{spin} = 0.059, r= 0.309)', 'FontSize', 14);
legend('show', 'FontSize', 10);
grid on;
hold off;

%{
        % Fit a linear regression model
lm = fitlm(dataTable.ACW, dataTable.GSCORR);

% Plot the scatter plot
scatter(dataTable.ACW, dataTable.GSCORR, 'filled', 'black', 'DisplayName', 'Data');
hold on;

% Plot the regression line
xValues = linspace(min(dataTable.Myelin), max(dataTable.Myelin), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'r=-0.3431, p_{spin} < 0.05');

% Plot confidence intervals
[yPred, delta] = predict(lm, xValues', 'Prediction', 'curve', 'Confidence', 'on');
plot(xValues, yPred, 'b--', 'DisplayName', '95% Confidence Interval');
plot(xValues, yPred + delta, 'b--', 'HandleVisibility', 'off');
plot(xValues, yPred - delta, 'b--', 'HandleVisibility', 'off');

% Customize the plot
xlabel('GSCORR','FontSize',10);
ylabel('Intracortical Myelin Content','FontSize',10);
title('Spearman Correlation (p_{spin} < 0.05, r= -0.3431)', 'FontSize', 14);
legend('show', 'FontSize', 10);
grid on;
hold off;

%}

print('E:/EIB/FIGURES/Global_Dist','-dpng','-r300');

%% Self and Nonself
INT = sort([108, 220, 120, 302, 286, 291, 148, 63, 258, 189]) ;
EXT = sort([138, 82, 78, 109, 318, 145, 116, 48, 292, 20, 297, 230, 57, 249]) ;
MENT = sort([241, 214, 291, 330, 248, 212, 312, 150, 107, 267, 78]) ;

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

%%
% Spin permutation testing for two cortical maps
perm_id_self = rotate_parcellation(CORD(SELF(1:15),:), CORD(SELF(16:end),:), 10000) ;

NONSELF = 1:360 ;
NONSELF = setdiff( NONSELF, SELF) ;
perm_id_nonself = rotate_parcellation(CORD(NONSELF(1:165),:), CORD(NONSELF(166:end),:), 10000) ;
 
% = + = + = + = + = + = + S E L F + = + = + = + = + = + = +
[GS_MY_s_p, GS_MY_s_d] = perm_sphere_p(dataTable.GSCORR(logical(dataTable.Class1),:),...
    dataTable.Myelin(logical(dataTable.Class1),:), perm_id_self, 'spearman') ; %GSCORR and MY

[ACW_MY_s_p, ACW_MY_s_d] = perm_sphere_p(dataTable.ACW(logical(dataTable.Class1),:),...
    dataTable.Myelin(logical(dataTable.Class1),:), perm_id_self, 'spearman') ; %ACW and MY

%[ACW_GSCORR_s_p, ACW_GSCORR_s_d] = perm_sphere_p(dataTable.ACW(logical(dataTable.Class1),:),...
%    dataTable.GSCORR(logical(dataTable.Class1),:), perm_id_self, 'spearman') ; %ACW and GSCORR

% R values 
[rvals.SELF_GSCORR_MY ,~] = corr(dataTable.GSCORR(logical(dataTable.Class1),:), ...
    dataTable.Myelin(logical(dataTable.Class1),:),'type','Spearman') ;

[rvals.SELF_ACW_MY ,~] = corr(dataTable.ACW(logical(dataTable.Class1),:), ...
    dataTable.Myelin(logical(dataTable.Class1),:),'type','Spearman') ;

%[rvals.ACW_GSCORR_s ,~] = corr(dataTable.ACW(logical(dataTable.Class1),:), ...
%    dataTable.GSCORR(logical(dataTable.Class1),:),'type','Spearman') ;

% = + = + = + = + = + = + N O N - S E L F + = + = + = + = + = + = +

[GS_MY_ns_p, GS_MY_ns_d] = perm_sphere_p(dataTable.GSCORR(~logical(dataTable.Class1),:),...
    dataTable.Myelin(~logical(dataTable.Class1),:), perm_id_nonself, 'spearman') ; %GSCORR and MY

[ACW_MY_ns_p, ACW_MY_ns_d] = perm_sphere_p(dataTable.ACW(~logical(dataTable.Class1),:),...
    dataTable.Myelin(~logical(dataTable.Class1),:), perm_id_nonself, 'spearman') ; %ACW and MY

%[ACW_GSCORR_ns_p, ACW_GSCORR_ns_d] = perm_sphere_p(dataTable.ACW(~logical(dataTable.Class1),:),...
%    dataTable.GSCORR(~logical(dataTable.Class1),:), perm_id_nonself, 'spearman') ; %ACW and GSCORR

% R values 
[rvals.NONSELF_GSCORR_MY ,~] = corr(dataTable.GSCORR(~logical(dataTable.Class1),:), ...
    dataTable.Myelin(~logical(dataTable.Class1),:),'type','Spearman') ;

[rvals.NONSELF_ACW_MY ,~] = corr(dataTable.ACW(~logical(dataTable.Class1),:), ...
    dataTable.Myelin(~logical(dataTable.Class1),:),'type','Spearman') ;

%[rvals.ACW_GSCORR_ns ,~] = corr(dataTable.ACW(~logical(dataTable.Class1),:), ...
%    dataTable.GSCORR(~logical(dataTable.Class1),:),'type','Spearman') ;


p_and_d =  cell2struct({[ACW_MY_s_p; ACW_MY_s_d], [GS_MY_s_p; GS_MY_s_d], ...
    [ACW_MY_ns_p; ACW_MY_ns_d], [GS_MY_ns_p; GS_MY_ns_d] }, ...
                       { 'SELF_ACW_MY', 'SELF_GSCORR_MY', ...
                       'NONSELF_ACW_MY', 'NONSELF_GSCORR_MY'
                        }, 2);


f = figure,
    set(gcf,'color','w');
    set(gcf,'units','normalized','position',[0 0 0.3 0.3])
    fns = fieldnames(p_and_d);

    for k = 1:(numel(fieldnames(rvals)))
        % Define plot colors
        if k <= 2; col = [0.66 0.13 0.11]; else; col = [0.2 0.33 0.49]; end

        % Plot null distributions
        t = sgtitle('Global');
        t.FontWeight = 'bold'   ;


        axs = subplot(2, 4, k); hold on
        h = histogram(p_and_d.(fns{k})(2:end), 50, 'Normalization', 'pdf', 'edgecolor', 'w', ...
                      'facecolor', col, 'facealpha', 1, 'linewidth', 0.5);
        l = line([rvals.(fns{k}) rvals.(fns{k})], get(gca, 'ylim'), 'linestyle', '--', ...
                 'color', 'k', 'linewidth', 1.5);
        xlabel(['Null correlations' newline '(' strrep(fns{k}, '_', ' ') ')'])
        ylabel('Density')
        legend(l,['{\it r}=' num2str(round(rvals.(fns{k}), 2)) newline ...
                  '{\it p}=' num2str(round(p_and_d.(fns{k})(1), 3))])
        legend boxoff

    end


xTitle = {'ACW', 'GSCORR'} ;
    for k = 1:(numel(fieldnames(rvals)))
        % Define plot colors
        if k <= 2; col = [0.66 0.13 0.11]; else; col = [0.2 0.33 0.49]; end

        % Plot null distributions
        t = sgtitle('Self and Non-Self');
        t.FontWeight = 'bold'   ;

        axs = subplot(2, 4, k+4); hold on
        r = rvals.(fns{k}) ;
        p = p_and_d.(fns{k})(1) ;
        index = mod(k-1, numel(xTitle)) + 1;
        if k == 1 || k == 2
        f = fig_subplot(dataTable.Myelin(logical(dataTable.Class1)), ...
            table2array(dataTable(logical(dataTable.Class1),index)), col, r, p, xTitle{index}) ;
        else
        f = fig_subplot(dataTable.Myelin(~logical(dataTable.Class1)), ...
            table2array(dataTable(~logical(dataTable.Class1),index)), col, r, p, xTitle{index}) ;
        end
    end

print('E:/EIB/FIGURES/Self_NonSelf','-dpng','-r300');

