% Assuming you have a 2D matrix 'fMRI_data' where rows represent grayordinates and columns represent time points
function [fisher_z_transform] = GSCORR(global_signal, fMRI_data)

% Step 1: Standardize (z-score) the fMRI signals across time points
z_scored_data = zscore(fMRI_data, 0, 2);  % Standardize along the second dimension

% Step 3: Calculate Pearson correlation (GSCORR) between the GS and each grayordinate
correlation_matrix = corr(z_scored_data', global_signal');

% Step 4: Fisher z transformation for statistical analyses
fisher_z_transform = 0.5 * log((1 + correlation_matrix) ./ (1 - correlation_matrix));

% Now 'correlation_matrix' contains Pearson correlation values, and 'fisher_z_transform'
% contains the transformed values for statistical analyses.
end
