function [final_clusters, final_epsilons, final_fracs] = clustering(jump_stats, jumps_measured, desired_fraction, eps_range, select_moments)
% Specify the desired fraction of data set as 
% desired_fraction = [.6 .3];
% Specify range of epsilon values to sweep over
% eps_range = linspace(0.000,1,1000);

% Color used for various plots
purple = [61 38 168]/255;   % purple
orange = [217 83 25]/255;   % orange
green = [34 177 76]/255;   % green
blue = [39 150 235]/255;  % blue
yellow = [255 201 14]/255; % yellow
black = [0 0 0]/255;   % black
cyan = [77 190 238]/255;  % cyan

cmap = [blue;
  green; ...
  orange; ...
  purple];
set(groot,'defaultAxesColorOrder',cmap)

% X axis = Feature 1, Y Axis = Feature 2
Feature1 = jump_stats(:,2); 
Feature2 = jump_stats(:,3);
Feature3 = jump_stats(:,4);
Feature4 = jump_stats(:,5);
Feature5 = jump_stats(:,8);

% NormFeature1 =  (Feature1 - min(Feature1)) / ( max(Feature1) - min(Feature1) );
% NormFeature2 =  (Feature2 - min(Feature2)) / ( max(Feature2) - min(Feature2) );
% NormFeature3 =  (Feature3 - min(Feature3)) / ( max(Feature3) - min(Feature3) );
% NormFeature4 =  (Feature4 - min(Feature4)) / ( max(Feature4) - min(Feature4) );
% NormFeature5 =  (Feature5 - min(Feature5)) / ( max(Feature5) - min(Feature5) );
NormFeature1 =  Normalize(Feature1, 25, 75);
NormFeature2 =  Normalize(Feature2, 25, 75);
NormFeature3 =  Normalize(Feature3, 25, 75);
NormFeature4 =  Normalize(Feature4, 25, 75);
NormFeature5 =  Normalize(Feature5, 25, 75);
MomentFeatures = [NormFeature1, NormFeature2, NormFeature3, NormFeature4, NormFeature5];


% Choosing cluster features and get k-dist for each point
% select_moments = [1,2];
X = MomentFeatures(:,select_moments);
MinPts = 2*size(X,2);

% Alternate approach to find epsilon is knn serach 
% and look for knee point for MinPts' nearest neighbors
% [MinPtsNN,PairwiseDistances] = knnsearch(X,X,'K',MinPts);

% Get clustering set for each choice of epsilon
% For each value in eps_range, store associated percent of data
eps_percent = zeros(length(eps_range),1);

% Iterate over different choices of epsilon
for eps = 1:length(eps_range)
    % Report the eps' set of clustering
    clusters_found = dbscan(X,eps_range(eps),MinPts);
    
    % Determine unique clusters found
    unique_clusters = unique(clusters_found);
   
    % Exclude the case in which the only cluster is the noise-cluster
    % and case where values are in one cluster
    if unique_clusters(end) == -1
        eps_percent(eps) = 0;
    elseif length(unique_clusters) == 1
        eps_percent(eps) = 100;
    else
        % Count clusters and find largest cluster
        [cluster_sizes, ~] = hist(clusters_found,unique_clusters);
        [max_size,max_idx] = max(cluster_sizes);
        
        % If the largest cluster is noise-cluster,
        % find second largest cluster
        if max_idx == 1
            
            cluster_sizes = sort(cluster_sizes, 'descend');
            max2_size = cluster_sizes(2);
            eps_percent(eps) = max2_size/length(jump_stats);
  
        else
            eps_percent(eps) = max_size/length(jump_stats);
        end
    end

   
end

figure;
scatter(eps_range,eps_percent, "filled")
set(gca, 'XScale', 'log');
ylabel('Fraction of points in largest cluster');
xlabel('Distance (\epsilon)');

final_epsilons = zeros(size(desired_fraction));
final_fracs = zeros(size(desired_fraction));

for frac_i = 1:length(desired_fraction)

    % Initialize parameters for epsilon search
    eps_search = eps_range(1);
    frac_search = eps_percent(1);
    
    % Search through eps_num to find point near desired fraction
    for i = 1:length(eps_percent)
        % If closer fraction is found, update final values
        if eps_percent(i) < desired_fraction(frac_i)
            frac_search = eps_percent(i);
            eps_search = eps_range(i);
        end
    end
    
    final_epsilons(frac_i) = eps_search;
    final_fracs(frac_i) = frac_search;
end

cmap = [orange;
  blue; ...
  green; ...
  cyan];

% Recolor plots based on multiple eps options

% Get clustering set for each choice of epsilon
final_clusters = zeros(length(jump_stats),1);
% Iterate over different choices of epsilon
for eps_i = 1:length(final_epsilons)
    cluster_i = dbscan(X,final_epsilons(eps_i), MinPts);

%   Find the biggest cluster for assignment
    [cnt_unique, unique_values] = hist(cluster_i,unique(cluster_i));
    % Find biggest cluster
    [~,max_idx] = max(cnt_unique);

    % If noise the biggest cluster, find second biggest cluster
    if max_idx == 1
        [~,max2_idx] = max(cnt_unique(2:end));
        cluster_label = unique_values(max2_idx+1);
    else 
        
        cluster_label = unique_values(max_idx);
    end
    
    
    for i = 1:length(jump_stats)
        if cluster_i(i) == cluster_label
            final_clusters(i) = eps_i;
        end
    end
end


set(groot,'defaultAxesColorOrder',cmap)



% Plot clusters in in feature space
figure;
% scatter(NormFeature1, NormFeature2, 15, final_clusters, "filled");
plot(NormFeature1(final_clusters==0),NormFeature2(final_clusters==0),'.','MarkerSize',14); hold on
plot(NormFeature1(final_clusters==1),NormFeature2(final_clusters==1),'.','MarkerSize',14); hold on
plot(NormFeature1(final_clusters==2),NormFeature2(final_clusters==2),'.','MarkerSize',14); hold on
% plot(NormFeature1(final_clusters==3),NormFeature2(final_clusters==3),'.','MarkerSize',14); hold on
xlabel('Standard deviation (normalized)');
ylabel('FWHM (normalized)');
colormap(cmap)

% Plot clusters in frequency space
figure;
plot(jumps_measured(final_clusters==0,5),jumps_measured(final_clusters==0,6),'.','MarkerSize',14); hold on
plot(jumps_measured(final_clusters==1,5),jumps_measured(final_clusters==1,6),'.','MarkerSize',14); hold on
plot(jumps_measured(final_clusters==2,5),jumps_measured(final_clusters==2,6),'.','MarkerSize',14); hold on
% plot(jumps_measured(final_clusters==3,2),jumps_measured(final_clusters==3,3),'.','MarkerSize',14); hold on
% scatter(jumps_measured(:,2), jumps_measured(:,3), 15, final_clusters, "filled");
xlabel('Relative frequency shift (Mode 1)');
ylabel('Relative frequency shift (Mode 2)');
% colormap(cmap)
