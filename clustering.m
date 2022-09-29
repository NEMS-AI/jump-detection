% Specify the desired fraction of data set as 
desired_fraction = .6;
% desired_fraction = .2;

% Specify range of epsilon values to sweep over
eps_range = linspace(0.000,.1,1000);

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
Feature1 = jump_stats(:,4); 
Feature2 = jump_stats(:,3);

NormFeature1 =  (Feature1 - min(Feature1)) / ( max(Feature1) - min(Feature1) );
NormFeature2 =  (Feature2 - min(Feature2)) / ( max(Feature2) - min(Feature2) );

% Choosing cluster features and get k-dist for each point
X = [NormFeature1 NormFeature2];
[Idx,D] = knnsearch(X,X,'K',4);

% Get clustering set for each choice of epsilon
final_clusters = zeros(length(Idx),1);
eps_num = zeros(length(eps_range),1);

% Iterate over different choices of epsilon
for eps = 1:length(eps_range)
    idx = dbscan(X,eps_range(eps),4);

    % Determine clusters dbscan found
    unique_clusters = unique(idx);
   
    % Exclude the the case in which the only cluster is the noise-cluster
    if unique_clusters(end) == -1
        eps_num(eps) = 0;

    else
        % Count clusters and find largest cluster
        [cnt_unique, unique_a] = hist(idx,unique_clusters);
        [M,I] = max(cnt_unique);
        
        % If the largest cluster is noise-cluster,
        % find second largest cluster
        if I == 1
            eps_num(eps) = max(cnt_unique(cnt_unique<max(cnt_unique)));
  
        else
            eps_num(eps) = M;
        end
    end
   
end

figure;
scatter(eps_range,eps_num/eps_num(end), "filled")
set(gca, 'XScale', 'log');
ylabel('Fraction of points in largest cluster');
xlabel('Distance (\epsilon)');

% Initialize parameters for epsilon search
final_eps = eps_range(1);
final_frac = eps_num(1)/eps_num(end);

% Search through eps_num to find point near desired fraction
for i = 1:length(eps_num)
    % If closer fraction is found, update final values
    if eps_num(i)/eps_num(end) < desired_fraction
        final_frac = eps_num(i)/eps_num(end);
        final_eps = eps_range(i);
    end
end

epsilon = final_eps;
pts = final_frac;

cmap = [orange;
  blue; ...
  green; ...
  cyan];

% Recolor plots based on multiple eps options

% Get clustering set for each choice of epsilon
final_clusters = zeros(length(Idx),1);
% Iterate over different choices of epsilon
for eps = 1:length(epsilon)
    idx = dbscan(X,epsilon(eps),4);

%   Find the biggest cluster for assignment
    [cnt_unique, unique_a] = hist(idx,unique(idx));
    [M,I] = max(cnt_unique);
    [M2,I2] = max(cnt_unique(cnt_unique<max(cnt_unique)));
    
    if unique_a(end) == -1
        continue;
    end
    for i = 1:length(idx)
        if I == 1
            if idx(i) == unique_a(I2+1)
                final_clusters(i) = eps;
            end
        elseif idx(i) == unique_a(I)
            final_clusters(i) = eps;
        end
    end
end


set(groot,'defaultAxesColorOrder',cmap)



% Plot clusters in in feature space
figure;
% scatter(NormFeature1, NormFeature2, 15, final_clusters, "filled");
plot(NormFeature1(final_clusters==0),NormFeature2(final_clusters==0),'.','MarkerSize',14); hold on
plot(NormFeature1(final_clusters==1),NormFeature2(final_clusters==1),'.','MarkerSize',14); hold on
% plot(NormFeature1(final_clusters==2),NormFeature2(final_clusters==2),'.','MarkerSize',14); hold on
% plot(NormFeature1(final_clusters==3),NormFeature2(final_clusters==3),'.','MarkerSize',14); hold on
xlabel('Standard deviation (normalized)');
ylabel('FWHM (normalized)');
colormap(cmap)

% Plot clusters in frequency space
figure;
plot(jumps_measured(final_clusters==0,5),jumps_measured(final_clusters==0,6),'.','MarkerSize',14); hold on
plot(jumps_measured(final_clusters==1,5),jumps_measured(final_clusters==1,6),'.','MarkerSize',14); hold on
% plot(jumps_measured(final_clusters==2,2),jumps_measured(final_clusters==2,3),'.','MarkerSize',14); hold on
% plot(jumps_measured(final_clusters==3,2),jumps_measured(final_clusters==3,3),'.','MarkerSize',14); hold on
% scatter(jumps_measured(:,2), jumps_measured(:,3), 15, final_clusters, "filled");
xlabel('Relative frequency shift (Mode 1)');
ylabel('Relative frequency shift (Mode 2)');
% colormap(cmap)


