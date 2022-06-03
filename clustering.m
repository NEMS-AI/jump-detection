% Parameter selection for colors
epsilon10xSNR = 0.055;
epsilon1xSNR = 0.045;
epsilonGROEL = .038;

epsilon = epsilon10xSNR;
epsilon_scale = [1, 1/3, 1/6];

Color1 = [0 0 0];
Color2 = [1 0 0];
Color3 = [0 1 0];
Color4 = [0 0 1];

cmap = [Color1;
  Color2; ...
  Color3; ...
  Color4];

% X axis = Feature 1, Y Axis = Feature 2
Feature1 = jump_stats(:,6); 
Feature2 = jump_stats(:,3);

NormFeature1 =  (Feature1 - min(Feature1)) / ( max(Feature1) - min(Feature1) );
NormFeature2 =  (Feature2 - min(Feature2)) / ( max(Feature2) - min(Feature2) );

% Choosing cluster features and get k-dist for each point
X = [NormFeature1 NormFeature2];
[Idx,D] = knnsearch(X,X,'K',10);

% Plot k-dist graph
figure;
KDistValues =  D(:,4);
plot(sort(KDistValues), Color=Color1)
hold on;

% Add epsilon lines to k-dist graph
yl1 = yline(epsilon*epsilon_scale(1), '--' ,'LineWidth',3);
yl1.Color = Color2;
hold on;
yl2 = yline(epsilon*epsilon_scale(2), '--' ,'LineWidth',3);
yl2.Color = Color3;
hold on;
yl3 = yline(epsilon*epsilon_scale(3), '--' ,'LineWidth',3);
yl3.Color = Color4;

% Get clustering set for each choice of epsilon
final_clusters = zeros(length(Idx),1);
% Iterate over different choices of epsilon
for eps = 1:length(epsilon_scale)
    idx = dbscan(X,epsilon*epsilon_scale(eps),4);

%   Find the biggest cluster for assignment
    [cnt_unique, unique_a] = hist(idx,unique(idx));
    [M,I] = max(cnt_unique);
    
    for i = 1:length(idx)
        if unique_a(I)>0 && idx(i) == unique_a(I)
            final_clusters(i) = eps;
        end
    
    end

end




% Plot clusters in in feature space
figure;
scatter(NormFeature1, NormFeature2, 50, final_clusters, "filled");
colormap(cmap)

% Plot clusters in frequency space
figure;
scatter(jumps_measured(:,2), jumps_measured(:,3), 50, final_clusters, "filled");
colormap(cmap)


