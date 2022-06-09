% Parameter selection for colors
epsilon = 1;

eps_range = linspace(0.000,.1,1000);

Color1 = [61 38 168]/255;
% Color1 = [0 0 0];
Color2 = [39 150 235]/255;
% Color3 = [128 203 88]/255;
Color3 = [34 177 76]/255;
Color4 = [255 201 14]/255;
% Color4 = [255 127 39]/255;

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


% Get clustering set for each choice of epsilon
final_clusters = zeros(length(Idx),1);
eps_num = zeros(length(eps_range),1);
% Iterate over different choices of epsilon
for eps = 1:length(eps_range)
    idx = dbscan(X,epsilon*eps_range(eps),4);

%   Find the biggest cluster for assignment
    unique_clusters = unique(idx);
   
    if unique_clusters(end) == -1
        eps_num(eps) = 0;
    else
        [cnt_unique, unique_a] = hist(idx,unique_clusters);
        [M,I] = max(cnt_unique);
        
        if I == 1
            eps_num(eps) = max(cnt_unique(cnt_unique<max(cnt_unique)));
  
        else
            eps_num(eps) = M;
        end
    end
    
end

scatter(eps_range, eps_num, "filled")
xlabel('Epsilon');
ylabel('# Points in largest cluster');




