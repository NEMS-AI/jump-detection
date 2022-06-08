% Parameter selection for colors
epsilon10xSNR = 0.055;
epsilon1xSNR = 0.045;
epsilonGROEL = .038;
epsilon = epsilonGROEL;

% epsilon_scale = [1, 1/3, 1/6];
epsilon_scale_10xSNR = [1 1/4 .14];
epsilon_scale_1xSNR = [1 .35 .2];
epsilon_scale_GROEL = [1 .5 .28];
epsilon_scale = epsilon_scale_GROEL;

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

% Plot k-dist graph
figure;
KDistValues =  D(:,4);
plot(sort(KDistValues),'Color',Color1,'LineWidth',2)
set(gca, 'YScale', 'log');
xlim([0 length(KDistValues)]);
ylim([1e-4 1]);
xlabel('Points');
ylabel('Distance (\epsilon)');
hold on;

% Add epsilon lines to k-dist graph
yl1 = yline(epsilon*epsilon_scale(1), '--','LineWidth',2);
yl1.Color = Color2;
hold on;
yl2 = yline(epsilon*epsilon_scale(2), '--','LineWidth',2);
yl2.Color = Color3;
hold on;
yl3 = yline(epsilon*epsilon_scale(3), '--','LineWidth',2);
yl3.Color = Color4;
legend('All data','Knee point','Restricted','More restricted','Location','northwest');

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


set(groot,'defaultAxesColorOrder',cmap)


% Plot clusters in in feature space
figure;
% scatter(NormFeature1, NormFeature2, 15, final_clusters, "filled");
plot(NormFeature1(final_clusters==0),NormFeature2(final_clusters==0),'.','MarkerSize',14); hold on
plot(NormFeature1(final_clusters==1),NormFeature2(final_clusters==1),'.','MarkerSize',14); hold on
plot(NormFeature1(final_clusters==2),NormFeature2(final_clusters==2),'.','MarkerSize',14); hold on
plot(NormFeature1(final_clusters==3),NormFeature2(final_clusters==3),'.','MarkerSize',14); hold on
xlabel('Standard deviation (normalized)');
ylabel('FWHM (normalized)');
colormap(cmap)

% Plot clusters in frequency space
figure;
plot(jumps_measured(final_clusters==0,2),jumps_measured(final_clusters==0,3),'.','MarkerSize',14); hold on
plot(jumps_measured(final_clusters==1,2),jumps_measured(final_clusters==1,3),'.','MarkerSize',14); hold on
plot(jumps_measured(final_clusters==2,2),jumps_measured(final_clusters==2,3),'.','MarkerSize',14); hold on
plot(jumps_measured(final_clusters==3,2),jumps_measured(final_clusters==3,3),'.','MarkerSize',14); hold on
% scatter(jumps_measured(:,2), jumps_measured(:,3), 15, final_clusters, "filled");
xlabel('Relative frequency shift (Mode 1)');
ylabel('Relative frequency shift (Mode 2)');
% colormap(cmap)


