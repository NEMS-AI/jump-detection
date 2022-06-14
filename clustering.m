
epsilon_10xSNR = [.029 .0125 .0035];
epsilon_1xSNR = [.0548 .01871 .00781];
epsilon_GROEL = [.0239 .012 .0069];
epsilon = epsilon_GROEL;

% Color = [255 201 14]/255; % yellow
% Color1 = [0 0 0]/255;   % black
% Color4 = [39 150 235]/255;  % blue
Color1 = [61 38 168]/255;   % purple
Color2 = [217 83 25]/255;   % orange
Color3 = [34 177 76]/255;   % green
Color4 = [77 190 238]/255;  % cyan

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
% plot(sort(KDistValues),'Color',Color1,'LineWidth',2)
scatter(eps_range, eps_num/eps_num(end), "filled")
set(gca, 'XScale', 'log');
% xlim([0 length(KDistValues)]);
% ylim([1e-4 1]);
ylabel('Fraction of points in largest cluster');
xlabel('Distance (\epsilon)');
hold on;

% Add epsilon lines to k-dist graph
xl3 = xline(epsilon(3), '--','LineWidth',2);
xl3.Color = Color4;
hold on
xl2 = xline(epsilon(2), '--','LineWidth',2);
xl2.Color = Color3;
xl1 = xline(epsilon(1), '--','LineWidth',2);
xl1.Color = Color2;
% legend('All data','Knee point','Restricted','More restricted','Location','northwest');
legend('All data','Selection 1','Selection 2','Selection 3','Location','northwest');

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


