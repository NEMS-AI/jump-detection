
% epsilon_10xSNR = [.02903 .012614 .003604];
% epsilon_1xSNR = [.05486 .01871 .00781];
% epsilon_GROEL = [.02393 .01202 .006907];
epsilon_10xSNR = [.019 .0037];
epsilon_1xSNR = [.0146 .007];
epsilon_GROEL = [.0142 .0072];
epsilon = epsilon_GROEL(2);

% pts_10xSNR = [.673 .578 .297];
% pts_1xSNR = [.804 .673 .465];
% pts_GROEL = [.728 .5422 .2804];
pts_10xSNR = [.62 .31];
pts_1xSNR = [.62 .31];
pts_GROEL = [.62 .31];

pts = pts_GROEL(2);

% Color = [255 201 14]/255; % yellow
% Color1 = [0 0 0]/255;   % black
Color2 = [39 150 235]/255;  % blue
% Color2 = [61 38 168]/255;   % purple
Color1 = [217 83 25]/255;   % orange
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
% figure;
% KDistValues =  D(:,4);
% plot(eps_range, eps_num/eps_num(end), '.','MarkerSize',14,'Color',Color1);
% set(gca, 'XScale', 'log');
% xlims = xlim;
% hold on;
% plot([epsilon(1) epsilon(1)],[0 pts(1)],'--','LineWidth',2,'Color',Color2);
% plot([epsilon(2) epsilon(2)],[0 pts(2)],'--','LineWidth',2,'Color',Color3);
% plot([epsilon(3) epsilon(3)],[0 pts(3)],'--','LineWidth',2,'Color',Color4);
% plot([xlims(1) epsilon(1)],[pts(1) pts(1)],'--','LineWidth',2,'Color',Color2);
% plot([xlims(1) epsilon(2)],[pts(2) pts(2)],'--','LineWidth',2,'Color',Color3);
% plot([xlims(1) epsilon(3)],[pts(3) pts(3)],'--','LineWidth',2,'Color',Color4);
% 
% use_eps1 = eps_range < epsilon(1);
% use_eps2 = eps_range < epsilon(2);
% use_eps3 = eps_range < epsilon(3);
% plot(eps_range(use_eps1), eps_num(use_eps1)/eps_num(end), '.','MarkerSize',14,'Color',Color2);
% plot(eps_range(use_eps2), eps_num(use_eps2)/eps_num(end), '.','MarkerSize',14,'Color',Color3);
% plot(eps_range(use_eps3), eps_num(use_eps3)/eps_num(end), '.','MarkerSize',14,'Color',Color4);
% ylabel('Fraction of points in largest cluster');
% xlabel('Distance (\epsilon)');
% 
% legend('All data','Restricted','More restricted','Most restricted','Location','northwest');

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
plot(jumps_measured(final_clusters==0,2),jumps_measured(final_clusters==0,3),'.','MarkerSize',14); hold on
plot(jumps_measured(final_clusters==1,2),jumps_measured(final_clusters==1,3),'.','MarkerSize',14); hold on
% plot(jumps_measured(final_clusters==2,2),jumps_measured(final_clusters==2,3),'.','MarkerSize',14); hold on
% plot(jumps_measured(final_clusters==3,2),jumps_measured(final_clusters==3,3),'.','MarkerSize',14); hold on
% scatter(jumps_measured(:,2), jumps_measured(:,3), 15, final_clusters, "filled");
xlabel('Relative frequency shift (Mode 1)');
ylabel('Relative frequency shift (Mode 2)');
% colormap(cmap)


