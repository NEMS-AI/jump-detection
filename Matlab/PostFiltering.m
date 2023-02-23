function PostFilter = PostFiltering(jumps_stats, percentile)
    thresh = prctile(jumps_stats(:,8),percentile);
    FWHMCondition = zeros(size(jumps_stats, 1),1);
    for i = 1:length(jumps_stats)
        if jumps_stats(i,8) > thresh
            FWHMCondition(i) = 1;
        end
    end
    PostFilter = FWHMCondition;
end