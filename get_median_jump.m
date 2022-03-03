function [med_rel_jump,med_Fstats] = get_median_jump(usejumps,rel_jump_ts_1,rel_jump_ts_2,rel_jump_ts_3,Fstats_ts)
%get_median_jump finds the median of the time series selected with usejumps
% if one of the modes isn't used, pass []
    med_rel_jump = [median(rel_jump_ts_1(usejumps,:)); median(rel_jump_ts_2(usejumps,:))];
    if ~isempty(rel_jump_ts_3)
        med_rel_jump = [med_rel_jump; median(rel_jump_ts_3(usejumps,:))];
    end
    med_Fstats = median(Fstats_ts(usejumps,:));
end

