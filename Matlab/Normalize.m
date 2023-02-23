function Normalized =  Normalize(Feature, lower_pct, higher_pct)
    lower = prctile(Feature,lower_pct);
    higher = prctile(Feature,higher_pct);
    Normalized = 2*(Feature - lower) / ( higher - lower )-1;
end