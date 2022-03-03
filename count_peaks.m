function npeaks = count_peaks(Fstats,ti1,ti2,Fdetect,Fmeas,hi_thresh,lo_thresh)
%count_peaks counts the number of peaks in Fstats between ti1 and ti2 based
% on given threshold criteria
%   Fstats: F statistic as a function of time for the entire dataset
%   ti1: time index where F stat crosses above detection threshold
%   ti2: time index where F stat cross back below detection threshold
%   Fdetect: detection threshold
%   Fmeas: measurement threshold (only count peaks that cross above this).
%   hi_thresh: fraction of hi, such that crossing below counts a peak.
%   lo_thresh: fraction of lo, such that crossing above starts looking for
%   a new peak.

    ti = ti1;
    Fstathi = Fstats(ti);
    Fstatlo = Fdetect*hi_thresh;
    npeaks = 0;
    peakfound = 0;
    
    tii = 0;
    while ti+tii < ti2 && tii < length(Fstats)
        tii = tii+1;
        if peakfound == 0
            if Fstats(ti+tii) > Fstathi
                Fstathi = Fstats(ti+tii);
            elseif Fstats(ti+tii) < Fstathi*hi_thresh
                npeaks = npeaks + 1;
                Fstatlo = Fstats(ti+tii);
                peakfound = 1;
            end
        elseif peakfound == 1
            if Fstats(ti+tii) < Fstatlo
                Fstatlo = Fstats(ti+tii);
            elseif Fstats(ti+tii) > Fstatlo*lo_thresh && Fstats(ti+tii) > Fmeas
                Fstathi = Fstats(ti+tii);
                peakfound = 0;
            end
        end 
    end
    if Fstats(ti2) > Fmeas && peakfound == 0
        npeaks = npeaks + 1;
    end    
end

