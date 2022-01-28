function [tvect,fvect] = load_data(filenames, data_type, nmodes)
% Loads data into time and frequency vectors
% Examples provided here for loading experimental and synthetic data based
% on data_selection
% Input:
%   filenames: string array with files containing data
%   data_type: selection from noise, data, synthetic
%   nmodes: number of modes
%   header_lines: number of header lines in data files, if applicable
% Output:
%   tvect: time series vector in seconds
%   fvect: frequency vector, each row is a different mode.

% different formats of data files were used based on # of modes
if nmodes == 2
    colformat = '%f %f %f %f %f';
    n_header_lines = 82;
elseif nmodes == 3
    colformat = '%f %f %f %f %f %f';
    n_header_lines = 1;
end

fvect = [];
if strcmp(data_type, 'noise') || strcmp(data_type, 'data')
    
    for fi=1:nmodes   
        fid = fopen(filenames(fi));
        nems_data = textscan(fid, colformat,'HeaderLines', n_header_lines );
        fclose(fid);
        if fi == 1
            tvect = nems_data{1};
            tvect=tvect-tvect(1);
        end
        
        freq_i = nems_data{4};
        % optionally scale noise for simulation purpose
        if strcmp(data_type, 'noise')
            noise_factor = 1;
            freq_i_start = freq_i(1);
            freq_i = (freq_i-freq_i_start)*noise_factor + freq_i_start;
        end
        fvect = [fvect; freq_i'];
    end
end

