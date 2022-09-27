function jumps_detected = preprocess(params, jumps_detected)
% Preprocess the detected jumps to exclude the ones
% that are likely noise outliers


% TODO move processing from measure function to here
ti = 1;
fprintf('Preprocessing Detected Jumps:         ');
while ti < params.tfin
    
    fprintf('\b\b\b\b\b\b\b%5.01f%%\n',ti/params.ttot*100); 
    ti = ti+1;
end
end
