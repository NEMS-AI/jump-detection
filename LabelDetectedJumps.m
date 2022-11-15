function classes = LabelDetectedJumps(detected_data, real_times, tmeas, tjump)
    % Every jump is either, fp(0), single(1) or multi jump(2)
    classes = zeros(size(detected_data(:,1)));
    % Iterate through detected jumps and assign them as either a false positive
    % or true negative
    
    for measured_event_i = 1:size(detected_data,1)
    
        % Build a range around each detected jump
        detected_event = detected_data(measured_event_i,1);
        lower_range = detected_data(measured_event_i,2);
        upper_range = detected_data(measured_event_i,3);
        
        % Iterate through real events and see if any fall within jump range
        for real_event_i=1:length(real_times)
            real_event = real_times(real_event_i);
            if lower_range <= real_event && real_event <= upper_range
                classes(measured_event_i) = 1;
                real_count = 0;
                % If jump found, check surounding area for >1 jumps
                for real_event_i=1:length(real_times)
                    real_event = real_times(real_event_i);
                    if lower_range-tmeas-tjump <= real_event && real_event <= upper_range+tmeas+2*tjump
                        real_count = real_count + 1;
                        if real_count > 1
                            classes(measured_event_i) = 2;
                        end
                    end
                end
            end
    
        end
    
    end
end
