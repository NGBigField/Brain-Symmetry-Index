function [t, f, h] = hht(imf, fs, options)
    arguments
        imf
        fs
        options.freqResolution = 0.01
        options.freqLimits = [0,40]
    end
    % hht analysis:
    [hs, freq_used, time_used] = hht(imf, fs, FrequencyResolution=options.freqResolution, FrequencyLimits=options.freqLimits );

    % Find elements:
    [non_zero_vec_i,  non_zero_vec_j, h] = find(hs);
    f = freq_used(non_zero_vec_i);
    t = time_used(non_zero_vec_j);
end