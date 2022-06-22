function y = band_pass_filter(x, fs)
    f_pass1 = 2;
    f_pass2 = 45;
    high_pass_filt = @(x) highpass(x, f_pass1, fs);
    low_pass_filt = @(x) lowpass(x, f_pass2, fs);
    band_pass_filt_func = @(x) high_pass_filt(low_pass_filt(x));    
    y = band_pass_filt_func(x);
end