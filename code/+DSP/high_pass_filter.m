function  filter_func = high_pass_filter(fs, f_stop, f_pass, options)
    arguments
        fs     (1,1) double
        f_stop (1,1) double
        f_pass (1,1) double        
        options.A_pass  (1,1) double = 0.5
        options.A_stop  (1,1) double = 40  
        options.method  (1,1) string {mustBeMember( options.method, [ ...
            "butter", "cheby1", "cheby2", "ellip", "equiripple", "firls", "freqsamp", "ifir", ...
            "iirlinphase", "iirlpnorm", "iirls", "fircls", "kaiserwin", "maxflat", "multistage", "window"] ...
        )} = "butter"
    end
    % Check inputs:
    assert(f_stop<f_pass)
    assert(f_pass<fs)   
    
    % Create Filter:
    A_pass = options.A_pass;
    A_stop = options.A_stop;
    d = fdesign.highpass('Fst,Fp,Ast,Ap',f_stop, f_pass, A_stop, A_pass, fs);
    lp_filt = design(d, options.method);

    % Create function handle:
    filter_func = @(sig) filter(lp_filt, sig);
end