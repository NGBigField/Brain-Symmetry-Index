function  filter_func = low_pass_filter(fs, f_pass, f_stop, options)
    arguments
        fs     (1,1) double
        f_pass (1,1) double
        f_stop (1,1) double
        options.A_pass  (1,1) double = 0.5
        options.A_stop  (1,1) double = 40          
        options.method  (1,1) string {mustBeMember( options.method, [ ...
            "butter", "cheby1", "cheby2", "ellip", "equiripple", "firls", "freqsamp", "ifir", ...
            "iirlinphase", "iirlpnorm", "iirls", "fircls", "kaiserwin", "maxflat", "multistage", "window"] ...
        )} = "butter"
    end
    % Check inputs:
    assert(f_pass<f_stop)
    assert(f_stop<fs)   
    
    % Create Filter:
    A_pass = options.A_pass;
    A_stop = options.A_stop;
    d = fdesign.lowpass('Fp,Fst,Ap,Ast',f_pass, f_stop, A_pass, A_stop, fs);
    lp_filt = design(d, options.method);

    % Create function handle:
    filter_func = @(sig) filter(lp_filt, sig);
end