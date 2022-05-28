function [X,f] = fft(x,Fs)
    arguments
        x  (:,1) double       
        Fs (1,1) double
    end
    
    N = length(x);
    % freqHz = (0:1:N-1)*Fs/N;

    [X, f] = my_fft(x, N,Fs);

end

function [X,f] = my_fft(x, n, Fs )

    arguments
        x    (:,1) double {mustBeNumeric, mustBeReal}  
        n    (1,1) double 
        Fs   (1,1) double 
    end
    % Constants:
    dim = 1;

    % fft + shift:
    X = fft(x,n,dim);
    X = fftshift(X,dim);
    
    % frequency vec: 
    N = size(X, dim); % number of samples
    delta_f = Fs / N;
    f = 0 : delta_f : delta_f*(N-1);
    f = fftshift(f);
    ind = f>=Fs/2-eps(Fs/2);
    f(ind) = f(ind)-Fs;
    f = shiftdim(f, 2-dim);
end
