function [sorted_freqs, h] = marginal_hilbert( tVec, fVec, HVec )
    sorted_freqs = sort( unique(fVec) ); 
    nFreqs = length(sorted_freqs);
    h = zeros(nFreqs, 1);
    for iF = 1 : nFreqs 
        f = sorted_freqs(iF);
        h(iF) = sum( HVec(fVec==f) );
    end
end