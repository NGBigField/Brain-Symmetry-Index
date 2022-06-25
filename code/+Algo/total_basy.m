function [over_all_sum, per_imf] = total_basy(BAsyMat, BAsyFreqs, options)
    arguments
        BAsyMat                 (:,:) double
        BAsyFreqs               (:,1) double
        options.freqsToInclude  (1,2) double          = [1, 25]  % Frequencies to take into account during calculation.
        options.imfsToInclude   (1,:) {mustBeInteger} = 1 : 100  % IMFs        to take into account during calculation.
    end
    
    %% Check inputs:
    assert(size(BAsyMat,1)==length(BAsyFreqs));
    givenIMFs = 1:size(BAsyMat,2);

    % derive Wanted freq indices:
    indices_within_freqs = BAsyFreqs>=options.freqsToInclude(1) & BAsyFreqs<=options.freqsToInclude(2);

    % init outputs:    
    over_all_sum = 0;
    per_imf = zeros(size(givenIMFs));

    % Iterate
    for iImf = givenIMFs    
        % calc:
        BAsyVec = BAsyMat(:,iImf);        
        BAsy_vals_within_freqs = BAsyVec( indices_within_freqs );
        imf_sum = sum(BAsy_vals_within_freqs, "all", "omitnan");
       
        % Keep data:
        per_imf(iImf) = imf_sum;        
        if ~ismember(iImf, options.imfsToInclude) 
            continue % Filter from overall sum:
        end
        over_all_sum = over_all_sum + imf_sum;
    end
end