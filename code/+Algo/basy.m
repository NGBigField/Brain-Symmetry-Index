% Brain Asymmetry:

function [BAsyMat, freqVec] = basy( HHTsCell, dF, freqLimits )
    arguments
        HHTsCell    (:,2) cell
        dF          (1,1) double
        freqLimits  (1,2) double
    end
    %% Parse data:
    nIMFs = size(HHTsCell ,1);   

    %% Create common freq interval vectors:
    freqVec = freqLimits(1):dF:freqLimits(2);    
    freqVec = freqVec.';
    nFreqs = length(freqVec);

    %% Calc marginal hilbert spectrum: H -> to h:
    h_cell = cell(size(HHTsCell));
    for iIMF = 1 : nIMFs
        for iSide = 1 : 2
            t_vec = HHTsCell{iIMF, iSide}.t;
            f_vec = HHTsCell{iIMF, iSide}.f;
            H_vec = HHTsCell{iIMF, iSide}.H;
            [f, h] = Algo.marginal_hilbert( t_vec, f_vec, H_vec );
            h_cell{iIMF, iSide} = struct(f=f, h=h);
        end
    end

    %% Preallocate Output:
    BAsyMat = zeros(nFreqs-1 ,nIMFs);

    %% Iterate: 
    for iIMF = 1 : nIMFs
        for iF = 1 : (nFreqs-1)
            freq_interval = [freqVec(iF), freqVec(iF+1)];
            left  = calc_sum_marginal_squared(iIMF, 1, freq_interval, h_cell);
            right = calc_sum_marginal_squared(iIMF, 2, freq_interval, h_cell);
            BAsy = (left-right)/(left+right);
            BAsyMat(iF, iIMF) = BAsy;
        end
    end

    %% Adjust output:
    freqVec = freqVec(1:end-1);
end
%% end

%% Subs:
function h_sum = calc_sum_marginal_squared(iIMF, iSide, freq_interval, h_cell)    
    h = h_cell{iIMF, iSide}.h;
    f = h_cell{iIMF, iSide}.f;
    indices_in_interval = f>=freq_interval(1) & f<=freq_interval(2);
    h_sum = sum(h(indices_in_interval).^2);
end