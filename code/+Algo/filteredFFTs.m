function [L, R, f] = filteredFFTs(l, r, options)
%function [L, R, f] = filteredFFTs(l, r, options) used for Brain Symmetry Index:
% A code based on:
% Continuous Quantitative EEG Monitoring in Hemispheric Stroke Patients Using the Brain Symmetry Index
    
    arguments
        l                (:,1) Classes.EEGSignal  % EEG of left side
        r                (:,1) Classes.EEGSignal  % EEG of right side
        options.isPlot   (1,1) logical            = false
        options.reqFreqs (1,2) double             = [1, 25] % [Hz]         
    end

    %% Check inputs:
    assert(l.placement.side==Classes.Side.Left )
    assert(r.placement.side==Classes.Side.Right)
    assert(l.props.frequency==r.props.frequency)
    Fs = l.props.frequency;
       
    %% FFT:
    [L, fL] = Algo.fft(l.recording, Fs);
    [R, fR] = Algo.fft(r.recording, Fs); 

    %% Filter Requested Frequencies:
    [L, fL] = filter_frequencies(L, fL, options.reqFreqs);
    [R, fR] = filter_frequencies(R, fR, options.reqFreqs);
    assert(all(fL==fR));
    f = fL;

    %% Plot:
    if options.isPlot
        figure;
        cellItems = {{L,"Left"},{R,"Right"}};
        for iC = 1:length(cellItems)
            cellItem = cellItems{iC};
            vec = cellItem{1};
            side = cellItem{2};
            hold on;
            plotH = plot(f, abs(vec));
            plotH.LineWidth = 1;
            plotH.DisplayName = side;
        end
        legend
        grid on; grid minor
        xlabel("freq [Hz]")
        ylabel("fft ")        
    end

    %% Compute:


    
end

%% End

%% Sub-Funcions:
function [X, f] = filter_frequencies(X, f, reqFreqs)
    i1 = find_closest_index(f, reqFreqs(1));
    i2 = find_closest_index(f, reqFreqs(2));
    X = X(i1:i2);
    f = f(i1:i2);
end

function res = find_closest_index(vec, val)
    arguments
        vec (:,1) double % Must be monotone growing
        val (1,1) double
    end
    distance = abs(vec-val);
    [~, minInd] = min(distance);
    res = minInd;
end