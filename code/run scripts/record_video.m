% How-To:
%
% Before running this script, have one tagged-Data loaded to your workspace.
% These can be found in folder: Data\tagged\
%
%% clear
clearvars('-except','Data');
addpath(genpath(pwd));
close all; clc;

% Assert user chose Data:
error_msg = "Before running this script, have one tagged-Data loaded to your workspace. These can be found in folder: Data\tagged\";
assert(exist("Data", "var"), error_msg)

%% Constants:
% Choose if to filter:
Constants.use_pre_bandpass_filter = true;
Constants.freq_limits = [0.01, 45];
% How to calculate brain Assymetry:
Constants.delta_freq = 0.2;  % Brain Assymetry frequency resolution:
Constants.BAsy_freqs = Classes.BrainWavesFreqs.Alpha.freqs;  % freqs to indclude in BAsy clac
Constants.BAsy_imfs  = [2,3,4]; % imfs to indclude in BAsy clac

% Define what not to plot:
filterOutData.placements = []; %["F7", "F8"];
filterOutData.conditions = [];

%% Prepare figure, axes and visual helpers:
% Prepare figures:
[figH, subplotHandles] = init_subplots(Data, Constants);
% My Fig Recorder:
figRecorder = Classes.FigureRecorder(); % Prepare Recorder

%% Iterate:
nData = length(Data.EEGs);
progBar = Classes.ProgressBar(nData, cancelable=true);

for iData = 1 : nData
    %% Step:
    disp("Overall Progress: "+StringUtils.num_out_of_num(iData,nData))
    cancel_requested = progBar.step();
    if cancel_requested
        break
    end
    
    %% Skip?
    eegData = Data.EEGs(iData);
    isSkip = filter_data(eegData, filterOutData);
    if isSkip        
        continue
    end 

    try 
        %% Plot:
        plot_iteration(iData, Data, Constants, subplotHandles)

        %% Capture video:
        drawnow()
        figRecorder.capture(figH);
    catch ME
        error_msg = string(ME.getReport());
        warning(error_msg)
    end

end

%% Close:
progBar.close();
figRecorder.finish();
Sounds.finish()
disp("End.")

%% End:

% === % === % == % === % === % == % === % === % == % === % === % == % === % === % == % === % === % == % === % === % == % === % === % == %

% === % === % == % === % === % == % === % === % == % === % === % == % === % === % == % === % === % == % === % === % == % === % === % == %

%% Subs:
function line_props(plotH, side, color) 
    arguments
        plotH (1,1) matlab.graphics.chart.primitive.Line
        side  (1,1) string
        color 
    end 
    plotH.DisplayName = side;
    plotH.Color = color;
    plotH.LineWidth = 1;
end
%%
function pretty_plot(axisH)
    legH = legend(axisH);
    legH.Interpreter = 'none';
    legH.FontSize = 12;
    % grid(axisH,"on");
end
%% 
function [Limits] = derive_plot_limits(EEGs, Constants)
    %% Init Global limits:
    Limits.time = [inf, -inf];
    Limits.freq = [inf, -inf];
    Limits.bsi  = [inf, -inf];
    Limits.hht  = [inf, -inf];
    Limits.BAsy = [inf, -inf];

    % Init crnt Limits:
    crntLims = struct;
    fieldNames = string(fieldnames(Limits));
    for iF = 1 : length(fieldNames)
        fieldName = fieldNames(iF);
        crntLims.(fieldName) = struct;
        crntLims.(fieldName).min = [];
        crntLims.(fieldName).max = [];
    end

    %% Iterate:
    nData = length(EEGs);
    progBar = Classes.PrintedProgressBar(nData, msg="Deriving Plot Limits:");
    for i = 1 : nData
        progBar.step();
        
        % get data:
        eegData = EEGs(i);
        left  = eegData.left;
        right = eegData.right;
        [L, R, ~] = Algo.filteredFFTs(left, right, reqFreqs=Constants.freq_limits);      
        l = left.recording;
        r = right.recording;
        bsi = Algo.singleBSI(L,R);
        hht = [0, 30];
        BAsy = [-2, 2];        
        %if Constants.use_pre_bandpass_filter
        %    fs = eegData.left.props.frequency;
        %    l = DSP.band_pass_filter(l, fs);
        %    r = DSP.band_pass_filter(r, fs);
        %end

        % Derive current limits:
        lr = [l; r];
        LR = [abs(L); abs(R)];
        crntLims.time.max = max(lr);
        crntLims.time.min = min(lr);
        crntLims.freq.max = max(LR);
        crntLims.freq.min = min(LR);
        crntLims.bsi.max  = bsi;
        crntLims.bsi.min  = bsi;
        crntLims.hht.max  = max(hht);
        crntLims.hht.min  = min(hht);
        crntLims.BAsy.max = max(BAsy);
        crntLims.BAsy.min = min(BAsy);

        % Update global limits:
        cellItems = {{"min", @min, 1}, {"max", @max, 2}};
        for iF = 1 : length(fieldNames)
            valFieldName = fieldNames(iF);

            for iC =  1 : length(cellItems)
                str    = cellItems{iC}{1};
                minmax = cellItems{iC}{2};
                index  = cellItems{iC}{3};
                Limits.(valFieldName)(index) = minmax( Limits.(valFieldName)(index) , crntLims.(valFieldName).(str) );
            end
        end % for iF

    end % for i 
    progBar.clear();
end % func
%%
function  isSkip = filter_data(crntEegData, filterLists)
    arguments
        crntEegData (1,1) Classes.EEGData
        filterLists (1,1) struct
    end
    % derive crnt:
    pLeft     = crntEegData.left.placement.name;
    pRight    = crntEegData.right.placement.name;
    consition = crntEegData.condition;

    % Decide:
    isSkip = false;
    if IndexUtils.is_member(pLeft, filterLists.placements) 
        isSkip = true;
    end 
    if IndexUtils.is_member(pRight, filterLists.placements)
        isSkip = true;
    end
    if IndexUtils.is_member(consition, filterLists.conditions)
        isSkip = true;
    end 

    % Print skipped:
    if isSkip
        disp("Skipped eeg data with condition '"+string(consition)+"' and placements ["+string(pLeft)+", "+string(pRight)+"]")
    end

end
%% 
function [sig, color, side] = deal_data(i, Sigs, cyclic_colors, Sides)
    sig = Sigs(:,i);
    side = Sides(i);
    color = cyclic_colors(i);
end
%%
function [imf_num, color, name] = deal_imf_props_func(iImf, IMFs_to_include, cyclic_colors)
    imf_num = IMFs_to_include(iImf);
    color = cyclic_colors(iImf);
    name = "imf-"+string(imf_num); 
end
%%
function indices = get_IMFs_to_include(N)
    arguments
        N (1,1) {mustBeInteger}
    end
    % constants:
    minimum_req_index = 1;

    % Check consition:
    if N < minimum_req_index
        error("Not enough IMFs");
    end    

    % calc
    first = minimum_req_index;
    last = max(first, N);
    indices = first : last ;
end
%%
function child = get_child_with_name(parent, name)
    children = parent.Children;
    child = [];
    for iC = 1 : length(children)
        crntChild = children(iC);
        if string(crntChild.DisplayName) == name
            child = crntChild;
            return
        end
    end    
end
%%
function Limits = get_plot_limits(Data, Constants)
    Limits = derive_plot_limits(Data.EEGs, Constants);
    Limits.time = Limits.time/2;
    Limits.freq = Limits.freq/4;
    Limits.freq(1) = 0;
    Limits.bsi(1) = 0;
    Limits.bsi(2) = 2;
    Limits.hht = [0, 20];
    Limits.BAsy = [0, 20];
    Limits.BAsyPlot.y = Constants.BAsy_freqs;
    Limits.BAsyPlot.x = [-2, 2];
end
%%
function [figH, subplotHandles] = init_subplots(Data, Constants)
    % constants:
    Sides = ["left", "right"];

    % Adjust plot limits:
    Limits = get_plot_limits(Data, Constants);

    % Create Figure:
    figH = figure;
    figH.Position(1) = 5;
    figH.Position(2) = 65;
    figH.Position(3) = 2.3*figH.Position(3); % longer figure;
    figH.Position(4) = 1.5*figH.Position(4); % longer figure;

    % arange locations in grid:
    m=2; n=17;
    ind = @(x,y) Visuals.grid_xy_to_index(m,n,x,y);
    subplotHandles = struct();
    subplotHandles.time     = subplot(m,n,ind(  1   ,  1:4  ) );
    subplotHandles.freq     = subplot(m,n,ind(  2   ,  1:4  ) );
    subplotHandles.hhtLeft  = subplot(m,n,ind(  1   ,  6:9  ) );
    subplotHandles.hhtRight = subplot(m,n,ind(  1   , 11:14 ) );
    subplotHandles.bsi      = subplot(m,n,ind( [1,2], 16    ) );
    subplotHandles.BAsyBar  = subplot(m,n,ind( [1,2], 17    ) );
    subplotHandles.BAsyPlot = subplot(m,n,ind(  2   , 6:7   ) );
    subplotHandles.BAsyAccu = subplot(m,n,ind(  2   ,  9:14 ) ); 
    subplotHandles_hht_array = [subplotHandles.hhtLeft, subplotHandles.hhtRight];

    % Shift subplots' size and positions:
    enlarge_factor = 1.3;
    reduce_factor  = 0.7;
    left_factor    = 0.06;
    Visuals.sequencialy_change_axes_widths( [subplotHandles.freq, subplotHandles.BAsyPlot, subplotHandles.BAsyAccu], enlarge_factor );
    Visuals.sequencialy_change_axes_widths( [subplotHandles.time, subplotHandles.hhtLeft, subplotHandles.hhtRight, subplotHandles.bsi, subplotHandles.BAsyBar], enlarge_factor );
    Visuals.sequencialy_change_axes_widths( [subplotHandles.bsi, subplotHandles.BAsyBar] , reduce_factor );
    keys = fieldnames(subplotHandles);
    for iK = 1 : length(keys)
        key = string(keys(iK));
        subplotHandles.(key).Position(1) = subplotHandles.(key).Position(1) - left_factor;
    end
    subplotHandles.hhtLeft.Position(4)  = 0.9  * subplotHandles.hhtLeft.Position(4);
    subplotHandles.hhtRight.Position(4) = 0.9  * subplotHandles.hhtRight.Position(4);
    subplotHandles.BAsyPlot.Position(3) = 0.95 * subplotHandles.BAsyPlot.Position(3);
    subplotHandles.bsi.Position(1) = subplotHandles.bsi.Position(1) - 0.02;

    % time:
    axisH = subplotHandles.time;
    pretty_plot(axisH)
    xlabel(axisH, "time [sec]",FontSize=15, Interpreter="none")
    ylabel(axisH, "$$x(t)$$"  ,FontSize=16, Interpreter="latex")
    title(axisH, "Time Analysis");
    ylim(axisH, Limits.time)

    % freq:
    axisH = subplotHandles.freq;
    pretty_plot(axisH)
    xlabel(axisH, "freq [Hz]", FontSize=15, Interpreter="none")
    ylabel(axisH, "$$|X^{F}(f)|$$"  ,FontSize=16, Interpreter="latex")
    title(axisH, "Frequency Analysis");
    ylim(axisH, Limits.freq)

    % hht:
    for iSide = 1 : 2
        axisH = subplotHandles_hht_array(iSide);
        side = Sides(iSide);
        pretty_plot(axisH)
        xlabel(axisH, "time [sec]", FontSize=15, Interpreter="none")
        ylabel(axisH, "freq [hz] ", FontSize=15, Interpreter="none")
        title(axisH, "HHT "+side);
        ylim(axisH, Limits.hht)
        axisH.Legend.Position(2) = axisH.Legend.Position(2) + 0.03;
        axisH.Legend.Position(1) = axisH.Legend.Position(1) + 0.03;
        axisH.Legend.FontSize = 10;
    end

    % bsi:
    axisH = subplotHandles.bsi;
    xlim( axisH, [0,1/2]);
    ylim( axisH, Limits.bsi);
    title(axisH, "BSI")

    % BAsyBar:
    axisH = subplotHandles.BAsyBar;
    xlim( axisH, [0,1/2]);
    ylim( axisH, Limits.BAsy);
    title(axisH, "BAsy")

    % BAsyPlot:
    axisH = subplotHandles.BAsyPlot;    
    xlabel(axisH, "BAsy(freq)", FontSize=15, Interpreter="none")
    ylabel(axisH, "freq [Hz]" , FontSize=15, Interpreter="none")    
    title(axisH, "IMFs' Brain Asymmetry");
    xlim(axisH, Limits.BAsyPlot.x)
    ylim(axisH, Limits.BAsyPlot.y)
    lineH = xline(axisH, 0, Color='black', Tag="keep");
    lineH.Annotation.LegendInformation.IconDisplayStyle = "off";

    % BAsyAccu:
    axisH = subplotHandles.BAsyAccu;
    pretty_plot(axisH)
    xlabel(axisH, "samples", FontSize=15, Interpreter="none")
    ylabel(axisH, "$$\int BAsy$$" , FontSize=15, Interpreter="latex")    
    title(axisH, "Accumulated Brain Asymmetry");    
    lineH = yline(axisH, 0, Color='black', Tag="keep");
    lineH.Annotation.LegendInformation.IconDisplayStyle = "off";
    axisH.Legend.Location = "southwest";
    axisH.Legend.Position(2) = axisH.Legend.Position(2) - 0.02;
    axisH.Legend.Position(1) = axisH.Legend.Position(1) - 0.02;
    axisH.Legend.FontSize = 10;
end
%%
function plot_iteration(iData, Data, Constants, subplotHandles)

    % Constants:
    Sides = ["left", "right"];
    delta_freq              = Constants.delta_freq;
    freq_limits             = Constants.freq_limits;
    freq_to_include_in_BAsy = Constants.BAsy_freqs;
    IMFs_to_include_in_BAsy = Constants.BAsy_imfs;
    use_pre_bandpass_filter = Constants.use_pre_bandpass_filter;

    % Parse Basic data:
    eegData = Data.EEGs(iData);    
    cyclic_colors = @(i) Visuals.get_default_color_cyclic(i);
    subplotHandles_hht_array = [subplotHandles.hhtLeft, subplotHandles.hhtRight];

    %% parse signal data:    
    left  = eegData.left;
    right = eegData.right;
    time  = eegData.time;
    fs    = eegData.left.props.frequency;
    Sigs  = [left.recording, right.recording];

    %% Filter Signal:
    if use_pre_bandpass_filter
        Sigs(:,1) = DSP.band_pass_filter(Sigs(:,1), fs);
        Sigs(:,2) = DSP.band_pass_filter(Sigs(:,2), fs);
    end

    %% Title with global info:
    super_title_str = string(eegData.condition)+newline+left.placement.name+"-"+right.placement.name; 
    sgtitle(super_title_str);
    
    %% Plot sigs in time:
    axisH = subplotHandles.time;
    Visuals.clear_old_lines(axisH)
    for iSig = [1, 2]
        [sig, color, side] = deal_data(iSig, Sigs, cyclic_colors, Sides);
        % plot:
        hold( axisH, "on");
        plotH = plot( axisH, time, sig);
        line_props(plotH, side, color)
    end

    %% Plot freqs:
    axisH = subplotHandles.freq;
    Visuals.clear_old_lines(axisH)
    [L, R, f] = Algo.filteredFFTs(left, right, reqFreqs=freq_limits);    
    cellItems = {{L,left},{R,right}};
    for iC = 1 : length(cellItems)
        X   = cellItems{iC}{1};
        sig = cellItems{iC}{2};
        [~, color, side] = deal_data(iC, Sigs, cyclic_colors, Sides);
        hold( axisH, "on");
        plotH = plot( axisH, f, abs(X));
        line_props(plotH, side, color)
    end

    %% Plot BSI:
    axisH = subplotHandles.bsi;
    Visuals.clear_old_lines(axisH);
    bsi = Algo.singleBSI(L,R);        
    hold(axisH, "on");
    bar( subplotHandles.bsi, 0, bsi,1, "stacked", "blue");

    %% Gather IMFs and HHT information:
    [imf_tensor, ~] = Algo.bemd(Sigs);
    num_components = size(imf_tensor, 3);    
    IMFs_to_include = get_IMFs_to_include(num_components);
    % Pepare iteration data:
    nImfs = length(IMFs_to_include);
    nSides = 2;
    HHTs_cell = cell(nImfs, nSides);
    % iterate to collect data:
    ppb = Classes.PrintedProgressBar(nSides*nImfs, msg="Computing HHT for all IMFs:", print_length=70);
    for iSide = 1 : nSides
        for iImf = 1 : nImfs
            ppb.step()
            imf_num = IMFs_to_include(iImf);
            imf = imf_tensor(:, iSide, imf_num);
            [t, f, H] = Algo.hht(imf, fs, freqLimits=freq_limits);
            HHTs_cell{iImf, iSide} = struct(H=H, t=t, f=f);
        end
    end
    ppb.clear()   

    % IMFs lines names and colors:
    deal_imf_props = @(iImf) deal_imf_props_func(iImf, IMFs_to_include, cyclic_colors);
    %% Plot HHT:
    time_domain_start = min(time);
    for iSide = 1 : 2
        axisH = subplotHandles_hht_array(iSide);
        Visuals.clear_old_lines(axisH);
        for iImf = 1 : length(IMFs_to_include)  
            [~, color, name] = deal_imf_props(iImf);
            hht_struct = HHTs_cell{iImf, iSide};
            [t, f] = deal(hht_struct.t, hht_struct.f);
            time_vec = time_domain_start + t;
            hold(axisH, "on")
            scatter(axisH, time_vec, f, 3, color, "filled", DisplayName=name)
        end
    end
    %% Compute BAsy:
    [BAsyMat, BAsyFreqs] = Algo.basy(HHTs_cell, delta_freq, freq_to_include_in_BAsy);
    [BAsyTotal, BAsyPerImf] = Algo.total_basy(BAsyMat, BAsyFreqs, freqsToInclude=freq_to_include_in_BAsy, imfsToInclude=IMFs_to_include_in_BAsy);

    %% Plot BAsy plot:
    axisH = subplotHandles.BAsyPlot;
    Visuals.clear_old_lines(axisH);  
    for iImf = 1 : length(IMFs_to_include)        
        [~, color, name] = deal_imf_props(iImf);        
        hold(axisH, "on")
        BAsyVec = BAsyMat(:,iImf);
        plot(axisH, BAsyVec, BAsyFreqs, DisplayName=name, Color=color, LineWidth=2, Marker="*", MarkerSize=4);
    end

    %% Plot BAsy Bar:
    axisH = subplotHandles.BAsyBar;
    Visuals.clear_old_lines(axisH);        
    % Decide color by negative\positive:
    if BAsyTotal>=0
        color = "blue";   % Higher left side
    else
        color = "red";    % Higher right side 
    end
    BAsyTotal = abs(BAsyTotal);
    % Adjust limits:
    if BAsyTotal > axisH.YLim(2)
        axisH.YLim(2) = BAsyTotal;
    end
    % Plot:
    hold(axisH, "on");
    bar( axisH, 0, BAsyTotal, 1, "stacked", color);  

    %% Plot Accumulated BAsy Lines:
    axisH = subplotHandles.BAsyAccu;
    for iImf = 1 : length(IMFs_to_include)
        [imf_num, color, name] = deal_imf_props(iImf);       
        existing_plot = get_child_with_name(axisH, name);
        add_val = BAsyPerImf(imf_num);
        if isempty(existing_plot)
            hold(axisH, "on");
            plot(axisH, 1, add_val, LineWidth=3, Marker="*", DisplayName=name, Color=color);
        else
            prev_val = existing_plot.YData(end);            
            new_val = prev_val + add_val;
            new_x = existing_plot.XData(end) + 1;
            existing_plot.XData = [existing_plot.XData, new_x  ];
            existing_plot.YData = [existing_plot.YData, new_val];
        end
        
    end


end