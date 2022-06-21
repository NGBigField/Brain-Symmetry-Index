%% clear
close all; clear all; clc;
addpath(genpath(pwd));

%% Constants:
Constants = struct;
Constants.epoch_duration = 50.0 ; % [sec]
Constants.epoch_overlap  = 10.0 ; % [sec]
Constants.Condition = Classes.Condition.Normal;

%% laod table
Table = readtable("data_mining_viable_data.xlsx");
disp(Table);

%% Parse Names:
eeg_file_names = unique( string( Table.file_name ) );

%% Init files to be saved:
Data = struct();
Data.date = datetime;
Data.EEGs = Classes.EEGData.empty();
Data.classDefFiles = load_files("+Classes");

%% Check file location:
fs = string(filesep);
folder = string(pwd)+fs+"Data"+fs+"tagged";
assert(isfolder(folder));
filename = StringUtils.time_str(Data.date);
fullpath = folder+fs+filename;

%% Get Data:
N = length(eeg_file_names);
prog_bar = Classes.ProgressBar(N);
for i = 1 : N
    prog_bar.step()

    eeg_file_name = eeg_file_names(i);      
        
    isSkip = check_skip(Constants, Table, eeg_file_name);
    if isSkip
        continue
    end

    % load file:
    [hdr, record] = edfread(proper_edf_file_name(eeg_file_name));    
    % derive requested data:
    [requested_data] = derive_requested_data(Constants, Table, hdr, eeg_file_name);

    for j = 1 : length(requested_data)
        label1 = requested_data(j).label1;
        label2 = requested_data(j).label2;

        [Data.EEGs] = get_data( Data.EEGs, Constants, hdr, record, [label1, label2], requested_data(j) );

    end % for j
end
prog_bar.close()

%% Save Data:
save(fullpath,"Data");
%% Finish
disp("End.");
%% End



%% Subs:
function [requested_data] = derive_requested_data(Constants, Table, hdr, target_name)
    arguments
        Constants   (1,1) struct
        Table             table
        hdr         (1,1) struct
        target_name (1,1) string
    end
    requested_data_cell = {};
    for i = 1 : size(Table,1)
        row = table2struct( Table(i,:) ); 
        row = StringUtils.to_string(row, 'string_types');
        name = row.file_name;

        % compare:
        if name == target_name        
            % add
            requested_data_cell{end+1,1} = row;
        end % if
 
    end % for i
    requested_data = struct_cell_2_struct_array(requested_data_cell);
end % func
%%
function [edf_file_name] = proper_edf_file_name(edf_file_name)
    [~,name,ext] = fileparts(edf_file_name);
    if ext == ""
        ext = "edf";
    end
    edf_file_name = name+"."+ext;
end


function [EEGs] = get_data( EEGs, Constants, hdr, record, labels, requested_data )
    arguments
        EEGs            (:,1) Classes.EEGData
        Constants       (1,1) struct
        hdr             (1,1) struct
        record          (:,:)
        labels          (:,1) string
        requested_data  (1,1) struct
    end
    
    % Derive basic input data:
    num_electrodes = size(record,1);
    % init indices tracking:
    signal_indices = [0,0]; 
    condition = Classes.Condition(requested_data.condition);

    %% Search requested electrodes indices:
    for i = 1 : num_electrodes
        % label:
        label = hdr.label{i};
        label = StringUtils.eeg_placement(label);    

        if ismember(label, labels)
            if     label == labels(1)
                signal_indices(1) = i;
            elseif label == labels(2)
                signal_indices(2) = i;
            end                        
        end % is label in [label1,label2]
    end % for i    
    assert(all(signal_indices~=0));

    %% Derive requested samples of data:
    % check matching frequencies 
    assert( hdr.frequency(signal_indices(1)) == hdr.frequency(signal_indices(2)) )
    Fs = hdr.frequency(signal_indices(1)); %[Hz == Samples/sec]
    % derive amount of samples per recording:
    samples_per_record = Constants.epoch_duration * Fs;
    samples_for_overlap = Constants.epoch_overlap * Fs;
    % prepare array of all sample times:
    requested_samples = double.empty(0,2);
    record_end = requested_data.start_samples;
    first = true;
    while record_end < requested_data.end_samples
        if first
            record_start = record_end;
        else
            record_start = record_end - samples_for_overlap ;
        end
        record_end = record_start + samples_per_record;
        requested_samples(end+1, [1,2]) = [record_start, record_end];
        first = false;
    end   


    %% Go over recording until end of requested samples
    
    for iS = 1 : size(requested_samples,1)
        % Init
        eeg_data = Classes.EEGData;
        [eeg_signal_left, eeg_signal_right, time] = get_signals( hdr, record, labels, signal_indices, requested_samples(iS,:) );

        eeg_data.left  = eeg_signal_left;
        eeg_data.right = eeg_signal_right;
        eeg_data.time  = time;
        eeg_data.condition = condition;

        eeg_data.header = hdr;
        eeg_data.props.indices_in_origin_recording = requested_samples(iS,:);

        % add to final array:
        EEGs(end+1,1) = eeg_data;

    end % while
    
end
%%
function [eeg_signal_left, eeg_signal_right, time] = get_signals( hdr, record, labels, signal_indices, samples_indices)

    % Derive inputs;
    N = length(labels);
    record_start = samples_indices(1);
    record_end   = samples_indices(2);
    indices = record_start : record_end;

    for i = 1 : N       
        label = labels(i);
        iSig = signal_indices(i);

        % Derive data:
        side = Classes.Side.from_label(label);
        Fs = hdr.frequency(iSig);
        time = indices / Fs;

        % Get signal:
        eeg_signal = Classes.EEGSignal();
        eeg_signal.recording = record(iSig,record_start:record_end);
        eeg_signal.placement = Classes.ElectrodePlacement();
        eeg_signal.placement.name = label;
        eeg_signal.placement.side = side;

        % Add extra info:
        eeg_signal.props = struct;
        eeg_signal.props.frequency = Fs;
        eeg_signal.props.prefilter = hdr.prefilter{iSig};

        % Assign to output:
        if side == Classes.Side.Left
            eeg_signal_left = eeg_signal;
        elseif side == Classes.Side.Right
            eeg_signal_right = eeg_signal;
        else
            error("Not a known side");
        end
    end % for i

end
%%
function out = load_files(folderName)
    % Constants:    
    s = @(c) string(c);
    fs = s(filesep);
    
    % init output;
    out = struct();

    % Go over all subfolders
    listing = dir(folderName);
    for i = 1 : length(listing)
        item = listing(i);
        name = item.name;
        if ismember( s(name), ["." ".."] )
            continue 
        end
        fullpath = s(item.folder)+fs+s(name);
        % Read file:
        text = string( textread(fullpath,"%s"));
        % save file:
        className = strsplit(s(name),'.');
        className = className(1);
        out.(className) = text;
    end

end
%%
function isSkip = check_skip(Constants, Table, eeg_file_name)
    if isempty(Constants.Condition)
        isSkip = false;
        return
    end
    fileRows = Table.file_name == eeg_file_name;
    conditions = unique(string(Table.condition(fileRows)));
    if ismember(string(Constants.Condition),conditions)
        isSkip = false;
    else
        isSkip = true;
    end


end