classdef EEGSignal
    properties
        recording (:,1) double
        placement (1,1) Classes.ElectrodePlacement
        props     (1,1) struct
    end
end