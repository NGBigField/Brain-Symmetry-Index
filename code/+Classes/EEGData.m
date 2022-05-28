classdef EEGData
    properties
        left        (1,1) Classes.EEGSignal
        right       (1,1) Classes.EEGSignal
        time        (:,1) double
        condition   (1,1) Classes.Condition
        header      (1,1) struct
        props       (1,1) struct
    end
end