classdef BrainWavesFreqs 
    enumeration
        Alpha %  8   to  12 - Relaxed but not drowsy
        Beta  % 13   to  30 - Thinking, Aware of self and surrounding
        Gamma % 30   to 100 - Higher Mental Activity, motor function
        Delta %  0.5 to   3 - Deep, dream less sleep
        Theta %  4   to   7 - Fantasy, imaginary dream
    end

    methods
        function res = freqs(obj)
            arguments
                obj Classes.BrainWavesFreqs
            end
            import Classes.BrainWavesFreqs
            switch obj
                case BrainWavesFreqs.Alpha 
                    res = [  8  ,  12 ];  %  Relaxed but not drowsy
                case BrainWavesFreqs.Beta  
                    res = [ 13  ,  30 ];  %  Thinking, Aware of self and surrounding
                case BrainWavesFreqs.Gamma 
                    res = [ 30  , 100 ];  %  Higher Mental Activity, motor function
                case BrainWavesFreqs.Delta 
                    res = [  0.5,   3 ];  %  Deep, dream less sleep
                case BrainWavesFreqs.Theta 
                    res = [  4  ,   7 ];  %  Fantasy, imaginary dream
                otherwise
                    error("Not a possible brain-wave-frequency");
            end
        end 
        %%
        function res = association(obj)
            arguments
                obj Classes.BrainWavesFreqs
            end
            import Classes.BrainWavesFreqs
            switch obj
                case BrainWavesFreqs.Alpha 
                    res = "Relaxed but not drowsy";
                case BrainWavesFreqs.Beta  
                    res = "Thinking, Aware of self and surrounding";
                case BrainWavesFreqs.Gamma 
                    res = "Higher Mental Activity, motor function";
                case BrainWavesFreqs.Delta 
                    res = "Deep, dream less sleep";
                case BrainWavesFreqs.Theta 
                    res = "Fantasy, imaginary dream";
                otherwise
                    error("Not a possible brain-wave-frequency");
            end
        end
        %%

    end
end