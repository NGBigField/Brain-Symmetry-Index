classdef Side
    enumeration
        Right
        Left
    end

    methods (Static)
        function side = from_label(str)
            arguments
                str (1,1) string
            end
            import Classes.Side
   
            % Left:
            if     ismember(str, ["Fp1" "FP1" "F3" "F7" "C3" "T7" "P3" "P7" "O1" "T1" "T3" "T5" "A1"])   
                side = Side.Left;            
            % Right:
            elseif ismember(str, ["Fp2" "FP2" "F4" "F8" "C4" "T8" "P4" "P8" "O2" "T2" "T4" "T6" "A2"])
                side = Side.Right;                   
            else
                error("Unkown electrode placement");
            end

        end
    end
end