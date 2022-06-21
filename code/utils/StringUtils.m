classdef StringUtils

    methods (Static)
        %%
        function out = eeg_placement(in)              
            % Constants
            expectedPrefix = "EEG";
            expectedSuffix = "REF";

            % Derive from imput:
            L = @(s) strlength(s);
            prefix = string(in( 1 : L(expectedPrefix) ));
            suffix = string(in( (end-L(expectedSuffix)+1) : end ));

            % Assign output:
            out = in;            

            % Check and remove:
            if prefix == expectedPrefix % Remove prefix
                out = out((L(expectedPrefix)+1):end);
            end
            if suffix == expectedSuffix % Remove suffix
                out = out(1:(end-L(expectedSuffix)));
            end

            % String:
            out = string(out);
        end
        %%
        function in_struct = to_string(in_struct, what)
            arguments
                in_struct   (1,1) struct
                what        (1,1) string {mustBeTextScalar, mustBeMember(what,["everything", "string_types"])} = "everything"
            end
            field_names = string(fieldnames(in_struct));
            for i = 1 : length(field_names)
                fieldname = field_names(i);
                switch what 
                    case "everything"
                        % pass
                    case "string_types"
                        if ~is_string_type(in_struct.(fieldname))
                            continue
                        end
                    otherwise
                        error("Not a supported argument")
                end                
                in_struct.(fieldname) = string(in_struct.(fieldname));
            end
        end
        %%
        function str = time_str(d)
            arguments
                d (1,1) datetime = datetime
            end
            chars = char(d);
            for iC = 1 : length(chars)
                c = chars(iC);
                if c == ':'
                    c = '-';
                end
                chars(iC) = c;
            end
            str = string(chars);
        end
        %%
        function obj = static_printer()
            obj = Classes.StaticPrinter()
        end
    end
end

%% subs:
function tf = is_string_type(input)
    tf = ischar(input) || isstring(input);
end