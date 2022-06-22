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
            obj = Classes.StaticPrinter();
        end
        %%
        function str = formatted(value, options)
            arguments
                value  (1,1)                
                options.signed   (1,1) logical = false
                options.pad      (1,1) string {mustBeMember(options.pad, ["spaces", "zeros"])} = "spaces"
                options.width    (1,1) int16 = nan
                options.decimal  (1,1) int16 = nan      
                options.notation (1,1) string {mustBeMember(options.notation, ["c","d","e","E","f","g","G","o","s","u","x","X","default"])} = "default"                         
            end
            % Helper functions:
            s = @(x) string(x);

            % Build format string:
            format = "%";
            if options.signed
                format = format + "+";
            end
            if options.pad == "zeros"
                format = format + "0";
            end
            if ~isnan(options.width)
                format = format + s(options.width) ;
            end
            if ~isnan(options.decimal)
                format = format + "." + s(options.decimal);
            end
            if options.notation == "default"
                if ~isnan(options.decimal)
                    options.notation = "f";
                else
                    options.notation = "d";
                end
            end
            format = format + options.notation ;
            
            % get string:
            str = sprintf(format, value);

        end
        %%
        function str = num_out_of_num(num1, num2, options)
            arguments
                num1 (1,1) int32
                num2 (1,1) int32
                options.separator (1,1) string = "/"
            end
            % Define helper funcs:
            width = strlength(string(num2));
            formatted = @(x) StringUtils.formatted(x, width=width);
            % get string:
            str = formatted(num1) + options.separator + formatted(num2);
        end
        %%
    end
end

%% subs:
function tf = is_string_type(input)
    tf = ischar(input) || isstring(input);
end