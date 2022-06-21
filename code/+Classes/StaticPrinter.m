classdef StaticPrinter < handle

    properties (Hidden)
        printed_length (1,1) uint64 = 0
        end_with_newline (1,1) logical = true
    end

    methods 

        %%
        function obj = StaticPrinter(options) 
            arguments
                options.end_with_newline (1,1) logical = true
            end
            obj.end_with_newline = options.end_with_newline;
        end
        %%
        function print(obj, str)
            arguments
                obj (1,1) Classes.StaticPrinter
                str (1,1) string
            end
            obj.clear()
            if obj.end_with_newline
                obj.printed_length = fprintf(str+newline);
            else
                obj.printed_length = fprintf(str);
            end
        end        
        %%        
        function add_print(obj, str)
            arguments
                obj (1,1) Classes.StaticPrinter
                str (1,1) string
            end
            obj.printed_length = obj.printed_length + fprintf(str);
        end
    end % methods

    methods 
        function [] = clear(obj)
            fprintf(repmat('\b',1,obj.printed_length))
        end
    end

end % classdef