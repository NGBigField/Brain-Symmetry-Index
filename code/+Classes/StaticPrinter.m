classdef StaticPrinter < handle

    properties (Access=protected)
        printed_length   (1,1) {mustBeInteger} = 0
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
                str = str + newline;
            end
            obj.printed_length = fprintf(str);           
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
            if obj.printed_length==0
                return
            end
            str_that_clears_print = repmat('\b',1,obj.printed_length);
            fprintf(str_that_clears_print);
            obj.printed_length = 0;
        end
    end

end % classdef