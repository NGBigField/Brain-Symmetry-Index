classdef PrintedProgressBar < handle

    properties 
        expected_end   (1,1) uint32        
        print_length   (1,1) uint32
        msg            (1,1) string
        mark           (1,1) string
        num_marks      (1,1) uint32 = inf
    end

    properties (SetAccess=protected)
        static_printer (1,1) Classes.StaticPrinter = Classes.StaticPrinter()
        crnt_count     (1,1) uint32 = 0
    end


    methods 
        %%
        function obj = PrintedProgressBar(expected_end, options )
            arguments
                expected_end         (1,1) uint32
                options.print_length (1,1) uint32 = 100
                options.msg          (1,1) string = "Working:"
                options.mark         (1,1) string = "#"
            end
            obj.expected_end = expected_end;
            obj.print_length = options.print_length;
            obj.msg          = options.msg;
            obj.mark         = options.mark;
        end
        %%
        function [] = step(obj)            
            obj.crnt_count = obj.crnt_count + 1;
            fraction = double(obj.crnt_count) / double(obj.expected_end) ;
            remained_print_space = obj.print_length - strlength(obj.msg+" []");            
            crnt_num_marks = round(remained_print_space*fraction);
            if crnt_num_marks == obj.num_marks
                return
            end
            [marks, spaces] = split_space(remained_print_space, crnt_num_marks, obj.mark);
            text = obj.msg+" ["+marks+spaces+"]";
            obj.static_printer.print(text)
            obj.num_marks = crnt_num_marks;            
        end

    end

end
%% 

%% Subs:
function [marks, spaces] = split_space(remained_print_space, num_marks, mark)
    if num_marks == 0
        marks  = "";
        spaces = strjoin( repmat(" ", 1, remained_print_space) );        
    elseif num_marks == remained_print_space
        marks  = strjoin( repmat(mark, 1, num_marks) );
        spaces = "";
    else
        marks  = strjoin( repmat(mark, 1, num_marks) );
        spaces = strjoin( repmat(" ", 1, remained_print_space-num_marks) );
    end
end