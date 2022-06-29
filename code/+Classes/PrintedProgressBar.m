classdef PrintedProgressBar < handle

    properties 
        expected_end   (1,1) uint32        
        print_length   (1,1) uint32
        msg            (1,1) string
        mark           (1,1) string
        num_marks      (1,1) uint32 = inf
    end

    properties (SetAccess=protected)
        static_printer (1,1) Classes.StaticPrinter 
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
            obj.expected_end   = expected_end;
            obj.print_length   = options.print_length;
            obj.msg            = options.msg;
            obj.mark           = options.mark;
            obj.static_printer = Classes.StaticPrinter();
        end
        %%
        function [] = step(obj)            
            % Count:
            obj.crnt_count = obj.crnt_count + 1;
            % Get num_out_of_num string:
            n_outof_n_str = StringUtils.num_out_of_num(obj.crnt_count, obj.expected_end);
            % Compute amount of marks and spaces to plot:
            fraction = double(obj.crnt_count) / double(obj.expected_end) ;
            remained_print_space = obj.print_length - strlength(obj.msg+" [] ()"+n_outof_n_str);            
            crnt_num_marks = round(remained_print_space*fraction);
            crnt_num_spaces = remained_print_space-crnt_num_marks;
            if crnt_num_marks == obj.num_marks
                return
            end
            % Get strings:
            [marks, spaces] = split_space(crnt_num_marks, crnt_num_spaces, obj.mark);
            text = obj.msg+" ["+marks+spaces+"] ("+n_outof_n_str+")";
            % Print:
            obj.static_printer.print(text)
            % Keep num_marks so we won't print the same thing twice.
            obj.num_marks = crnt_num_marks;            
        end
        %% 
        function [] = clear(obj)
            obj.static_printer.clear();
        end

    end

end
%% 

%% Subs:
function [marks, spaces] = split_space(nMarks, nSpaces, mark)
    marks  = StringUtils.repeat(mark, nMarks);
    spaces = StringUtils.repeat(" ", nSpaces);
end
