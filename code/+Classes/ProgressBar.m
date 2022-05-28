classdef ProgressBar < handle
    properties
        endCount   (1,1) uint64
        crntCount  (1,1) uint64 = 0
        waitbar    (1,1) matlab.ui.Figure
    end

    methods 
        %%    
        function obj = ProgressBar(numExpectedCalls, msg)
            arguments
                numExpectedCalls    (1,1) uint64
                msg                 (1,1) string = "working..."
            end
            obj.endCount = numExpectedCalls;            
            obj.waitbar = init_figure(msg);
        end
        %%
        function [] = step(obj)
            obj.crntCount = 1 + obj.crntCount;
            obj.update_plot()
        end
        %% 
        function [] = close(obj)
            close(obj.waitbar)
        end
    end

    methods (Hidden)
        function [] = update_plot(obj)
            compFrac = obj.completion_fraction();
            waitbar( compFrac, obj.waitbar ) %#ok<*CPROP> 
        end
        function compFrac = completion_fraction(obj)
            compFrac = double(obj.crntCount) / double(obj.endCount);
            compFrac = min(1.00, compFrac);
        end
    end

end


function f = init_figure(msg)
    f = waitbar(0,msg);
    f.Color = [0.3, 0.5, 0.9];
    f.Name = "Progress Bar";
end