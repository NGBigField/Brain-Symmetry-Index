classdef ProgressBar < handle

    properties (SetAccess=protected)
        endCount        (1,1) uint64
        crntCount       (1,1) uint64 = 0
        waitBar         (1,1) matlab.ui.Figure
        cancelRequested (1,1) logical = false        
    end
    
    properties (Constant)
        closingFigureCausesCancelRequest (1,1) logical = false
    end

    methods 
        %%    
        function obj = ProgressBar(numExpectedCalls, options)
            arguments
                numExpectedCalls    (1,1) uint64
                options.msg         (1,1) string  = "Progress Bar"
                options.cancelable  (1,1) logical = false;
            end
            obj.endCount = numExpectedCalls;            
            obj.waitBar = obj.init_figure(options.msg, options.cancelable);
        end
        %%
        function cancel_requested = step(obj)
            cancel_requested = obj.cancelRequested;
            if cancel_requested
                return
            end               
            obj.crntCount = 1 + obj.crntCount;
            try
                obj.update_plot()
            catch ME
                error_msg = string(ME.getReport);
                warning(error_msg);
            end
        end
        %%
        function [] = request_cancle(obj, src, event)                              
            if isa(src, "matlab.ui.control.UIControl") %% User pressed on cancel:
                 obj.cancelRequested = true;
            elseif isa(src, "matlab.ui.Figure") %% User closed the figure:
                if obj.closingFigureCausesCancelRequest
                    obj.cancelRequested = true;
                end
                obj.close_figure();
            else
                error("No such option");
            end
        end
        %% 
        function close(obj)
            obj.close_figure()
        end
        %%
        function delete(obj)
            obj.close_figure();
        end

    end

    methods (Access=protected)
        function close_figure(obj)
            figH = obj.waitBar;
            if isgraphics(figH)
                close(figH, "force");
            end
        end
        %%
        function [] = update_plot(obj)
            compFrac = obj.completion_fraction();
            progText = obj.prog_text();
            waitbar( compFrac, obj.waitBar, progText  ) 
        end
        %%
        function compFrac = completion_fraction(obj)
            compFrac = double(obj.crntCount) / double(obj.endCount);
            compFrac = min(1.00, compFrac);
        end
        %% 
        function  str = prog_text(obj)
            str = StringUtils.num_out_of_num(obj.crntCount, obj.endCount);
        end
        %%
        function f = init_figure(obj, msg, cancelable)
            progText = obj.prog_text();
            if cancelable
                f = waitbar(0, progText, "CreateCancelBtn", @(src,event) obj.request_cancle(src,event) );
                button = f.Children(1);
                button.BackgroundColor = [0.9, 0.2, 0.3];
                button.ForegroundColor = [0.9, 0.9, 0.9];
            else
                f = waitbar(0, progText );
            end
            f.Name = msg;
            % Colors:
            f.Color = [0.3, 0.5, 0.9];
            
        end
    end % methods
end % class

