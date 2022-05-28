classdef FigureRecorder
    properties %(Hidden)
        writer VideoWriter = VideoWriter.empty();
    end
    
    methods 
        function obj = FigureRecorder()
            % obj.writer = VideoWriter("tempVideo.avi", "Uncompressed AVI");
            obj.writer = VideoWriter("tempVideo.avi");
            obj.writer.FrameRate = 3;
            obj.writer.Quality = 100;
            obj.writer.open();
        end
        %%
        function [] = capture(obj, figH)
            arguments  
                obj  (1,1) Classes.FigureRecorder
                figH (1,1) = gcf
            end
            frame = getframe(figH);
%             [X, Map] = frame2im(frame);
            obj.writer.writeVideo(frame)
            
        end
        %%
        function [] = finish(obj)
            obj.writer.close()
        end

    end

end