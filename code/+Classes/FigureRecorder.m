classdef FigureRecorder
    properties %(Hidden)
        writer VideoWriter = VideoWriter.empty();
    end
    
    methods 
        function obj = FigureRecorder()
            % Parse full path:
            name = "tempVideo.avi";
            folder = get_desktop_path();
            fullpath = folder+filesep+name;
            % Save objects data:            
            obj.writer = VideoWriter(fullpath);
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
            obj.writer.writeVideo(frame)            
        end
        %%
        function [] = finish(obj)
            obj.writer.close()
        end
        %%
        function delete(obj)
            obj.finish();
        end

    end

end