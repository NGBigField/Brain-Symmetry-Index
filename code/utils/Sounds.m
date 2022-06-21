classdef Sounds
    methods (Static)
        function gong(vol,frq,dur)
            % gong: sounds gong
            % John Gooden (2022). Gong (https://www.mathworks.com/matlabcentral/fileexchange/21231-gong), MATLAB Central File Exchange. Retrieved June 19, 2022.
            % 
            % call gong
            % call gong(vol)
            % call gong(vol,frq)
            % call gong(vol,frq,dur)
            %
            % input arguments (optional, if 0 then default taken)
            % vol = volume (default = 1)
            % frq = base frequency (default = 440 Hz)
            % dur = duration (default = 1 s)
            fb  = 440;
            td  = 1;
            vl  = 1;
            if nargin>=1
                if vol>0 vl = vol; end
            end
            if nargin>=2
                if frq>0 fb = frq; end
            end
            if nargin>=3
                if dur>0 td = dur; end
            end
                
            t   =[0:8192*td]'/8192;
            env = exp(-5*t/td);
            f   = fb;
            vol = 0.3*vl;
            tpft = 2*pi*f*t;
            sl  = sin(tpft)+0.1*sin(2*tpft)+0.3*sin(3*tpft);
            sl  = vol*sl;
            sr  = [sl(100:end);sl(1:99)];
            vl  = cos(20*t);
            vr  = 1-vl;
            y(:,1) = vl.*env.*sl;
            y(:,2) = vr.*env.*sr;
            sound(y)
        end 
    end
end