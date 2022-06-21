classdef Visuals
    %VISUALS Summary of this class goes here
    %   Detailed explanation goes here    
    
    methods (Static)
        function standart_bemd_fig_size(figH)
            arguments
                figH (1,1) matlab.ui.Figure
            end
            screen_size = Visuals.get_screen_size();
            screen_height = screen_size(2);
            figH.Position(2) = 1;
            figH.Position(4) = screen_height;
            figH.Position(1) = 1;
            figH.Position(3) = screen_height;
        end
        %%
        function result = get_screen_size()
            prev_units = get(0,'units');
            set(0,'units','pixels')
            result = get(0,'screensize');
            set(0,'units', prev_units)
            result = result(3:4);
        end
        %%
        function [] = link_axes(figH, dim)
            arguments                
                figH (1,1) matlab.ui.Figure
                dim (1,1) string {mustBeMember(dim, ["xyz" , "x", "y", "z", "xy", "xz", "yz", "off" ] )} = "xyz"
            end
            children = figH.Children;
            axes_list = [];
            for iC = 1 : length(children)
                child = children(iC);
                if string(child.Type) == "axes"
                    axes_list = [axes_list, child]; %#ok<AGROW> 
                end
            end
            linkaxes(axes_list, dim)
        end
        %%
        function colors_out = color_vector(values, options)
            arguments
                values (:,1) {mustBeNumeric}
                options.colormap (1,1) string {mustBeMember(options.colormap, [ ...
                    "parula", "turbo", "hsv", "hot", "cool", "spring", "summer", ...
                    "autumn", "winter", "gray", "bone", "copper", "pink", "jet", ...
                    "lines", "colorcube", "prism", "flag", "white"])} = "parula"
                options.value_limits (2,1) {mustBeNumeric} = [nan, nan]
            end
            
            % Derive inputs:
            c_map = colormap(options.colormap);
            num_colors = size(c_map, 1);

            if all(isnan(options.value_limits))
                options.value_limits(1) = min(values, [], 'all');
                options.value_limits(2) = max(values, [], 'all');
            end

            % Derive color indices:
            mapped_vals = linspace( options.value_limits(1), options.value_limits(2), num_colors);

            % Init output:
            Len = length(values);
            colors_out = zeros(Len, 3);
            
            % Assign colors to values:
            for i = 1 : Len
                val = values(i);
                if val <= options.value_limits(1)
                    color = c_map(1,:);
                elseif val >= options.value_limits(2)
                    color = c_map(end,:);
                else
                    index = IndexUtils.find_closer(val, mapped_vals);
                    color = c_map(index,:);
                end
                colors_out(i,:) = color;
            end            
            
            
        end % func
    end
end

