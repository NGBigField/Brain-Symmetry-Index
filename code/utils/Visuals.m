classdef Visuals
    %VISUALS Summary of this class goes here
    %   Detailed explanation goes here    
    
    methods (Static)
        function colors = distinguishable_colors(n_colors)
            colors = distinguishable_colors_full(n_colors);
        end
        %%
        function colors = default_colors()
            colors = colororder();
        end
        %%
        function color = get_default_color_cyclic(ind_in)
            colors = Visuals.default_colors();
            N = size(colors, 1);
            ind_used = mod(ind_in-1,N)+1;
            color = colors(ind_used, :);
        end
        %%
        function clear_old_lines(axisH)
            arguments
                axisH (1,1) matlab.graphics.axis.Axes
            end
            children = axisH.Children;
            for i = 1 : length(children)
                child = children(i);
                if ismember( string(child.Type), ["line", "bar", "scatter"]) ;
                    delete(children(i))
                end
            end % for i
        end % func
        %%
        function standart_bemd_fig_size(figH)
            arguments
                figH (1,1) matlab.ui.Figure
            end
            heigt = 760;
            figH.Position(2) = 1;
            figH.Position(4) = heigt;
            figH.Position(1) = 1;
            figH.Position(3) = heigt;
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
                    "lines", "colorcube", "prism", "flag", "white", "none" ])} = "none"
                options.value_limits (2,1) {mustBeNumeric} = [nan, nan]
                options.color  % Any
            end
            
            
            if options.colormap == "none"
                colors_out = cons_color(options.color, length(values) );
            else
                colors_out = colors_from_map(values, options.colormap, options.value_limits);
            end % if
                       
        end % func
        %%
        function save_fig(figH, name)
            arguments
                figH (1,1) matlab.ui.Figure
                name (1,1) string
            end
            folder = get_desktop_path();
            fullname = folder+filesep+name;
            saveas(figH, fullname, "svg")
        end
        %%
        function save_all_figs()
            figHandles = findobj('Type', 'figure');
            sp = Classes.StaticPrinter();
            N = length(figHandles);
            for i = 1 : N
                prog_str = StringUtils.num_out_of_num(i, N, separator=" out of ");
                sp.print("Saving figure " + prog_str );
                figH = figHandles(i);
                name = "temp"+string(i);
                Visuals.save_fig(figH, name)
            end
        end
        %%
        function subaxes = get_subaxes(figH)
            subaxes = [];
            children = figH.Children;
            for iC = 1 : length(children)
                child = children(iC);
                if string(child.Type) == "axes"
                    subaxes = [child, subaxes];
                end
            end % for iC
        end
        %%
        function change_axis_width(axis, axes_to_the_right, factor)
            % calc width:
            prevWidth = axis.Position(3);
            newWidth = prevWidth * factor;
            % Assign:
            axis.Position(3) = newWidth;
            % adjust next axes:
            dW = newWidth-prevWidth;
            for i = 1 : length(axes_to_the_right)
                drawnow;
                axis_to_the_right = axes_to_the_right(i);
                axis_to_the_right.Position(1) = axis_to_the_right.Position(1) + dW;
            end
        end
        %%
        function sequencialy_change_axes_widths(sorted_axes_array, factor)
            for i = 1 : length(sorted_axes_array)
                Visuals.change_axis_width(sorted_axes_array(i), sorted_axes_array(i+1:end), factor);
            end
        end
        %% 
        function index = grid_xy_to_index(m,n,y,x)
            arguments
                m (1,1) {mustBeInteger} % grid dimensions 1
                n (1,1) {mustBeInteger} % grid dimensions 2
                y (1,:) {mustBeInteger} % axis location 1 
                x (1,:) {mustBeInteger} % axis location 2 
            end
            index = n*(y-1)+x;
        end
        %%

    end % methods
end % class

%% End class

%% Subs:

function colors_out = colors_from_map(values, colormap, value_limits)
    % Derive inputs:
    c_map = colormap(colormap);
    num_colors = size(c_map, 1);

    if all(isnan(value_limits))
        value_limits(1) = min(values, [], 'all');
        value_limits(2) = max(values, [], 'all');
    end

    % Derive color indices:
    mapped_vals = linspace( value_limits(1), value_limits(2), num_colors);

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
    end % for
end
%%
function colors_out = cons_color(color, len)
    colors_out = repmat(color, len, 1);    
end
%%
function colors = distinguishable_colors_full(n_colors,bg,func)
% DISTINGUISHABLE_COLORS: pick colors that are maximally perceptually distinct
%
% Credit to: Tim Holy (2022). 
% Generate maximally perceptually-distinct colors (https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors), 
% MATLAB Central File Exchange. 
% Retrieved May 16, 2022
%
% When plotting a set of lines, you may want to distinguish them by color.
% By default, Matlab chooses a small set of colors and cycles among them,
% and so if you have more than a few lines there will be confusion about
% which line is which. To fix this problem, one would want to be able to
% pick a much larger set of distinct colors, where the number of colors
% equals or exceeds the number of lines you want to plot. Because our
% ability to distinguish among colors has limits, one should choose these
% colors to be "maximally perceptually distinguishable."
%
% This function generates a set of colors which are distinguishable
% by reference to the "Lab" color space, which more closely matches
% human color perception than RGB. Given an initial large list of possible
% colors, it iteratively chooses the entry in the list that is farthest (in
% Lab space) from all previously-chosen entries. While this "greedy"
% algorithm does not yield a global maximum, it is simple and efficient.
% Moreover, the sequence of colors is consistent no matter how many you
% request, which facilitates the users' ability to learn the color order
% and avoids major changes in the appearance of plots when adding or
% removing lines.
%
% Syntax:
%   colors = distinguishable_colors(n_colors)
% Specify the number of colors you want as a scalar, n_colors. This will
% generate an n_colors-by-3 matrix, each row representing an RGB
% color triple. If you don't precisely know how many you will need in
% advance, there is no harm (other than execution time) in specifying
% slightly more than you think you will need.
%
%   colors = distinguishable_colors(n_colors,bg)
% This syntax allows you to specify the background color, to make sure that
% your colors are also distinguishable from the background. Default value
% is white. bg may be specified as an RGB triple or as one of the standard
% "ColorSpec" strings. You can even specify multiple colors:
%     bg = {'w','k'}
% or
%     bg = [1 1 1; 0 0 0]
% will only produce colors that are distinguishable from both white and
% black.
%
%   colors = distinguishable_colors(n_colors,bg,rgb2labfunc)
% By default, distinguishable_colors uses the image processing toolbox's
% color conversion functions makecform and applycform. Alternatively, you
% can supply your own color conversion function.
%
% Example:
%   c = distinguishable_colors(25);
%   figure
%   image(reshape(c,[1 size(c)]))
%
% Example using the file exchange's 'colorspace':
%   func = @(x) colorspace('RGB->Lab',x);
%   c = distinguishable_colors(25,'w',func);
% Copyright 2010-2011 by Timothy E. Holy
  % Parse the inputs
  if (nargin < 2)
    bg = [1 1 1];  % default white background
  else
    if iscell(bg)
      % User specified a list of colors as a cell aray
      bgc = bg;
      for i = 1:length(bgc)
	bgc{i} = parsecolor(bgc{i});
      end
      bg = cat(1,bgc{:});
    else
      % User specified a numeric array of colors (n-by-3)
      bg = parsecolor(bg);
    end
  end
  
  % Generate a sizable number of RGB triples. This represents our space of
  % possible choices. By starting in RGB space, we ensure that all of the
  % colors can be generated by the monitor.
  n_grid = 30;  % number of grid divisions along each axis in RGB space
  x = linspace(0,1,n_grid);
  [R,G,B] = ndgrid(x,x,x);
  rgb = [R(:) G(:) B(:)];
  if (n_colors > size(rgb,1)/3)
    error('You can''t readily distinguish that many colors');
  end
  
  % Convert to Lab color space, which more closely represents human
  % perception
  if (nargin > 2)
    lab = func(rgb);
    bglab = func(bg);
  else
    C = makecform('srgb2lab');
    lab = applycform(rgb,C);
    bglab = applycform(bg,C);
  end
  % If the user specified multiple background colors, compute distances
  % from the candidate colors to the background colors
  mindist2 = inf(size(rgb,1),1);
  for i = 1:size(bglab,1)-1
    dX = bsxfun(@minus,lab,bglab(i,:)); % displacement all colors from bg
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
  end
  
  % Iteratively pick the color that maximizes the distance to the nearest
  % already-picked color
  colors = zeros(n_colors,3);
  lastlab = bglab(end,:);   % initialize by making the "previous" color equal to background
  for i = 1:n_colors
    dX = bsxfun(@minus,lab,lastlab); % displacement of last from all colors on list
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
    [~,index] = max(mindist2);  % find the entry farthest from all previously-chosen colors
    colors(i,:) = rgb(index,:);  % save for output
    lastlab = lab(index,:);  % prepare for next iteration
  end
end
function c = parsecolor(s)
  if ischar(s)
    c = colorstr2rgb(s);
  elseif isnumeric(s) && size(s,2) == 3
    c = s;
  else
    error('MATLAB:InvalidColorSpec','Color specification cannot be parsed.');
  end
end
function c = colorstr2rgb(c)
  % Convert a color string to an RGB value.
  % This is cribbed from Matlab's whitebg function.
  % Why don't they make this a stand-alone function?
  rgbspec = [1 0 0;0 1 0;0 0 1;1 1 1;0 1 1;1 0 1;1 1 0;0 0 0];
  cspec = 'rgbwcmyk';
  k = find(cspec==c(1));
  if isempty(k)
    error('MATLAB:InvalidColorString','Unknown color string.');
  end
  if k~=3 || length(c)==1,
    c = rgbspec(k,:);
  elseif length(c)>2,
    if strcmpi(c(1:3),'bla')
      c = [0 0 0];
    elseif strcmpi(c(1:3),'blu')
      c = [0 0 1];
    else
      error('MATLAB:UnknownColorString', 'Unknown color string.');
    end
  end
end
