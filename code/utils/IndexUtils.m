classdef IndexUtils

    methods (Static)
        %%
        function out = index_mat_or_cell(in, index, val)
            arguments
                in {mustBeNonempty}
                index {mustBeInteger,mustBeGreaterThan(index,0)}
                val = [] % optional
            end

            if isempty(val)
                %% setting a value to input:

            else 
                %% getting a value from input:

            end
        end
        %%
        function index = find_closer(val, ref)
            arguments
                val (1,1) {mustBeNumeric}
                ref (:,1) {mustBeNumeric}
            end
            distance = abs(val-ref);
            [~, index] = min(distance);
        end
        %%
        function tf = is_member(val, list)
            tf = false;
            if isempty(list)
                return                
            end
            if iscell(list)
                for i = 1 : length(list)
                    elem = list{i};
                    if elem == val
                        tf = true;
                        return
                    end
                end % for i
            end % if cell
            tf = ismember(val, list);
        end
    end

end % class def

%%

%% Subs:
function tf = is_scalar(x)
    [r,c] = size(x);
    if (r > 1) || (c > 1)
        tf = 0;
    else
        tf = 1;
    end
end