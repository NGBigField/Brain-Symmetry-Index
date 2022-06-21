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
    end

end