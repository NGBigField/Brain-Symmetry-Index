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
    end

end