function [a] = struct_cell_2_struct_array(c)
    arguments
        c (:,1) cell
    end
    %% Derive Fileds:
    first = c{1};
    field_names = string(fieldnames(first));
    nC = length(c);

    %% Init Array:
    a = struct.empty(nC,0);

    %% Copy data:
    for iC = 1 : length(c)
        for iF = 1 : length(field_names)
            field = field_names(iF);            
            a(iC,1).(field) = c{iC}.(field);
        end
    end
end