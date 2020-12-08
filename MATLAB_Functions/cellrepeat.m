function repeat = cellrepeat(cellarray, varargin)
    
    opt = 0;
    if nargin > 1
        opt = cell2mat(varargin);
    end
    
    repeat = [ ];
    for i = 1:length(cellarray)
        for j = i+1:length(cellarray)
            if strcmpi(cellarray(i), cellarray(j))
                repeat = [repeat; [i, j]];
            end
        end
    end
    
    if isempty(repeat) && ~opt
        disp('No Repeat!')
    end
end