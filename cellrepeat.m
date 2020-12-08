function repeat = cellrepeat(cellarray)
    repeat = [ ];
    for i = 1:length(cellarray)
        for j = i+1:length(cellarray)
            if strcmpi(cellarray(i), cellarray(j))
                repeat = [repeat; [i, j]];
            end
        end
    end
    
    if isempty(repeat)
        fprintf('>> No Repeat! \n')
    end
end