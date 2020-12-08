function index = cellempty(cellarray)
    index = [ ];
    for k = 1:length(cellarray)
        if isempty(cellarray{k})
            index = [index; k];
        end
    end
end