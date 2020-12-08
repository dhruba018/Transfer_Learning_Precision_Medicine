function index = cellnan(cellarray)
    index = [ ];
    for k = 1:length(cellarray)
        if isnan(cellarray{k})
            index = [index; k];
        end
    end
end