function index = cellnumeric(cellarray)
    index = [ ];
    for k = 1:length(cellarray)
        if isnumeric(cellarray{k})
            index = [index; k];
        end
    end
end