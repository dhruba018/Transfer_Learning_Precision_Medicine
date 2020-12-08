function index = cellempty(cellarray, varargin)
    
    opt = 0;
    if nargin > 1
        opt = cell2mat(varargin);
    end
    
    index = [ ];
    for k = 1:length(cellarray)
        if isempty(cellarray{k})
            index = [index; k];
        end
    end
    
    if isempty(index) && ~opt
        disp('No empty cell!')
    end
end