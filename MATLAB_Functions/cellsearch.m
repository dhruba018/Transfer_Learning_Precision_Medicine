function index = cellsearch(cellarray, string, varargin)      % cellsearch(y_class, '1')
    
    opt = 0;            % print option w/ no match
    if nargin > 2
        opt = cell2mat(varargin);
    end
    
    index = [ ];
    if (size(cellarray, 1) == 1 || size(cellarray, 2) == 1)        
        for k = 1:length(cellarray)           
            if ~isnumeric(cellarray{k})
                if strcmpi(cellarray{k}, string)
                    index = [index; k];
                end
            else
                if cellarray{k} == string
                    index = [index; k];
                end
            end
        end
        
    else
        for i = 1:size(cellarray, 1)
            for j = 1:size(cellarray, 2)               
                if ~isnumeric(cellarray{i, j})
                  if strcmpi(cellarray{i, j}, string)
                        index = [index; [i, j]];
                  end
                else
                    if cellarray{i, j} == string
                        index = [index; [i, j]];
                    end
                end
            end
        end
    end
    
    if isempty(index) && ~opt
        disp('No Match!')
    end
    
end