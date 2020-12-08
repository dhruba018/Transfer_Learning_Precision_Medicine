function [C, varargout] = intersectn(varargin)
% INTERSECTN performs intersection among N input sets.
    % NOTATION: function [C, I1, I2,...] = INTERSECTN(A1, A2,... , setOrder)
    %   Inputs:
    %       A1, A2,... => Input sets
    %       setOrder (optional) => Order of the intersected set. Use either 'sorted' or 'stable'. 
    %                                           Default output is 'sorted'.
    %   Outputs:
    %       C => Intersected set. Order defined by setOrder input.
    %       I1, I2,... => Indices of the elements in C in A1, A2,... Order defined by setOrder input.
    %
    % (c) 2018 S.R.Dhruba
    
    % Check inputs...
    sortOrders = {'sorted', 'stable'};
    if nargin < 2
        error('Insufficient inputs!!! At least two input sets are required.')
    else
        orderFlag = 0;
        if and(ischar(class(varargin{end})), length(varargin{end}) == 6)
            orderFlag = sum(strcmpi(sortOrders, varargin{end}));
        end
        numSets = nargin - orderFlag;                      % #input sets
        inClasses = cellfun(@class, varargin(1 : numSets), 'uniformoutput', 0);
        switch (sum(cell2mat(regexpi(inClasses, inClasses(1)))) == numSets)     % Check if input classes are same
            case 1
                if orderFlag > 0
                    setOrder = varargin{end};
                else
                    setOrder = sortOrders{1};
                end
                A = varargin(1 : numSets);
            case 0
                error('All input sets must have the same class!!!')
        end
    end
    
    % Perform intersection & retrieve indices...
    C = intersect(A{1}, A{2});
    for k = 3 : numSets
        C = intersect(C, A{k});                                  % Get the overall intersection
    end
    C = intersect(A{1}, C, setOrder);                      % Get the ordered intersection set
    
    Aidx = cell(1, numSets);
    for k = 1 : numSets
        [~, ~, Aidx{k}] = intersect(C, A{k}, setOrder);  % Get the contributing element indices for each set
    end
    varargout = Aidx;
    
end

