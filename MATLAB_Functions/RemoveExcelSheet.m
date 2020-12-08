function RemoveExcelSheet(ExcelFileName, SheetName)
    % RemoveExcelSheet - removes the sheets that are automatically added to excel file. 
    % When MATLAB writes data to a new Excel file, the Excel software automatically 
    % creates 3 sheets (the names are depended on the user language). This appears 
    % even if the user has defined the sheet name to be added. 
    %
    % Usage:
    % RemoveExcelSheet(ExcelFileName) 
    %       removes "Sheet1", "Sheet2", "Sheet3" from the  excel file. ExcelFileName 
    %       is a string of the Excel file name.
    % RemoveExcelSheet(ExcelFileName, SheetName) 
    %       enables the user to enter the sheet name when the language is other than 
    %       English. SheetName is the default sheet name, WITHOUT the number.
    %
    %                       Written by Noam Greenboim
    %                       www.perigee.co.il
    %
    
    %% Check Input Arguments.
    if nargin < 1 || isempty(ExcelFileName)
        error('Filename must be specified.')
    end
    
    if ~ischar(ExcelFileName)
        error('Filename must be a string.')
    end
    
    try
        ExcelFileName = ValidPath(ExcelFileName);
    catch 
        error('File not found.')
    end
    
    if nargin < 2
        SheetName = 'Sheet';            % EN: Sheet, DE: Tabelle, HE: גיליון , etc. (Language dependent)
    else
        if ~ischar(SheetName)
            error('Default sheet name must be a string.')
        end
    end

    %%
    % Open Excel File.
    objExcel = actxserver('Excel.Application');
    objExcel.Workbooks.Open(ExcelFileName);         % Full path is necessary!
    
    % Delete sheets.
    try
        % Throws an error if the sheets do not exist.
        objExcel.ActiveWorkbook.Worksheets.Item([SheetName, '1']).Delete;
        fprintf('\nsheet #1 - deleted.')
        objExcel.ActiveWorkbook.Worksheets.Item([SheetName, '2']).Delete;
        fprintf('\nsheet #2 - deleted.')
        objExcel.ActiveWorkbook.Worksheets.Item([SheetName, '3']).Delete;
        fprintf('\nsheet #3 - deleted.\n')
    catch
        fprintf('\n')
        O = objExcel.ActiveWorkbook.Worksheets.get;
        if O.Count == 1
            error('Can''t delete the last sheet! Excel file must containt at least one sheet.')
        else
            % warning('Problem occured. Check excel file.')
            warning('Sheets don''t exist! Check excel file.')
        end
    end
    
    % Save, close and clean up.
    objExcel.ActiveWorkbook.Save;
    objExcel.ActiveWorkbook.Close;
    objExcel.Quit;
    objExcel.delete;
end

function FileNameOut = ValidPath(FileName)
    % ValidPath builds a full path from a partial path specification
    % FileName = ValidPath(FileName) returns a string vector containing full path to a file. 
    % FileName is string vector containing a partial path ending in a file or directory name. 
    % May contain ..\  or ../ or \\. The current directory (pwd) is prepended to create a full 
    % path if necessary. On UNIX, when the path starts with a tilde, '~', then the current 
    % directory is not prepended.
    %
    % See also XLSREAD, XLSWRITE, XLSFINFO.
    
    % Copyright 1984-2012 The MathWorks, Inc.
    
    % First check for wild cards, since that is not supported.
    if strfind(FileName, '*') > 0
        error(message('MATLAB:xlsread:Wildcard', FileName));
    end
    
    % Break partial path in to file path parts.
    [Directory, File, Ext] = fileparts(FileName);
    if ~isempty(Ext)
        FileNameOut = getFullName(FileName);
    else
        ExtIn = matlab.io.internal.xlsreadSupportedExtensions;
        for i = 1:length(ExtIn)
            try                                                           %#ok<TRYNC>
                FileNameOut = getFullName(fullfile(Directory, [File, ExtIn{i}]));
                return
            end
        end
        error(message('MATLAB:xlsread:FileDoesNotExist', FileName))
    end
end

function AbsolutePath = AbsPath(PartialPath)
    % Parse partial path into path parts
    [PathName, FileName, Ext] = fileparts(PartialPath);
    
    % No path qualification is present in partial path; assume parent is pwd, except
    % when path string starts with '~' or is identical to '~'.
    if isempty(PathName) && strncmp('~', PartialPath, 1)
        Directory = pwd;
    elseif (isempty(regexp(PartialPath, '(.:|\\\\)', 'once')) && ...
            ~strncmp('/', PartialPath, 1) && ...
            ~strncmp('~', PartialPath, 1))
        % Path did not start with any of drive name, UNC path or '~'.
        Directory = [pwd, filesep, PathName];
    else
        % Path content present in partial path; assume relative to current directory or absolute.
        Directory = PathName;
    end
    
    % Construct absolute filename.
    AbsolutePath = fullfile(Directory, [FileName, Ext]);
end

function FileName = getFullName(FileName)
    FileOnPath = which(FileName);
    if isempty(FileOnPath)
        % Construct full path to source file
        FileName = AbsPath(FileName);
        if isempty(dir(FileName)) && ~isdir(FileName)
            % File does not exist. Terminate importation of file.
            error(message('MATLAB:xlsread:FileDoesNotExist', FileName))
        end
    else
        FileName = FileOnPath;
    end
end