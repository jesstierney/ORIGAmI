function showplot(Files,bounds)
    
    % showplot plots the contents of each file. The files can either be
    % .csv or .txt files. You can include only a portion of the data by
    % defining time boundaries in the variable 'bounds'.
    
    % Inputs:
    % Files - names of files to plot. If there is only one
    %       file, the input can be a string. If there are multiple files,
    %       the input should be a cell array of strings.
    %           e.g. {'GC27_1.csv','GC27_3.csv','GC27_4.csv',...
    %                 'KNR197_1.csv','KNR197_2.csv','KNR197_3.csv',...
    %                 'P178_1.csv','P178_2.csv','P178_3.csv'};
    % bounds - [T1, T2], where T1 is the minimum time to include in the
    % plot, and T2 is the maximum time.
    %           Default: Plot all values.
    
    
    % Example: showplot('1022.txt');
    % Example: showplot({'1022.txt';'1036.txt'});
    % Example: showplot({'1022.txt';'1036.txt'},[10 30]);
    
    
    if ischar(Files)
        Files={Files};
    end
    
    for f=1:length(Files)
        
        File=Files{f};
        
        % Import Data
        try
            Data = textscanu(File,'UTF16-LE',32,13);
            X=str2double(Data(:,1));
            Y=str2double(Data(:,2));
                        
        catch
            Data=importdata(File);
            X=Data.data(:,1);
            Y=Data.data(:,2);
        end
        
        % Clip Data
        if exist('bounds','var')
            b=min(bounds);
            e=max(bounds);
            
            [~,beginning]=min(abs(X-b));
            [~,ending]=min(abs(X-e));
            X=X(beginning:ending);
            Y=Y(beginning:ending);
        else
            b=min(X);
            e=max(X);
        end
        
        % Plot
        figure; hold on
        plot(X,Y);
        set(gca,'xlim',[b e])
        title(File,'Interpreter','none')
        
    end
function C = textscanu(filename, encoding, del_sym, eol_sym)
    
    % C = textscanu(filename, encoding) reads Unicode
    % strings from a file and outputs a cell array of strings.
    %
    % Syntax:
    % -------
    % filename - string with the file's name and extension
    %                 example: 'unicode.txt'
    % encoding - encoding of the file
    %                 default: UTF-16LE
    %                 examples: UTF16-LE (little Endian), UTF8.
    %                 See http://www.iana.org/assignments/character-sets
    %                 MS Notepad saves in UTF-16LE ('Unicode'),
    %                 UTF-16BE ('Unicode big endian'), UTF-8 and ANSI.
    % del_sym - column delimitator symbol in ASCII numeric code
    %                 default: 9 (tabulator)
    % eol_sym - end of line delimitator symbol in ASCII numeric code
    %                 default: 13 (carriage return) [Note: line feed=10]
    %
    % Example:
    % -------
    % C = textscanu('unicode.txt', 'UTF8', 9, 13);
    % Reads the UTF8 encoded file 'unicode.txt', which has
    % columns and lines delimited by tabulators, respectively
    % carriage returns.
    %
    % Note:
    % -------
    % Matlab's textscan function doesn't seem to handle
    % properly multiscript Unicode files. Characters
    % outside the ASCII range are given the \u001a or
    % ASCII 26 value, which usually renders on the
    % screen as a box.
    %
    % Additional information at "Loren on the Art of Matlab":
    % http://blogs.mathworks.com/loren/2006/09/20/
    % working-with-low-level-file-io-and-encodings/#comment-26764
    %
    % Bug:
    % -------
    % When inspecting the output with the Array Editor,
    % in the Workspace or through the Command Window,
    % boxes might appear instead of Unicode characters.
    % Type C{1,1} at the prompt: you will see the correct
    % string. Also: in Array Editor click on C then C{1,1}.
    %
    % Matlab version: starting with R2006b
    %
    % Revisions:
    % -------
    % 2008.02.27 - function creation
    %
    % Created by: Vlad Atanasiu / atanasiu@alum.mit.edu
    
    switch nargin
        case 4
        case 3
            eol_sym = 13;
        case 2
            eol_sym = 13;   % end of line symbol (CR=13, LF=10)
            del_sym = 9;    % column delimitator symbol (TAB=9)
        case 1
            eol_sym = 13;
            del_sym = 9;
            encoding = 'UTF16-LE';
    end
    warning off MATLAB:iofun:UnsupportedEncoding;
    
    % read input
    fid = fopen(filename, 'r', 'l', encoding);
    S = fscanf(fid, '%c');
    fclose(fid);
    
    % remove Byte Order Marker and add an
    % end of line mark at the end of the file
    S = [S(2:end) char(eol_sym)];
    
    % locates column delimitators and end of lines
    del = find(abs(S) == del_sym);
    eol = find(abs(S) == eol_sym);
    
    % get number of rows and columns in input
    row = numel(eol)-1;
    col = 1 + numel(del) / row;
    C = cell(row,col); % output cell array
    
    % catch errors in file
    if col - fix(col) ~= 0
        error(['Error: The file has an odd number of columns ',...
            'or line ends are malformed.'])
    end
    
    m = 1;
    n = 1;
    sos = 1;
    
    % parse input
    if col == 1
        % single column input
        for r = 1:row
            eos = eol(n) - 1;
            C(r,col) = {S(sos:eos)};
            n = n + 1;
            sos = eos + 3;
        end
    else
        % multiple column input
        for r = 1:row
            for c = 1:col-1
                eos = del(m) - 1;
                C(r,c) = {S(sos:eos)};
                sos = eos + 2;
                m = m + 1;
            end
            % last string in the row
            sos = eos + 2;
            eos = eol(n) - 1;
            C(r,col) = {S(sos:eos)};
            n = n + 1;
            sos = eos + 3;
        end
    end
