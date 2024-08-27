function Output=getHPLCAreas(varargin)
    
    % getHPLC calculates the peak areas of high-performance liquid chromatography
    % data. We assume that each sample has its own folder containing all
    % relevant EIC files.We assume that the files are .txt files.
    
    % Inputs:
    % varargin specifies optional parameter name/value pairs.
    % Folders - names of folders to calculate indices for. For one folder,
    %           this input can be a single string. To calculate the areas
    %           for multiple files at once, use a cell array of strings,
    %           e.g. {'ARCTICSTD';'RR B1-A1'}
    %           Default: All folders in the current directory
    % Files - names of individual files to calculate areas for. The files
    %           should be a cell array of strings,
    %           e.g. {'654.txt';'744.txt';'1302.txt';'1300.txt';...
    %               '1298.txt';'1296.txt';'1292.txt';'1050.txt';...
    %               '1048.txt';'1036.txt';'1034.txt';'1022.txt';...
    %               '1020.txt';'1018.txt'};
    %           Default: All files in the folder that have a retention time
    %           defined in the variable 'RT_Sample'
    % RT - Estimated retention times of the peaks of interest in each file.
    %           Default: {16.85;21.54;20.5;19.44;[15.29,16.60];13.84;12.72;11.58;10.66},
    %                     corresponding to the files:
    %               {'744';'1022';'1036';'1050';'1292';'1296';'1298';'1300';'1302'}
    %           Note: You can change this default value to match your mass
    %           spec values by changing the variable 'RT_Sample'
    % ShowGraphs - 0 or 1, where 0 indicates that you do not want to
    %           produce graphs that show the EIC trace and the Gaussian fit
    %           Default: 0
    % Window -  Peaks must fall within a window of length w that is
    %           centered around the retention times defined by RT. To
    %           improve the speed of the calculation, make the window
    %           smaller. To improve the accuracy, you may need to increase
    %           the window size.
    %           Default: 5
    
    % Output:
    % Folders - List of the folders included, as input
    % Files - List of the files included, as input
    % InputTimes - List of the retention times, as input
    % PeakTimes - Retention times of each peak for each folder
    % PeakAreas - Areas of each peak for each folder
    
    
    % Example: Output=getHPLCAreas;
    % Example: Output=getHPLCAreas('Folders','ARC_1','Files',{'1022.txt';'1036.txt';'1050.txt';'1292.txt';'1296.txt';'1298.txt';'1300.txt';'1302.txt'},'RT',{21.54; 20.5; 19.44; [15.29,16.60];  13.84; 12.72; 11.58; 10.66});
    % Example: Output=getHPLCAreas('Folders',{'ARC_1';'OCE205_1'},'Files',{'1022.txt','1036.txt','1050.txt','1292.txt','1296.txt','1298.txt','1300.txt','1302.txt'},'RT',{{21.54, 20.5, 19.44, [15.29,16.60],  13.84, 12.72, 11.58, 10.66}; {21.98, 21.006, 20.002, [15.623,16.934], 14.144 ,13.001, 11.829, 10.881}},'Window',5);

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %         Sort Inputs           %
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    % Define Defaults
    options = struct('folders',[],'files',[],'showgraphs',0,...
        'window',5,'rt',[]);
    
    % Read the Acceptable Names
    optionNames = fieldnames(options);
    
    % Count arguments
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
        error('getPeaks needs propertyName/propertyValue pairs')
    end
    
    listedrt=0;
    listedfiles=0;
    newrt=0;
    for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
        inpName = lower(pair{1}); % make case insensitive
        
        if strcmpi(inpName,'rt')
            listedrt=1;
        end
        if strcmpi(inpName,'files')
            listedfiles=1;
        end
        if any(strcmpi(inpName,optionNames))
            % Overwrite Options
            options.(inpName) = pair{2};
        else
            error('%s is not a recognized parameter name',pair{1})
        end
    end
    
    % Folders
    Folders=options.folders;
    if isempty(Folders)
        % Find All Folders in Current Directory
        d=dir;
        Folders={d([d(:).isdir]).name};
        Folders=Folders(3:end); % Remove . and ..
        
    elseif ischar(Folders)
        % Use Only Specified Folders
        Folders=cellstr(Folders);
    end
    
    % Files
    Files=options.files;
    if ischar(Files)
        Files={Files};
    end
    if listedfiles
        
        if ~iscell(Files);
            Files=num2cell(Files);
        end
        p1=size(Files);
        if p1==1
            Files=Files';
        end
        if size(Files,1)~=length(Folders)
            Files=repmat({Files},length(Folders),1);
        end
        for f=1:length(Files)
            if ~iscell(Files{f})
                Files{f}=Files(f);
            end
        end
    end
    OriginalFiles=Files;
    
    % Retention Times
    RT=options.rt;
    if listedrt
        if ~iscell(RT);
            RT=num2cell(RT);
        end
        [p1,~]=size(RT);
        if p1==1
            RT=RT';
        end
        if size(RT,1)~=length(Folders)
            RT=repmat({RT},length(Folders),1);
        end
        for f=1:length(RT)
            if ~iscell(RT{f})
                RT{f}=RT(f);
            end
        end
        if sum(cellfun(@numel,RT))<2
            pad_rt=1;
        else
            pad_rt=0;
        end
    else
        pad_rt=0;
    end
    
    OriginalRT=RT;
    
    
    
    % Other
    ShowGraphs=options.showgraphs;
    Window=options.window;
    
    CD=cd;
    
    % Preallocate Arrays
    Left=cell(length(Folders),1);
    Right=cell(length(Folders),1);
    XGauss=cell(length(Folders),1);
    YGauss=cell(length(Folders),1);
    AUC=cell(length(Folders),1);
    PeakTime=cell(length(Folders),1);
    Height=cell(length(Folders),1);
    TotalTime=cell(length(Folders),1);
    
    Areas=cell(length(Folders),1);
    PeakTimes=cell(length(Folders),1);
    
    % Sample RT based on ARC_1
    % This should be updated periodically to reflect the current state of
    % the mass spec.
    RT_Sample={'744', 16.85;
        '1022', 21.54;
        '1036',20.5;
        '1050',19.44;
        '1292',[15.29,16.60];
        '1296',13.84;
        '1298',12.72;
        '1300',11.58;
        '1302',10.66};
    
    
    for f=1:length(Folders)
        cd([CD '/' Folders{f}]);
        
        if listedfiles
            Files=OriginalFiles{f};
        else
            % Find All Files in Current Directory
            Files=dir('*.txt');
            Files={Files.name};
            Files=Files(~cellfun(@isempty,regexp(Files,'\d'))); % Remove TIC files
            
            for a=1:length(Files)
                Files{a} = Files{a}(1:end-4);
            end
            
            [~,ia,ib]=intersect(Files,RT_Sample(:,1));
            Files=Files(ia);
            RT_Sample=RT_Sample(ib,:);
            
            for a=1:length(Files)
                Files{a} = [Files{a} '.txt'];
            end
            
            if listedrt
                file_num=length(RT);
                if length(Files)~=file_num
                    error('Please include file names to match retention times.');
                end
            else
                RT=RT_Sample(:,2);
                newrt=1;
            end
        end
        
        if ischar(Files)
            % Use Only Specified Folders
            Files=cellstr(Files);
        end
        
        if newrt
            OriginalFiles{f}=Files;
            OriginalRT{f}=RT;
        end
        
        if listedrt
            RT=OriginalRT;
            RT=RT{f};
            
            if ~iscell(RT)
                RT={RT};
            end
        end

        % If fewer than 2 retention times were specified, use other files to
        % pad the list of retention times. This will improve the accuracy of
        % the relative retention times.
        if pad_rt
            n=2-length(OriginalRT); % Number of new files to add
            
            % List all options for new files
            TempFiles=dir('*.txt');
            TempFiles={TempFiles.name};
            TempFiles=TempFiles(~cellfun(@isempty,regexp(TempFiles,'\d'))); % Remove TIC files
            TempFiles=sort_nat(TempFiles);
            
            % Use only temporary files that have a pre-defined RT
            for j=1:length(TempFiles)
                TempFiles{j}=TempFiles{j}(1:end-4);
            end
            TempFiles=intersect(TempFiles,RT_Sample(:,1));
            for j=1:length(TempFiles)
                TempFiles{j}=[TempFiles{j} '.txt'];
            end
            
            % Find the indices of the files that were defined by the user
            first_ind=regexp(TempFiles,Files{1});
            first_ind=~cellfun(@isempty,first_ind);
            first_ind=find(first_ind);
            
            last_ind=regexp(TempFiles,Files{end});
            last_ind=~cellfun(@isempty,last_ind);
            last_ind=find(last_ind);
            
            % Use nearby files to pad the list of retention times
            if last_ind+n<=length(TempFiles)
                NewFiles=TempFiles(last_ind+1:last_ind+n);
            else
                NewFiles=TempFiles(last_ind+1:end);
                n=n-length(NewFiles);
                NewFiles=[TempFiles(first_ind-n:first_ind-1), NewFiles];
            end
            
            OldFiles=Files;
            Files=[OldFiles,NewFiles];
            Files=sort(Files);
            old_file_ind=find(strcmp(Files,OldFiles)); % Indices of the original RTs.
            RT=[]; % Delete old retention times because they will not match up with the new files
            
        end
        
        if ~listedrt || isempty(RT)
            RT=cell(length(Files),1);
            
            for i=1:length(Files)
                file=Files{i};
                gdgt=file(1:regexp(file,'.txt')-1);
                right_num=regexp(RT_Sample(:,1),gdgt);
                right_num=~cellfun(@isempty,right_num);
                if sum(right_num)
                    RT(i)=RT_Sample(right_num,2);
                else
                    error(['Unknown retention time for ' num2str(gdgt) '. Please include retention times as input.']);
                end
            end
            OriginalRT{f}=RT;
        end
        
        for i=1:length(RT)
            if size(RT{i},1) ~= length(RT{i})
                RT{i}=RT{i}';
            end
        end
        
        % Import Data
        Xs=cell(length(Files),1);
        Ys=cell(length(Files),1);
        
        for num=1:length(Files)
            file=Files{num};
            data = textscanu(file,'UTF16-LE',32,13);
            
            X=sscanf(strjoin(data(:,1)'),'%f');
            Y=sscanf(strjoin(data(:,2)'),'%f');
            
            Xs{num}=X;
            Ys{num}=Y;
        end
        
        % Find Relative Retention Times in Files
        FinalPeaks=findPeaks_RTCell(Xs,Ys,RT,Window);
        
        % Remove data from files that were only used to pad the list of RTs
        if pad_rt
            FinalPeaks=FinalPeaks(old_file_ind);
            Xs=Xs(old_file_ind);
            Ys=Ys(old_file_ind);
            Files=Files(old_file_ind);
            RT=RT(old_file_ind);
        end
        
        % Preallocate Cell
        AUC{f}=cell(length(Files),1);
        PeakTime{f}=cell(length(Files),1);
        Height{f}=cell(length(Files),1);
        Left{f}=cell(length(Files),1);
        Right{f}=cell(length(Files),1);
        XGauss{f}=cell(length(Files),1);
        YGauss{f}=cell(length(Files),1);
        
        for i=1:length(Files)
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
            % Baseline Drift Correction and Smoothing %
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
            X=Xs{i};
            Y=Ys{i};
            
            BaseLine=median(Y);
            Y=Y-BaseLine;
            
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
            %              Find Peaks                 %
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
            ind_peaks=zeros(length(FinalPeaks{i}),1);
            
            for p=1:length(FinalPeaks{i})
                [~,ind]=min(abs(X-FinalPeaks{i}(p)));
                ind_peaks(p)=ind;
            end
            
            try
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
                %            Find Boundaries              %
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
                [left, right, xgauss, ygauss]=findboundariesTEX(X, Y, ind_peaks);
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
                %              Find Area                  %
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
                AUC{f}{i,1}=zeros(length(ind_peaks),1);
                PeakTime{f}{i,1}=zeros(length(ind_peaks),1);
                Height{f}{i,1}=zeros(length(ind_peaks),1);
                TotalTime{f}{i,1}=zeros(length(ind_peaks),1);
                
                for j=1:length(ind_peaks)
                    AUC1=trapz(xgauss(left(j):right(j),j),ygauss(left(j):right(j),j)); % Find total Area Under Curve
                    AUC2=trapz(xgauss(left(j):right(j),j),median(ygauss(:,j)).*ones(size(xgauss(left(j):right(j))))); % Find Area under Baseline
                    AUC{f}{i,1}(j)=AUC1-AUC2;
                    AUC{f}{i,1}(j)=round(60*AUC{f}{i,1}(j)); %Convert to minutes
                    
                    PeakTime{f}{i,1}(j)=roundn(X(ind_peaks(j)),-3);
                    
                    Height{f}{i,1}(j)=roundn((Y(ind_peaks(j))-min(ygauss(left(j):right(j),j))),-1);
                    TotalTime{f}{i,1}(j)=xgauss(right(j),j)-xgauss(left(j),j);
                    
                    Left{f}{i,1}(j)=left(j);
                    Right{f}{i,1}(j)=right(j);
                    
                end
                
                XGauss{f}{i,1}=xgauss;
                YGauss{f}{i,1}=ygauss;
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
                %                Plot                     %
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
                if ShowGraphs
                    figure; hold on
                    
                    Colors=[0 .5 1];
                    
                    for j=1:length(ind_peaks)
                        area(xgauss(left(j):right(j),j),ygauss(left(j):right(j),j),'FaceColor',Colors,'EdgeColor',Colors); % Fill in AUC
                    end
                    
                    plot(X,Y,'k');
                    title([Folders{f} ' ' Files{i}],'Interpreter','None')
                    %Change y-axis limits
                    yl=get(gca,'ylim');
                    set(gca,'ylim',[0 yl(2)]);
                    
                    %Override scientific notation for the x-axis
                    n=get(gca,'YTick');
                    set(gca,'YTickLabel',strread(num2str(n),'%s'));
                    set(gca,'Layer','top');
                    hold off
                end
                
            catch err
                save('Error.mat','err')
                disp(['Folder: ' Folders{f} '; File: ' Files{i}])
                disp(err.message);
                
                AUC{f}{i,1}=NaN*ones(length(ind_peaks),1);
                PeakTime{f}{i,1}=NaN*ones(length(ind_peaks),1);
                Height{f}{i,1}=NaN*ones(length(ind_peaks),1);
                TotalTime{f}{i,1}=NaN*ones(length(ind_peaks),1);
                
                XGauss{f}{i,1}=NaN;
                YGauss{f}{i,1}=NaN;
                
                Left{f}{i,1}=NaN;
                Right{f}{i,1}=NaN;
                
                
                if ShowGraphs
                    figure; hold on
                    
                    plot(X,Y,'k');
                    %Change y-axis limits
                    yl=get(gca,'ylim');
                    set(gca,'ylim',[0 yl(2)]);
                    title(num2str(RT{i}))
                    
                    %Override scientific notation for the x-axis
                    n=get(gca,'YTick');
                    set(gca,'YTickLabel',strread(num2str(n),'%s'));
                    set(gca,'Layer','top');
                    hold off
                end
                
            end
        end
        
        Areas{f}=vertcat(AUC{f}{:}); % Convert to double array
        PeakTimes{f}=vertcat(PeakTime{f}{:}); % Convert to double array
        
    end
    cd(CD)
    
    % Create Double Array of Data
    RT=OriginalRT;
    Files=OriginalFiles;
    N_RTs=0;
    for f=1:length(Folders)
        N_RTs=max(N_RTs, sum(cellfun(@numel,RT{f})));
    end
    
    AllRTs=NaN*ones(length(Folders),N_RTs);
    PeakAreas=NaN*ones(length(Folders),N_RTs);
    PeakTimes=NaN*ones(length(Folders),N_RTs);
    AllFiles = cell(length(Folders),N_RTs);
    
    for f=1:length(Folders)
        AllRTs(f,1:sum(cellfun(@numel,RT{f}))) = horzcat(RT{f}{:});
        PeakAreas(f,1:sum(cellfun(@numel,RT{f})))=vertcat(AUC{f}{:});
        PeakTimes(f,1:sum(cellfun(@numel,RT{f})))=vertcat(PeakTime{f}{:});
        
        num_areas=cellfun(@numel,AUC{f});
        tempFiles=cell(sum(num_areas),1);
        count=0;
        for i=1:length(Files{f})
            for j=1:num_areas(i)
                count=count+1;
                Period=regexp(Files{f}{i},'\.');
                tempFiles{count}=Files{f}{i}(1:(Period-1)); % Remove extension
            end
        end
        
        AllFiles(f,1:sum(cellfun(@numel,RT{f}))) = tempFiles;
    end
    
    % Sort Output alphabetically
    [~,sort_index]=sort_nat(Folders);
    Folders=Folders(sort_index);
    AllFiles=AllFiles(sort_index,:);
    AllRTs=AllRTs(sort_index,:);
    PeakAreas=PeakAreas(sort_index,:);
    PeakTimes=PeakTimes(sort_index,:);
    
    % If all folders use the same files, include only one listing of the
    % files
    for f=1:length(Folders)
        A=AllFiles(1,:);
        B=AllFiles(f,:);
        
        if sum(cellfun(@isempty,B))
            break
        end
        
        if ~isempty(setdiff(A,B))
            break
        end
    end
    if f==length(Folders)
        AllFiles=AllFiles(1,:);
        
        % If all folders use the same retention times, include only one
        % listing of the retention times
        if ~sum(diff(AllRTs))
            AllRTs=AllRTs(1,:);
        end
    end
    
    Output.Folders=Folders;
    Output.Files=AllFiles;
    Output.InputTimes=AllRTs;
    Output.PeakTimes=PeakTimes;
    Output.Areas=PeakAreas;
function FinalPeaksCell=findPeaks_RTCell(Xs,Ys,RT,Window)
    % This function finds the peaks based on the retention time.
    % Xs is a cell array. Each cell contains the x-axis data from the files
    % you are interested in
    % Ys is a cell array. Each cell contains the y-axis data from the
    % files.
    % RT is a cell array. Each cell contains the an array of the retention
    % times you are interested in for that file.
    if ~exist('Window','var')
        Window=(max(Xs{1})-min(Xs{1}))./10;
    end
    
    Peaks=cell(length(Xs),1);
    Heights=cell(length(Xs),1);
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %   Find All Possible Peaks  %
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
    for i=1:length(Xs)
        X=Xs{i};
        Y=Ys{i};
        
        % Find peaks in smoothed data. Peaks must be separated by at least 0.1
        % seconds.
        
        [~,ind_peaks_smooth]=findpeaks(smooth(Y),'MinPeakDistance',ceil(0.1./mean(diff(X))));
        [~,ind_peaks]=findpeaks(Y,'MinPeakDistance',ceil(0.1./mean(diff(X))));
        
        for r=1:length(ind_peaks_smooth)
            [~,rightpeak]=min(abs(ind_peaks - ind_peaks_smooth(r)));
            ind_peaks_smooth(r)=ind_peaks(rightpeak);
        end
        
        ind_peaks=ind_peaks_smooth;
        
        % Finds all distinct peaks that are withing the window of each RT.
        peaks=[];
        
        for r=1:length(RT{i})
            
            minWindow=RT{i}(r)-(Window/2);
            maxWindow=RT{i}(r)+(Window/2);
            
            [~,a]=min(abs(X-minWindow));
            [~,b]=min(abs(X-maxWindow));
            
            % Fit a Baseline
            f1 = fit(X(a:b),Y(a:b),'poly1', 'Exclude', Y(a:b)>median(Y(a:b)));
            Median=f1.p1*X+f1.p2;
            
            % Peaks must be within the window about a certain height above
            % baseline
            tempPeaks=ind_peaks;
            tempPeaks=tempPeaks(tempPeaks>a & tempPeaks<=b);
            tempPeaks=tempPeaks(Y(tempPeaks) > 0.02*(max(Y(a:b) - Median(a:b)))+Median(tempPeaks));
            tempPeaks=tempPeaks(Y(tempPeaks) > 0.02*(max(Y - Median))+Median(tempPeaks));
            
            peaks=[peaks;tempPeaks];
        end
        
        peaks=peaks(X(peaks)>(min(RT{i})-Window));
        peaks=peaks(X(peaks)<(max(RT{i})+Window));
        peaks=unique(peaks);
        
        % If there aren't enough peaks for the number of retention times,
        % lower the height threshold until there are enough peaks.
        multiple=0.02;
        while length(peaks)<length(RT{i})
            multiple=multiple-0.005;
            for r=1:length(RT{i})
                
                
                minWindow=RT{i}(r)-(Window/2);
                maxWindow=RT{i}(r)+(Window/2);
                
                [~,a]=min(abs(X-minWindow));
                [~,b]=min(abs(X-maxWindow));
                
                
                
                % Peaks must be within the window about a certain height above
                % baseline
                tempPeaks=ind_peaks;
                tempPeaks=tempPeaks(tempPeaks>a & tempPeaks<=b);
                tempPeaks=tempPeaks(Y(tempPeaks) > multiple*(max(Y(a:b) - Median(a:b)))+Median(tempPeaks));
                tempPeaks=tempPeaks(Y(tempPeaks) > multiple*(max(Y - Median))+Median(tempPeaks));
                
                peaks=[peaks;tempPeaks];
            end
            
            peaks=peaks(X(peaks)>(min(RT{i})-Window));
            peaks=peaks(X(peaks)<(max(RT{i})+Window));
            peaks=unique(peaks);
            
        end
        
        % If there are too many peaks, take only the highest peaks. This
        % reduces the calculation time.
        if length(peaks)>5*length(RT{i})
            
            [~,ind_sorted]=sort(Y(peaks),'descend');
            peaks=peaks(ind_sorted(1:5*length(RT{i})));
        end
        
        % Find Critical Points
        [critical_points]=findCritPoints(X,Y,RT{i},Window);
        peaks=[peaks; critical_points];
        peaks=unique(peaks);
        
        Peaks{i,1}=X(peaks);
        Heights{i,1}=Y(peaks);
        
    end
    
    % Find all peaks with the same retention time. For repeating retention
    % times, keep only the highest peak.
    peaks=vertcat(Peaks{:});
    heights=Heights;
    Heights=vertcat(Heights{:});
    
    [peaks, is]=sort(peaks);
    Heights=Heights(is);
    
    i=1;
    while i<length(peaks)
        value=peaks(i);
        repeats=find(peaks==value);
        max_height=max(Heights(repeats));
        
        peaks(repeats(Heights(repeats) ~= max_height))=[];
        Heights(repeats(Heights(repeats) ~= max_height))=[];
        
        repeats=find(peaks==value);
        if length(repeats)>1
            peaks(repeats(2:end))=[];
            Heights(repeats(2:end))=[];
        end
        
        i=i+1;
        
    end
    
    % Find Difference in RTs
    tempRT=vertcat(RT{:});
    tempRT=sort(tempRT);
    DiffRT=diff(tempRT);
    
    % Find all possible combinations
    
    % Sort the retention time
    tempRT=vertcat(RT{:});
    [~, FullIX]=sort(tempRT);
    [~,RevIX]=sort(FullIX);
    
    tempRT=tempRT(FullIX);
    Choices=cellfun(@numel,RT);
    
    Options_Combos=cell(size(Peaks));
    Temp=cell(size(Peaks));
    for i=1:size(Peaks,1)
        Options_Combos{i} = nchoosek(Peaks{i},Choices(i));
        Temp{i}=1:size(Options_Combos{i},1);
    end
    
    while prod(cellfun(@numel,Temp))>10000000
        
        [~,ind_max]=max(cellfun(@numel,Temp));
        [~,ind_min]=min(heights{ind_max});
        
        Peaks{ind_max}(ind_min)=[];
        heights{ind_max}(ind_min)=[];
        Options_Combos{ind_max} = nchoosek(Peaks{ind_max},Choices(ind_max));
        Temp{ind_max}=1:size(Options_Combos{ind_max},1);
        
    end
    
    CombOfCombs=allcomb(Temp{:});
    
    combinations=zeros(size(CombOfCombs));
    count=0;
    for i=1:size(Peaks,1)
        
        combinations(:,count+1:count+size(RT{i},1))=Options_Combos{i}(CombOfCombs(:,i),:);
        count=count+size(RT{i},1);
    end
    
    combinations=combinations(:,FullIX);
    
    if min(diff(tempRT))>1
        % Remove combinations with out-of-order retention times
        [~, ISSORTED]=sort(combinations,2);
        ISSORTED=floor(sum(ISSORTED==repmat(1:length(tempRT),size(combinations,1),1),2)./size(combinations,2));
        ISSORTED=logical(ISSORTED);
        combinations=combinations(ISSORTED,:);
        
        % Remove combinations with repeating retention times
        temp_comb=sort(combinations,2);
        diff_combinations=diff(temp_comb,1,2);
        [r,~]=find(diff_combinations==0);
        r=unique(r);
        combinations(r,:)=[];
    end
    
    if isempty(combinations)
        error('No possible combinations. Increase window size.');
    end
    
    % Find the difference between the relative retention times for each
    % combnination and the defined relative retention times. Choose the
    % combination that has the minimum sum of the absolute values of these
    % retention times.
    DiffPeakComb=diff(combinations,1,2);
    temp=DiffPeakComb-repmat(DiffRT',size(combinations,1),1);
    temp=temp.^2;
    temp=sum(temp,2);
    
    [~,best]=min(temp);
    combinations=combinations(best,:);
    
    % Sort retention times into the original order and format
    FinalPeaks = combinations;
    FinalPeaks=FinalPeaks(RevIX);
    
    FinalPeaksCell=RT;
    count=1;
    for r=1:size(RT,1)
        FinalPeaksCell{r}=FinalPeaks(count:count+length(RT{r})-1);
        count=count+length(RT{r});
    end
function [left, right, XGAUSS, YGAUSS, smoothY, type, gof, A, B]=findboundariesTEX(X, Y, ind_peaks)
    
    % This function finds the boundaries on the left and right side of each
    % peak.
    
    right=zeros(length(ind_peaks),1);
    left=zeros(length(ind_peaks),1);
    
    N=10000;
    XGAUSS=zeros(N,length(ind_peaks));
    YGAUSS=zeros(N,length(ind_peaks));
    
    % Compute Derivatives
    smoothY=perfectSmoother(Y,1);
    velocity=forward_derivative(X,smoothY);
    velocity=velocity./max(abs(velocity));
    
    smooth_velo=perfectSmoother(velocity,1);
    
    [~,ind_peaks_sort]=sort(ind_peaks);
    ind_peaks=ind_peaks(ind_peaks_sort);
    unsort=zeros(size(ind_peaks));
    for i=1:length(unsort)
        unsort(i)=find(ind_peaks_sort==i);
    end
    
    for j=1:length(ind_peaks)
        try
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
            %           Fit Only The Main Peak                 %
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
            dt=ceil(0.025./mean(diff(X)));
            
            % 'A' is the critical point to the left of the peak
            % 'B' is the critical point to the right of the peak
            % 'C' are the local minima in the trace
            % 'D' is a left-hand shoulder
            % 'E' is a right-hand shoulder
            
            % Find Critical Point
            A=find(smooth_velo(1:ind_peaks(j)-2*dt)<0.005,1,'last');
            B=find(smooth_velo(ind_peaks(j)+ 2*dt :end)>-0.005,1,'first');
            
            if isempty(A)
                A=1;
            else
                A=A+1;
            end
            if isempty(B)
                B=length(X);
            else
                B=B+(ind_peaks(j)+2*dt -1);
            end
            
            smooth_velo=smooth_velo./max(abs(smooth_velo(A:B)));
            [~,C]=findpeaks(-smoothY,'threshold',1);
            leftC=C(C<ind_peaks(j));
            rightC=C(C>ind_peaks(j));
            
            % Find local mins in velocity to the left of the peak.
            [~,D]=findpeaks(-smooth_velo);
            D=D(smooth_velo(D)>0.015);
            [~,ind]=max(smooth_velo(A:ind_peaks(j)));
            ind=ind+A;
            D=D(D<ind);
            
            % Find local maxes in velocity to the left of the peak.
            [~,E]=findpeaks(smooth_velo);
            E=E(smooth_velo(E)<-0.015);
            [~,ind]=min(smooth_velo(ind_peaks(j):B));
            ind=ind+ind_peaks(j);
            E=E(E>ind);
            
            % 'A' is the left-hand boundary to be fit.
            % 'B' is the right-hand boundary to be fit.
            A=max([leftC; A; D]);
            B=min([rightC; B; E]);
            
            [fitresult, gof, type] = GaussFit(X(A:B), smoothY(A:B));
            
            a=fitresult.a1;
            b=fitresult.b1;
            c=fitresult.c1;
            
            XGauss=linspace(min(X),max(X),N);
            YGauss=a*exp(-((XGauss-b)/c).^2);
            
            RSquare=gof.rsquare;
            
            % If the fit is bad, re-fit with stricter boundaries
            if max(YGauss)>1.05*Y(ind_peaks(j)) || RSquare<0.9
                
                Lower=[0 0 0];
                Upper=[1.05*Y(ind_peaks(j)) max(X(A:B)) 1];
                StartPoint = [Y(ind_peaks(j)) X(ind_peaks(j)) 0.05];
                
                [fitresult, gof, type] = GaussFit(X(A:B), smoothY(A:B),Lower,Upper,StartPoint);
                
                a=fitresult.a1;
                b=fitresult.b1;
                c=fitresult.c1;
                
                XGauss=linspace(min(X),max(X),N);
                YGauss=a*exp(-((XGauss-b)/c).^2);
                
            end
            
            if strcmp(type,'asymmetric')
                disp(['A peak in ' file ' is asymmetrical. Consider hand-picking this trace.'] )
            end
            
            [~,ind_gausspeaks]=findpeaks(YGauss);
            [~,best_peak]=min(abs(X(ind_peaks(j)) - XGauss(ind_gausspeaks)));
            ind_gausspeaks=ind_gausspeaks(best_peak);
            
            % Find Where Peak Reaches Baseline
            MinBoundHeight=0.00001*YGauss(ind_gausspeaks); % Lower cut off point is 1% of peak height above baseline
            
            Y1=YGauss(1:ind_gausspeaks)-MinBoundHeight;
            leftgauss=find(Y1<0,1,'last');
            
            if isempty(leftgauss)
                leftgauss=1;
            end
            
            Y1=YGauss(ind_gausspeaks:end)-MinBoundHeight;
            rightgauss=find(Y1<0,1,'first') + ind_gausspeaks-1;
            if isempty(rightgauss)
                rightgauss=length(YGauss);
            end
            
            left(j)=leftgauss;
            right(j)=rightgauss;
            
            XGAUSS(:,j)=XGauss;
            YGAUSS(:,j)=YGauss;
            
        catch err
            disp(err.message)
            left(j)=NaN;
            right(j)=NaN;
            
            XGAUSS(:,j)=NaN;
            YGAUSS(:,j)=NaN;
        end
        
    end
    
    left=left(unsort);
    right=right(unsort);
    XGAUSS=XGAUSS(:,unsort);
    YGAUSS=YGAUSS(:,unsort);
function [fitresult, gof, type] = GaussFit(X, Y, Lower, Upper,StartPoint)
    
    % Set up fittype and options.
    [xData, yData] = prepareCurveData( X, Y );
    
    % Find the minimum height that is reached on both sides of the peak.
    % Use only the values above this minimum height on both sides to test the asymmetry.
    Height=max(Y);
    
    [~,peak]=max(Y);
    left_min_val=min(Y(1:peak));
    right_min_val=min(Y(peak:end));
    
    leftHeight=(Height-left_min_val)./Height;
    rightHeight=(Height-right_min_val)./Height;
    
    if left_min_val>right_min_val
        temp=Y(peak:end);
        right_min=find(temp<left_min_val,1,'first');
        
        testing_X=X(1:(peak+right_min-1));
        testing_Y=Y(1:(peak+right_min-1));
    else
        temp=Y(1:peak);
        left_min=find(temp<right_min_val,1,'last');
        
        testing_X=X(left_min:end);
        testing_Y=Y(left_min:end);
    end
    
    % Determine asymmetry of peak by comparing velocity on both sides of
    % peak
    [~,peak]=max(testing_Y);
    
    try
        velocity=forward_derivative(testing_X,testing_Y);
        velocity=velocity./max(abs(velocity));
        
        leftvelo=max(velocity(1:peak));
        rightvelo=max(-velocity(peak:end));
    catch
        leftHeight=0;
        rightHeight=0;
        
        velocity=forward_derivative(X,Y);
        velocity=velocity./max(abs(velocity));
        
        leftvelo=max(velocity(1:peak));
        rightvelo=max(-velocity(peak:end));
    end
    
    if (leftHeight<0.2) || (rightHeight<0.2)% Assume Symmetry (Small Peak)
        type='symmetric';
    elseif (leftvelo<0.25)
        type='asymmetric';
    elseif(rightvelo<0.25)
        type='asymmetric';
    else
        type='symmetric';
    end
    
    ft = fittype('gauss1');
    
    % Set boundaries for coefficients
    if ~exist('Lower','var')
        Lower=[0 0 0];
    end
    if ~exist('Upper','var')
        Upper = [2*max(Y) 2*max(X) 5];
    end
    if ~exist('StartPoint','var')
        [~,ind]=max(Y);
        StartPoint = [Y(ind) X(ind) 0.05];
    end
    
    % Fit model to data.
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'notify';
    opts.StartPoint=StartPoint;
    opts.Upper=Upper;
    opts.Lower=Lower;
    
    [fitresult, gof] = fit( xData, yData, ft, opts );
function [critical_points]=findCritPoints(X,Y,RT,Window)
    
    % Find Critical Points
    smoothY=perfectSmoother(Y,15);
    velocity=forward_derivative(X,smoothY);
    velocity=velocity./max(abs(velocity));
    velocity=perfectSmoother(velocity,30);
    
    minWindow=min(RT)-(Window/2);
    maxWindow=max(RT)+(Window/2);
    
    [~,a]=min(abs(X-minWindow));
    [~,b]=min(abs(X-maxWindow));
    
    velocity=velocity./max(abs(velocity(a:b)));
    
    criticalps=find(abs(velocity)<0.01);
    clusters=findclusters(criticalps,1);
    clusters=criticalps(clusters);
    
    critical_points=zeros(size(clusters,1),1);
    SameSign=zeros(size(clusters,1),1);
    
    for i=1:size(clusters,1)
        Sign=sign(velocity(clusters(i,1):clusters(i,2)));
        SameSign(i)=sum(diff(Sign)~=0);
        
        [~,m]=min(abs(velocity(clusters(i,1):clusters(i,2))));
        critical_points(i)=m+clusters(i,1)-1;
        
    end
    
    for i=1:size(critical_points,1)
        if SameSign(i)>1
            critical_points(i)=0;
            SameSign(i)=0;
        end
    end
    critical_points(critical_points==0)=[];
    
    f1 = fit(X,Y,'poly1', 'Exclude', Y>median(Y));
    Median=f1.p1*X+f1.p2;
    critical_points=critical_points(Y(critical_points)>Median(critical_points));
    
    [~,mins]=findpeaks(-smoothY,'MinPeakDistance',ceil(0.1./mean(diff(X))));
    [~,maxs]=findpeaks(smoothY,'MinPeakDistance',ceil(0.1./mean(diff(X))));
    
    % Remove any critical points that are close to a peak. They are
    % describing the same option.
    dt=ceil(0.2./mean(diff(X)));
    p=1;
    while p<=length(critical_points)
        DiffCrit1=abs(mins-critical_points(p));
        DiffCrit2=abs(maxs-critical_points(p));
        if min([DiffCrit1; DiffCrit2])<dt
            critical_points(p)=[];
        else
            p=p+1;
        end
    end
    
    critical_points=critical_points(critical_points>a & critical_points<b);
    critical_points=sort(critical_points);
    
    % Remove any critical points that are close to other critical points.
    dt=ceil(0.03./mean(diff(X)));
    p=1;
    while p<length(critical_points)
        DiffCrit=abs(critical_points-critical_points(p));
        DiffCrit=DiffCrit(DiffCrit>0);
        if min(DiffCrit)<=dt
            
            [~,maxY]=max(Y(critical_points(p) : critical_points(p+1)));
            critical_points(p)=critical_points(p)+maxY;
            critical_points(p+1)=[];
            p=p+1;
        else
            p=p+1;
        end
    end
function clusters=findclusters(X,p)
    % This function finds clusters in the array X. A cluster is defined as
    % a subset of X where each element is no more than p away from the
    % next.
    % 'clusters' is a nX2 array, where n is the number of clusters. The
    % first column contains the starting index of each cluster. The second
    % column contains the ending index of each cluster.
    % Ex: X=[1,1.5,3,5,5.7,5.8]; p = 1;
    % clusters = [1,2;3,3;4,6];
    % X(clusters(1,1):clusters(1,2)) = [1,1.5];
    % X(clusters(2,1):clusters(2,2)) = 3;
    % X(clusters(3,1):clusters(3,2)) = [5,5.7,5.8];
    
    % Sort Inputs
    % Default of 'p' is 1.
    if ~exist('p','var')
        original_p=1;
    else
        original_p=p;
    end
    
    % Find the number of clusters
    n_clust=1;
    p=original_p;
    
    while p<=length(X)-1
        
        num=1;
        while (p+num)<=length(X) && ...
                (X(p+num) - X(p+num-1))<=original_p
            
            num=num+1;
        end
        
        n_clust=n_clust+1;
        p=p+num;
    end
    
    % Find the indices of the start and end of each clusters.
    clusters=zeros(n_clust-1,2);
    p=original_p;
    n_clust=1;
    
    while p<=length(X)-1
        
        num=1;
        while (p+num)<=length(X) && ...
                (X(p+num) - X(p+num-1))<=original_p
            
            num=num+1;
        end
        
        clusters(n_clust,1)=p;
        clusters(n_clust,2)=p+num-1;
        
        n_clust=n_clust+1;
        p=p+num;
    end
function [FD1, X1, FD2, X2, FD3, X3, FD4, X4]=forward_derivative(x,y)
    % This function estimates the derivatives of the function y(x).
    
    % First Derivative
    % f'(x_k) = (f(x_k_next) - f(x_k))/(x_k_next - x_k)
    FD1 = diff(y)./diff(x);
    X1=x(1:end-1);
    
    % Second Derivative
    %  f''(x_k) = (f'(x_k_next) - f'(x_k))/(x_k_next - x_k)
    FD2 = diff(FD1)./diff(X1);
    X2=X1(1:end-1);
    
    % Third Derivative
    %  f'''(x_k) = (f''(x_k_next) - f''(x_k))/(x_k_next - x_k)
    FD3 = diff(FD2)./diff(X2);
    X3=X2(1:end-1);
    
    % Fourth Derivative
    %  f'''(x_k) = (f''(x_k_next) - f''(x_k))/(x_k_next - x_k)
    FD4 = diff(FD3)./diff(X3);
    X4=X3(1:end-1);
function A = allcomb(varargin)
    
    % ALLCOMB - All combinations
    %    B = ALLCOMB(A1,A2,A3,...,AN) returns all combinations of the elements
    %    in the arrays A1, A2, ..., and AN. B is P-by-N matrix is which P is the product
    %    of the number of elements of the N inputs. This functionality is also
    %    known as the Cartesian Product. The arguments can be numerical and/or
    %    characters, or they can be cell arrays.
    %
    %    Examples:
    %       allcomb([1 3 5],[-3 8],[0 1]) % numerical input:
    %       % -> [ 1  -3   0
    %       %      1  -3   1
    %       %      1   8   0
    %       %        ...
    %       %      5  -3   1
    %       %      5   8   1 ] ; % a 12-by-3 array
    %
    %       allcomb('abc','XY') % character arrays
    %       % -> [ aX ; aY ; bX ; bY ; cX ; cY] % a 6-by-2 character array
    %
    %       allcomb('xy',[65 66]) % a combination
    %       % -> ['xA' ; 'xB' ; 'yA' ; 'yB'] % a 4-by-2 character array
    %
    %       allcomb({'hello','Bye'},{'Joe', 10:12},{99999 []}) % all cell arrays
    %       % -> {  'hello'  'Joe'        [99999]
    %       %       'hello'  'Joe'             []
    %       %       'hello'  [1x3 double] [99999]
    %       %       'hello'  [1x3 double]      []
    %       %       'Bye'    'Joe'        [99999]
    %       %       'Bye'    'Joe'             []
    %       %       'Bye'    [1x3 double] [99999]
    %       %       'Bye'    [1x3 double]      [] } ; % a 8-by-3 cell array
    %
    %    ALLCOMB(..., 'matlab') causes the first column to change fastest which
    %    is consistent with matlab indexing. Example:
    %      allcomb(1:2,3:4,5:6,'matlab')
    %      % -> [ 1 3 5 ; 1 4 5 ; 1 3 6 ; ... ; 2 4 6 ]
    %
    %    If one of the arguments is empty, ALLCOMB returns a 0-by-N empty array.
    %
    %    See also NCHOOSEK, PERMS, NDGRID
    %         and NCHOOSE, COMBN, KTHCOMBN (Matlab Central FEX)
    
    % for Matlab R2011b
    % version 4.0 (feb 2014)
    % (c) Jos van der Geest
    % email: jos@jasen.nl
    
    % History
    % 1.1 (feb 2006), removed minor bug when entering empty cell arrays;
    %     added option to let the first input run fastest (suggestion by JD)
    % 1.2 (jan 2010), using ii as an index on the left-hand for the multiple
    %     output by NDGRID. Thanks to Jan Simon, for showing this little trick
    % 2.0 (dec 2010). Bruno Luong convinced me that an empty input should
    % return an empty output.
    % 2.1 (feb 2011). A cell as input argument caused the check on the last
    %      argument (specifying the order) to crash.
    % 2.2 (jan 2012). removed a superfluous line of code (ischar(..))
    % 3.0 (may 2012) removed check for doubles so character arrays are accepted
    % 4.0 (feb 2014) added support for cell arrays
    
    error(nargchk(1,Inf,nargin)) ;
    
    NC = nargin ;
    
    % check if we should flip the order
    if ischar(varargin{end}) && (strcmpi(varargin{end},'matlab') || strcmpi(varargin{end},'john')),
        % based on a suggestion by JD on the FEX
        NC = NC-1 ;
        ii = 1:NC ; % now first argument will change fastest
    else
        % default: enter arguments backwards, so last one (AN) is changing fastest
        ii = NC:-1:1 ;
    end
    
    % check for empty inputs
    if any(cellfun('isempty',varargin(ii))),
        warning('ALLCOMB:EmptyInput','Empty inputs result in an empty output.') ;
        A = zeros(0,NC) ;
    elseif NC > 1
        isCellInput = cellfun(@iscell,varargin) ;
        if any(isCellInput)
            if ~all(isCellInput)
                error('ALLCOMB:InvalidCellInput', ...
                    'For cell input, all arguments should be cell arrays.') ;
            end
            % for cell input, we use to indices to get all combinations
            ix = cellfun(@(c) 1:numel(c), varargin,'un',0) ;
            
            % flip using ii if last column is changing fastest
            [ix{ii}] = ndgrid(ix{ii}) ;
            
            A = cell(numel(ix{1}),NC) ; % pre-allocate the output
            for k=1:NC,
                % combine
                A(:,k) = reshape(varargin{k}(ix{k}),[],1) ;
            end
        else
            % non-cell input, assuming all numerical values or strings
            % flip using ii if last column is changing fastest
            [A{ii}] = ndgrid(varargin{ii}) ;
            % concatenate
            A = reshape(cat(NC+1,A{:}),[],NC) ;
        end
    elseif NC==1,
        A = varargin{1}(:) ; % nothing to combine
        
    else % NC==0, there was only the 'matlab' flag argument
        A = zeros(0,0) ; % nothing
    end
function z = perfectSmoother(y,lambda)
    
    m = length(y);
    
    if m>1000
        
        E = speye(m);
        D = diff(E);
        C = chol(E + lambda * D' * D);
        z = C\(C'\y);
    else
        
        E = eye(m);
        D = diff(E);
        z = (E + lambda * D' * D)\y;
        
    end
function [cs,index] = sort_nat(c,mode)
    %sort_nat: Natural order sort of cell array of strings.
    % usage:  [S,INDEX] = sort_nat(C)
    %
    % where,
    %    C is a cell array (vector) of strings to be sorted.
    %    S is C, sorted in natural order.
    %    INDEX is the sort order such that S = C(INDEX);
    %
    % Natural order sorting sorts strings containing digits in a way such that
    % the numerical value of the digits is taken into account.  It is
    % especially useful for sorting file names containing index numbers with
    % different numbers of digits.  Often, people will use leading zeros to get
    % the right sort order, but with this function you don't have to do that.
    % For example, if C = {'file1.txt','file2.txt','file10.txt'}, a normal sort
    % will give you
    %
    %       {'file1.txt'  'file10.txt'  'file2.txt'}
    %
    % whereas, sort_nat will give you
    %
    %       {'file1.txt'  'file2.txt'  'file10.txt'}
    %
    % See also: sort
    
    % Version: 1.4, 22 January 2011
    % Author:  Douglas M. Schwarz
    % Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
    % Real_email = regexprep(Email,{'=','*'},{'@','.'})
    
    
    % Set default value for mode if necessary.
    if nargin < 2
        mode = 'ascend';
    end
    
    % Make sure mode is either 'ascend' or 'descend'.
    modes = strcmpi(mode,{'ascend','descend'});
    is_descend = modes(2);
    if ~any(modes)
        error('sort_nat:sortDirection',...
            'sorting direction must be ''ascend'' or ''descend''.')
    end
    
    % Replace runs of digits with '0'.
    c2 = regexprep(c,'\d+','0');
    
    % Compute char version of c2 and locations of zeros.
    s1 = char(c2);
    z = s1 == '0';
    
    % Extract the runs of digits and their start and end indices.
    [digruns,first,last] = regexp(c,'\d+','match','start','end');
    
    % Create matrix of numerical values of runs of digits and a matrix of the
    % number of digits in each run.
    num_str = length(c);
    max_len = size(s1,2);
    num_val = NaN(num_str,max_len);
    num_dig = NaN(num_str,max_len);
    for i = 1:num_str
        num_val(i,z(i,:)) = sscanf(sprintf('%s ',digruns{i}{:}),'%f');
        num_dig(i,z(i,:)) = last{i} - first{i} + 1;
    end
    
    % Find columns that have at least one non-NaN.  Make sure activecols is a
    % 1-by-n vector even if n = 0.
    activecols = reshape(find(~all(isnan(num_val))),1,[]);
    n = length(activecols);
    
    % Compute which columns in the composite matrix get the numbers.
    numcols = activecols + (1:2:2*n);
    
    % Compute which columns in the composite matrix get the number of digits.
    ndigcols = numcols + 1;
    
    % Compute which columns in the composite matrix get chars.
    charcols = true(1,max_len + 2*n);
    charcols(numcols) = false;
    charcols(ndigcols) = false;
    
    % Create and fill composite matrix, comp.
    comp = zeros(num_str,max_len + 2*n);
    comp(:,charcols) = double(s1);
    comp(:,numcols) = num_val(:,activecols);
    comp(:,ndigcols) = num_dig(:,activecols);
    
    % Sort rows of composite matrix and use index to sort c in ascending or
    % descending order, depending on mode.
    [unused,index] = sortrows(comp);
    if is_descend
        index = index(end:-1:1);
    end
    index = reshape(index,size(c));
    cs = c(index);
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