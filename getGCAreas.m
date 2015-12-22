function Output=getGCAreas(varargin)
    % getGCAreas calculates the peak areas of gas chromatography data. We
    % assume that the samples are .csv files.
    
    % Inputs:
    % varargin specifies optional parameter name/value pairs.
    % Files - names of files to calculate indices for. If there is only one
    %       file, the input can be a string. If there are multiple files,
    %       the input should be a cell array of strings.
    %       Default: All files in the current directory
    %           e.g. {'GC27_1.csv','GC27_3.csv','GC27_4.csv',...
    %                 'KNR197_1.csv','KNR197_2.csv','KNR197_3.csv',...
    %                 'P178_1.csv','P178_2.csv','P178_3.csv'};
    % ?RT? - Estimated retention times of the peaks of interest in each file.  
    %       Default is something like alkenone peaks: [39.75,40.38,43.08,43.34,43.75,43.96]
	%       We recommend plotting a test sample and picking out the RTs of the peaks you want to integrate, then saving those values as a double for input.
    % ShowGraphs - 0 or 1, where 0 indicates that you do not want to
    %           produce graphs that show the EIC trace and the Gaussian fit
    %           Default: 0
    % Window -  Peaks must fall within a window of length w that is
    %           centered around the retention times defined by RT. To
    %           improve the speed of the calculation, make the window
    %           smaller. To improve the accuracy, you may need to increase
    %           the window size.
    %           Default: 2
    % SmoothParam - Amount to smooth the data before integration. Use lower values (e.g., 2) for larger peaks to avoid a reduction in peak height
    %           Default: 15
    
    % Output:
    % Files - List of the files included, as input
    % InputTimes - List of the retention times, as input
    % PeakTimes - Retention times of each peak for each file
    % Areas - Areas of each peak for each file
    
    % Example: Output=getGCAreas;
    % Example: Output=getGCAreas('Files','GC27_2.csv','ShowGraphs',1,'SmoothParam',10);
    % Example: Output=getGCAreas('RT',[36.62, 37.25, 39.78, 40.05, 40.47, 40.65],'ShowGraphs',1);
    % Example: Output=getGCAreas('Files',{'GC27_1.csv','KNR197_1.csv'},'RT',{[36.62, 37.25, 39.78, 40.05, 40.47, 40.65];[33.63, 36.66, 37.19, 39.78, 40.05, 40.42, 40.62]});
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    %         Sort Inputs           %
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    % Define Defaults
    options = struct('files',[],'showgraphs',0,...
        'window',2,'rt',[],'smoothparam',15);
    
    % Read the Acceptable Names
    optionNames = fieldnames(options);
    
    % Count arguments
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
        error('getPeaks needs propertyName/propertyValue pairs')
    end
    
    listedrt=0;
    for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
        inpName = lower(pair{1}); % make case insensitive
        
        if strcmpi(inpName,'rt')
            listedrt=1;
        end
        
        if any(strcmpi(inpName,optionNames))
            % Overwrite Options
            options.(inpName) = pair{2};
        else
            error('%s is not a recognized parameter name',pair{1})
        end
    end
    
    % Files
    Files=options.files;
    if ischar(Files)
        Files={Files};
    end
    
    if isempty(Files)
        % Find All Files in Current Directory
        Files=dir('*.csv');
        Files={Files.name};
    elseif ischar(Files)
        % Use Only Specified Folders
        Files=cellstr(Files);
    end
    
    % Retention Times
    RTs=options.rt;
    % Sample RT based on 1234std3_052414.csv
    % This should be updated periodically to reflect the current state of
    % the mass spec.
    RT_Sample=[39.75,40.38,43.08,43.34,43.75,43.96];
    
    if listedrt
        if ~iscell(RTs);
            RTs=num2cell(RTs);
        end
        [p1,~]=size(RTs);
        if p1==1
            RTs=RTs';
        end
        if size(RTs,1)~=length(Files)
            RTs=repmat({RTs},length(Files),1);
        end
        for f=1:length(RTs)
            if iscell(RTs{f})
                RTs{f}=horzcat(RTs{f}{:});
            end
        end
    else
        RTs=repmat({RT_Sample},length(Files),1);
    end
    
    OriginalRT=RTs;
    
    % Other
    ShowGraphs=options.showgraphs;
    Window=options.window;
    smoothParam=options.smoothparam;
    
    % Preallocate Arrays
    AUC=cell(length(Files),1);
    Height=cell(length(Files),1);
    PeakTime=cell(length(Files),1);
    TotalTime=cell(length(Files),1);
    Xs=cell(length(Files),1);
    Ys=cell(length(Files),1);
    SmoothYs=cell(length(Files),1);
    XGausses=cell(length(Files),1);
    YGausses=cell(length(Files),1);
    MedianValues=cell(length(Files),1);
    LeftBounds=cell(length(Files),1);
    RightBounds=cell(length(Files),1);
    
    
    for f=1:length(Files)
        RTs=OriginalRT;
        RT=RTs{f};
        
        
        % Import Data
        Data=importdata(Files{f});
        X=Data.data(:,1);
        Y=Data.data(:,2);
        
        % Clip Data
        b=min(RT)-Window;
        e=max(RT)+Window;
        
        [~,beginning]=min(abs(X-b));
        [~,ending]=min(abs(X-e));
        X=X(beginning:ending);
        Y=Y(beginning:ending);
        
        % Find Peaks
        [Peaks,IsShoulder]=findPeaks_RT(X,Y,RT,Window,Files{f},smoothParam);
        
        % Preallocate Arrays
        ind_peaks=Peaks;
        left=zeros(length(ind_peaks),1);
        right=zeros(length(ind_peaks),1);
        As=zeros(length(ind_peaks),1);
        Bs=zeros(length(ind_peaks),1);
        XGauss=zeros(10000,length(ind_peaks));
        YGauss=zeros(10000,length(ind_peaks));
        medians=zeros(length(ind_peaks),2);
        LeftBounds{f}=zeros(length(ind_peaks),1);
        RightBounds{f}=zeros(length(ind_peaks),1);
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        %              Find Peaks                 %
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        
        for j=1:length(ind_peaks)
            try
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
                %        Baseline Drift Correction        %
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
                [~,a]=min(abs(X-(X(ind_peaks(j))-2)));
                [~,b]=min(abs(X-(X(ind_peaks(j))+2)));
                
                sortY=sort(Y(a:b));
                sortY=sortY(ceil(0.33*length(sortY)));
                f1 = fit(X(a:b),Y(a:b),'poly1', 'Exclude', Y(a:b)>sortY);
                Median=f1.p1*X+f1.p2;
                
                
                medians(j,1)=f1.p1;
                medians(j,2)=f1.p2;
                
                
                y=Y-Median;
                
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
                %            Find Boundaries              %
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
                [lefts, rights, XGAUSS, YGAUSS, smoothY, ~, A, B]=findboundaries(X, y, ind_peaks(j),Files{f},IsShoulder(j),smoothParam);
                
                Median2=f1.p1*XGAUSS+f1.p2;
                YGAUSS=YGAUSS+Median2;
                smoothY=smoothY+Median;
                left(j)=lefts;
                right(j)=rights;
                As(j)=A;
                Bs(j)=B;
                XGauss(:,j)=XGAUSS;
                YGauss(:,j)=YGAUSS;
                LeftBounds{f}(j)=A;
                RightBounds{f}(j)=B;
                
            catch err
                save('Error.mat','err')
                disp(['Folder: ' Files{f}])
                disp(err.message);
                
                left(j)=NaN;
                right(j)=NaN;
                XGauss(:,j)=NaN;
                YGauss(:,j)=NaN;
            end
            
        end
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        %              Find Area                  %
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        Xs{f}=X;
        Ys{f}=Y;
        SmoothYs{f}=smoothY;
        AUC{f,1}=zeros(length(ind_peaks),1);
        Height{f,1}=zeros(length(ind_peaks),1);
        TotalTime{f,1}=zeros(length(ind_peaks),1);
        XGausses{f}= XGauss;
        YGausses{f}= YGauss;
        MedianValues{f}=medians;
        
        for j=1:length(Peaks)
            
            try
                AUC1=trapz(XGauss(left(j):right(j),j),YGauss(left(j):right(j),j)); % Find total Area Under Curve
                AUC2=trapz(XGauss(left(j):right(j),j),XGauss(left(j):right(j),j)*medians(j,1)+medians(j,2));
                AUC{f,1}(j)=AUC1-AUC2;
                AUC{f,1}(j)=round(60*AUC{f,1}(j)); %Convert to minutes
                PeakTime{f,1}(j)=roundn(X(Peaks(j)),-3);
                Height{f,1}(j)=roundn(Y(Peaks(j))-min(YGauss(left(j):right(j))),-1);
                TotalTime{f,1}(j)=XGauss(right(j))-XGauss(left(j));
                
            catch
                AUC{f,1}(j)=NaN;
                AUC{f,1}(j)=NaN;
                PeakTime{f,1}(j)=NaN;
                Height{f,1}(j)=NaN;
                TotalTime{f,1}(j)=NaN;
            end
            
        end
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        %                Plot                     %
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
        if ShowGraphs
            
            figure; hold on
            plot(X,Y,'b');
            YLim=ylim;
            
            for j=1:length(Peaks)
                try
                    area(XGauss(left(j):right(j),j),YGauss(left(j):right(j),j),'FaceColor',[205/256 229/255 255/255] ,'EdgeColor','w');% Fill in AUC
                    area(XGauss(left(j):right(j),j),XGauss(left(j):right(j),j)*medians(j,1)+medians(j,2),'FaceColor','w','EdgeColor','none');%Subtract out area under baseline
                end
            end
            plot(X,Y,'b')
            plot(X,smoothY,'k')
            set(gca,'Ylim',YLim);
            set(gca,'FontSize',16)
            title(Files{f}(1:end-4),'Interpreter','None','FontSize',16)
            
            %Override scientific notation for the x-axis
            n=get(gca,'YTick');
            set(gca,'YTickLabel',strread(num2str(n),'%s'));
            set(gca,'Layer','top');
            hold off
        end
    end
    
    % Create Double Array of Data
    N_RTs=0;
    for f=1:length(Files)
        N_RTs=max(N_RTs, length(RTs{f}));
    end
    
    AllRTs=NaN*ones(length(Files),N_RTs);
    PeakAreas=NaN*ones(length(Files),N_RTs);
    PeakTimes=NaN*ones(length(Files),N_RTs);
    
    for f=1:length(Files)
        AllRTs(f,1:length(RTs{f})) = RTs{f};
        PeakAreas(f,1:length(RTs{f}))=AUC{f};
        PeakTimes(f,1:length(RTs{f}))=PeakTime{f};
    end
    
    % Sort Output alphabetically
    [~,sort_index]=sort_nat(Files);
    Files=Files(sort_index);
    PeakAreas=PeakAreas(sort_index,:);
    PeakTimes=PeakTimes(sort_index,:);
    
    [p1,~]=size(Files);
    if p1==1
        Files=Files';
    end
    Output.Files=Files;
    Output.InputTimes=AllRTs;
    Output.PeakTimes=PeakTimes;
    Output.Areas=PeakAreas;
function [FinalPeaks, IsShoulder]=findPeaks_RT(X,Y,RT,Window,File,smoothParam)
    
    % This function finds the peaks based on the retention time.
    if ~exist('Window','var')
        Window=(max(X)-min(X))./5;
    end
    
    [RT,IX] = sort(RT);
    
    % Find peaks in smoothed data. Peaks must be separated by at least 0.1
    % seconds.
    smoothY=perfectSmoother(Y,smoothParam);
    
    [~,peaks]=findpeaks(smoothY,'MinPeakDistance',ceil(0.01./mean(diff(X))));
    
    % Finds all distinct peaks that are within the window of each RT.
    Peaks=[];
    
    for r=1:length(RT)
        
        minWindow=RT(r)-(Window/2);
        maxWindow=RT(r)+(Window/2);
        
        [~,a]=min(abs(X-minWindow));
        [~,b]=min(abs(X-maxWindow));
        
        f1 = fit(X(a:b),Y(a:b),'poly1', 'Exclude', Y(a:b)>median(Y(a:b)));
        Median=f1.p1*X+f1.p2;
        
        tempPeaks=peaks;
        tempPeaks=tempPeaks(tempPeaks>a & tempPeaks<=b);
        tempPeaks=tempPeaks(Y(tempPeaks) > 0.01*(max(Y(a:b) - Median(a:b)))+Median(tempPeaks));
        
        Peaks=[Peaks;tempPeaks];
        
    end
    
    Peaks=Peaks(X(Peaks)>(min(RT)-Window));
    Peaks=Peaks(X(Peaks)<(max(RT)+Window));
    Peaks=unique(Peaks);
    
    % Find Critical Points
    CriticalPoints=findCritPoints(X,smoothY,RT,Window,Peaks,smoothParam);
    
    IsShoulder=[zeros(length(Peaks),1); ones(length(CriticalPoints),1)];
    Peaks=[Peaks; CriticalPoints];
    [~,I]=sort(Peaks);
    Peaks=Peaks(I);
    IsShoulder=IsShoulder(I);
    
    DiffRT=diff(RT);
    
    if triangle(length(Peaks),length(RT))>=triangle(25,9)
        
        t=triangle(length(RT):length(Peaks),length(RT));
        t=t(t<triangle(25,9));
        
        new_n=length(t) + length(RT) - 1;
        
        [~,largest_peaks]=sort(Y(Peaks),'descend');
        Peaks=Peaks(largest_peaks(1:new_n));
        [~,largest_peaks]=sort(Peaks,'ascend');
        Peaks=Peaks(largest_peaks);
        
    end
    
    combinations=combnk(1:length(Peaks),length(RT));
    PeakComb=X(Peaks(combinations));
    if size(PeakComb,2)~=length(RT)
        PeakComb=PeakComb';
    end
    DiffPeakComb=diff(PeakComb,1,2);
    
    % Find the difference between the relative retention times for each
    % combination and the defined relative retention times. Choose the
    % combination that has the minimum sum of the absolute values of these
    % retention times.
    try
        temp=DiffPeakComb-repmat(DiffRT,size(combinations,1),1);
        temp=temp.^2;
        temp=sum(temp,2);
        [~,best]=min(temp);
        
        add_nan=0;
        while min(temp)>length(RT)
            disp(['Cannot identify all peaks for file ' File '. Fitting only available peaks.'])
            temp=DiffPeakComb-repmat(DiffRT,size(combinations,1),1);
            temp=temp.^2;
            
            if sum(sum(temp>1,1) == size(temp,1));
                missing=find( sum(temp>1,1) == size(temp,1),1,'first');
            else
                [~, missing]=max(sum(temp>1,1));
            end
            
            RT(missing+1)=[];
            IX(missing+1)=[];
            IX=[IX,IX(missing+1)];
            add_nan=add_nan+1;
            
            DiffRT=diff(RT);
            
            combinations=combnk(1:length(Peaks),length(RT));
            PeakComb=X(Peaks(combinations));
            if size(PeakComb,2)~=length(RT)
                PeakComb=PeakComb';
            end
            DiffPeakComb=diff(PeakComb,1,2);
            
            temp=DiffPeakComb-repmat(DiffRT,size(combinations,1),1);
            temp=temp.^2;
            temp=sum(temp,2);
            [~,best]=min(temp);
        end
        
    end
    
    combinations=combinations(best,:);
    
    FinalPeaks = Peaks(combinations);
    IsShoulder = IsShoulder(combinations);
    
    if exist('missing','var')
        FinalPeaks=[FinalPeaks;NaN*ones(add_nan,1)];
        IsShoulder=[IsShoulder;NaN*ones(add_nan,1)];
    end
    
    if isempty(combinations)
        warning('Not enough peaks were found. Consider increasing the window size.');
        FinalPeaks=NaN*ones(size(RT));
        IsShoulder=NaN*ones(size(RT));
    else
        [~,RevIX]=sort(IX);
        FinalPeaks=FinalPeaks(RevIX);
        IsShoulder=IsShoulder(RevIX);
    end
function [critical_points]=findCritPoints(X,Y,RT,Window,Peaks,smoothParam)
    smoothY=Y;
    velocity=forward_derivative(X,smoothY);
    velocity=velocity./max(abs(velocity));
    velocity=perfectSmoother(velocity,smoothParam*2);
    
    minWindow=RT(1)-(Window/2);
    maxWindow=RT(end)+(Window/2);
    
    [~,a]=min(abs(X-minWindow));
    [~,b]=min(abs(X-maxWindow));
    
    velocity=velocity./max(abs(velocity(a:b)));
    
    criticalps=find(abs(velocity)<0.01);
    clusters=findclusters(criticalps,1);
    clusters=criticalps(clusters);
    critical_points=zeros(size(clusters,1),1);
    
    % The critical points cluster together. Choose only one from each
    % cluster (the point with the velocity closest to zero).
    SameSign=zeros(size(clusters,1),1);
    
    for i=1:size(clusters,1)
        Sign=sign(velocity(clusters(i,1):clusters(i,2)));
        SameSign(i)=sum(diff(Sign)~=0);
        
        [~,m]=min(abs(velocity(clusters(i,1):clusters(i,2))));
        critical_points(i)=m+clusters(i,1)-1;
    end
    
    % Make sure the velocity of the points of the cluster changes signs no
    % more than once. If it changes signs often, it probably isn't a
    % 'shoulder' in the sense that we are looking for.
    i=1;
    while i<=size(critical_points,1)
        if SameSign(i)>1
            critical_points(i)=[];
            SameSign(i)=[];
            
        else
            i=i+1;
        end
        
    end
    
    % Shoulders must be above the median
    f1 = fit(X,Y,'poly1', 'Exclude', Y>median(Y));
    Median=f1.p1*X+f1.p2;
    critical_points=critical_points(Y(critical_points)>Median(critical_points));
    
    [~,mins]=findpeaks(-smoothY,'MinPeakDistance',ceil(0.1./mean(diff(X))));
    [~,maxs]=findpeaks(smoothY,'MinPeakDistance',ceil(0.1./mean(diff(X))));
    
    % Remove any critical points that are close to a peak. They are
    % describing the same option.
    dt=ceil(0.005./mean(diff(X)));
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
    
    minWindow=RT(1)-(Window/2);
    maxWindow=RT(end)+(Window/2);
    
    [~,a]=min(abs(X-minWindow));
    [~,b]=min(abs(X-maxWindow));
    
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
    
    % Remove any critical points that are close to a peak. They are
    % describing the same option.
    dt=ceil(0.01./mean(diff(X)));
    p=1;
    while p<=length(critical_points)
        DiffCrit=abs(Peaks-critical_points(p));
        if min(DiffCrit)<=dt
            critical_points(p)=[];
        else
            p=p+1;
        end
    end
    
    % Finds all distinct peaks that are withing the window of each RT.
    CriticalPoints=[];
    
    for r=1:length(RT)
        
        minWindow=RT(r)-(Window/2);
        maxWindow=RT(r)+(Window/2);
        
        [~,a]=min(abs(X-minWindow));
        [~,b]=min(abs(X-maxWindow));
        
        f1 = fit(X(a:b),Y(a:b),'poly1', 'Exclude', Y(a:b)>median(Y(a:b)));
        Median=f1.p1*X+f1.p2;
        
        tempCPs=critical_points;
        tempCPs=tempCPs(tempCPs>a & tempCPs<=b);
        tempCPs=tempCPs(Y(tempCPs) > 0.01*(max(Y(a:b) - Median(a:b)))+Median(tempCPs));
        
        CriticalPoints=[CriticalPoints;tempCPs];
        
    end
    CriticalPoints=unique(CriticalPoints);
    critical_points=CriticalPoints;
function [left, right, XGAUSS, YGAUSS, smoothY, gof, A, B]=findboundaries(X, Y, ind_peaks, file,IsShoulder,smoothParam)
    
    % This function finds the boundaries on the left and right side of each
    % peak.
    
    right=zeros(length(ind_peaks),1);
    left=zeros(length(ind_peaks),1);
    
    N=10000;
    XGAUSS=zeros(N,length(ind_peaks));
    YGAUSS=zeros(N,length(ind_peaks));
    
    % Compute Derivatives
    smoothY=perfectSmoother(Y,smoothParam);
    velocity=forward_derivative(X,smoothY);
    velocity=velocity./max(abs(velocity));
    smooth_velo=perfectSmoother(velocity,smoothParam*2);
    
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
            smoothY=perfectSmoother(Y,smoothParam);
            SignVelo=sign(smooth_velo);
            
            % 'A' is the critical point to the left of the peak
            % 'B' is the critical point to the right of the peak
            % 'C' are the local minima in the trace
            % 'D' is a left-hand shoulder
            % 'E' is a right-hand shoulder
            
            if IsShoulder(j)
                A=find(SignVelo(1:ind_peaks(j)-2*dt)<1,1,'last');
                B=find(SignVelo(ind_peaks(j)+ 2*dt :end)>-1,1,'first');
                
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
                
                dt1=ceil(0.1./mean(diff(X)));
                leftHeight=mean(Y(ind_peaks(j)-dt1:ind_peaks(j)));
                rightHeight=mean(Y(ind_peaks(j):ind_peaks(j)+dt1));
                
                if leftHeight>rightHeight
                    A=ind_peaks(j);
                    
                else
                    B=ind_peaks(j);
                end
                
            else
                
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
                
                [~,C]=findpeaks(-smoothY,'threshold',1);
                leftC=C(C<ind_peaks(j));
                rightC=C(C>ind_peaks(j));
                
                % 'A' is the left-hand boundary to be fit.
                % 'B' is the right-hand boundary to be fit.
                A=max([leftC; A]);
                B=min([rightC; B]);
                
            end
            
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
                Upper=[1.05*Y(ind_peaks(j)) max(X(A:B)) 3];
                StartPoint = [max(Y(A:B)) median(X(A:B)) 1];
                
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
    elseif leftvelo<0.25 || rightvelo<0.25
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
function z = perfectSmoother(y,lambda)
    
    m = length(y);
    
    if m>1000
        
        E = speye(m);
        D = diff(E);
        C = chol(E + lambda * (D' * D));
        z = C\(C'\y);
    else
        
        E = eye(m);
        D = diff(E);
        z = (E + lambda * (D' * D))\y;
        
    end
function t=triangle(n,k)
    % This function returns Pascal's triangle, where 'n' is the row or rows
    % and 'k' is the kth entry or entries in that row.
    
    % Ex: t=triangle(5,0:5);
    % t = [1,5,10,10,5,1]
    
    % Ex: t=triangle(0:4,0:4)
    % t=[1,  NaN,  NaN,  NaN,  NaN;
    %    1,  1,    NaN,  NaN,  NaN;
    %    1,  2,    1,    NaN,  NaN;
    %    1,  3,    3,    1,    NaN;
    %    1,  4,    6,    4,    1];
    
    % fact(n)./(fact(k)*fact(n-k))
    
    t=zeros(length(n),length(k));
    
    for N=1:length(n)
        for K=1:length(k)
            if (n(N)-k(K))>=0
                t(N,K)=prod(k(K)+1:n(N))./factorial(n(N)-k(K));
            else
                t(N,K)=NaN;
            end
        end
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
    [~,index] = sortrows(comp);
    if is_descend
        index = index(end:-1:1);
    end
    index = reshape(index,size(c));
    cs = c(index);
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