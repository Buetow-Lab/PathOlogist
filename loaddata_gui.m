function varargout = loaddata_gui(varargin)
% LOADDATA_GUI M-file for loaddata_gui.fig
%      LOADDATA_GUI, by itself, creates a new LOADDATA_GUI or raises the existing
%      singleton*.
%
%      H = LOADDATA_GUI returns the handle to a new LOADDATA_GUI or the handle to
%      the existing singleton*.
%
%      LOADDATA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOADDATA_GUI.M with the given input arguments.
%
%      LOADDATA_GUI('Property','Value',...) creates a new LOADDATA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before loaddata_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to loaddata_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help loaddata_gui

% Last Modified by GUIDE v2.5 26-Jun-2010 16:59:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @loaddata_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @loaddata_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function loaddata_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to loaddata_gui (see VARARGIN)

% import data
handles.output = hObject;
u = get(hObject, 'UserData');
[handles.DataStruct handles.SampleData handles.ChipData handles.ButtonEnable] = deal(u{:});

%set up vars
handles.cleared = 0;
handles.next = 0;
handles.Buttons = [handles.rma_calcButton, handles.find_rmaButton, handles.udp_calcButton, ...
    handles.find_udpButton, handles.sample_labelsButton, handles.find_hitsButton, ...
    handles.chipMenu, handles.find_actButton, handles.find_conButton];

%update screen
for ii = [2 4 6 8 9]
    update_loaddata(handles, handles.Buttons(ii));
end

%enable/disable buttons
for ii = 1:length(handles.ButtonEnable)
   set(handles.Buttons(ii),'Enable',handles.ButtonEnable{ii});
end
set(handles.chipMenu,'String', handles.ChipData.Names, 'Value', handles.ChipData.Selection);

% Update handles structure
guidata(hObject, handles);
uiwait

function varargout = loaddata_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.DataStruct;
varargout{3} = handles.SampleData;
varargout{4} = handles.ChipData;
varargout{5} = handles.ButtonEnable;
varargout{6} = handles.cleared;
varargout{7} = handles.next;

function chipMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chipMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_loaddata(handles, button)
%determine data type
Indices = get(button, 'UserData');
DataType = handles.DataStruct(Indices(1)).Data(Indices(2)).Type;
eval(['TextHandle = handles.' DataType 'Text;']);

%if successfully loaded, update button
    if handles.DataStruct(Indices(1)).Data(Indices(2)).Loaded
        set(button,'String', [DataType, ' Loaded']);
        set(button, 'BackgroundColor', 'yellow');
        set(TextHandle, 'String', handles.DataStruct(Indices(1)).Data(Indices(2)).Path);
    end
    
    button_enable(button,handles);

function find_fileButton_Callback(hObject, eventdata, handles)
% hObject    handle to find_udpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file_name, file_dir] = uigetfile('*.txt', 'Pick a text file');
if file_name
    loadData(hObject, file_name, file_dir, handles);
end

function loadData(button, file_name, file_dir, handles)
% hObject    handle to load_dataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

BLabel = get(button, 'String');
set(button, 'String', 'Loading...');
set(handles.RMA_progressText, 'String','');
drawnow expose
error_strings = {'Error: Column labels do not match previously loaded data'; 'Error: Row labels do not match previously loaded data'; 'Unknown Error'};

try
    Indices = get(button, 'UserData');

    Directory = cd;
    cd(file_dir);

    err = 0;

    newData1 = importdata(file_name, '\t', 1); % not use importdata
    % function due to NaN (2019-05-12)

%     % use textscan function to read all data w/ NaN (2019-05-12)
%     fid=fopen(file_name);
%     text=textscan(fid, '%s','Delimiter','\n');
%     fclose(fid);
% 
%     tmp=strsplit(text{1}{1},'\t');
%     s = size(tmp,2);
%     textdata(1,2:(s+1))=tmp;
% 
%     for i = 2:cellfun('length',text)
%         tmp=strsplit(text{1}{i},'\t');
%         textdata(i,1)=tmp(1);
%         s = size(tmp,2);
%         raw_data((i-1),:)=cellfun(@str2double,tmp(2:s));
%     end;
% 
%     newData1.data=raw_data; % assign to newData1 as one of cell array (2019-05-12)
%     newData1.textdata=textdata; % assign to newData1 as one of cell array (2019-05-12)

    if Indices(2)~=3
        data = newData1.data;
        handles.DataStruct(Indices(1)).Data(Indices(2)).Data=data;
    else
        data = newData1.data(:,2:end);
        handles.DataStruct(Indices(1)).Data(Indices(2)).Data{2}=data;
        newData1.textdata(2:size(data,1)+1,1) = num2cell(newData1.data(:,1));
    end

    [r c] = size(newData1.textdata);

    if r == size(data,1) + 1
        probe_names = newData1.textdata(2:end,1);
    else
        err = 1;
    end
    
    if c == size(data,2) + 1
        sample_names = newData1.textdata(1,2:end);
    else
        sample_names = textread(file_name,'%s',size(data,2)+1, 'delimiter', '\t');
        sample_names = sample_names(2:end)';
    end
    
    cd(Directory);

    Rows = handles.DataStruct(Indices(1)).Rows;
    Cols = handles.DataStruct(Indices(1)).Cols;
    if ~err
        if isempty(Cols)
            handles.DataStruct(Indices(1)).Cols = sample_names;
        else
            if length(Cols) ~= length(sample_names) || sum(strcmp(sample_names, Cols))~=length(sample_names)
                err = 1;%error: Col labels dont match
            end
        end

        if Indices(2)~=3
            if isempty(Rows) || sum(Indices==[1 2])==2
                handles.DataStruct(Indices(1)).Rows = probe_names;
            else
                if length(Rows) ~= length(probe_names) || sum(strcmp(probe_names, Rows))~=length(probe_names)
                    err = 2; %error: Row labels dont match
                end
            end
        else
            handles.DataStruct(Indices(1)).Data(Indices(2)).Data{1}=probe_names;
        end

        handles.DataStruct(Indices(1)).Data(Indices(2)).Loaded= true;
        handles.DataStruct(Indices(1)).Data(Indices(2)).Path = [file_dir file_name];

        
        if Indices(1)==2
            handles.SampleData.IDs = sample_names;
        end
        
        guidata(handles.figure1, handles);
        update_loaddata(handles, button);

    else
        set(handles.RMA_progressText, 'String', error_strings(err));
    end

catch ME
    set(button, 'String', BLabel);
    set(handles.RMA_progressText, 'String', 'Error: Make sure you have selected a valid text file.');
end

function clear_data1Button_Callback(hObject, eventdata, handles)
% hObject    handle to clear_data1Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
d1 = struct('Type', {'RMA', 'UDP', 'Hits'}, 'Data', {[] [] cell(1,2)}, 'Path', [], 'Loaded', false);
fields = {'Rows' 'Cols' 'RowSelect', 'ColSelect'};
handles.DataStruct(1).Data = d1;
for ff = 1:length(fields)
    handles.DataStruct = setfield(handles.DataStruct, {1}, fields{ff}, []);
end
handles.cleared=1;
handles.ButtonEnable(1:7) = {'on' 'on' 'off' 'on' 'off' 'on' 'off'};
guidata(handles.figure1, handles);
uiresume

%{
for ii = [2 4 6]
    button = handles.Buttons(ii);
    Indices = get(button, 'UserData');
    DataType = handles.DataStruct(Indices(1)).Data(Indices(2)).Type;
    set(button, 'String', ['Find ' DataType ' File'], 'BackgroundColor', 'white', 'Enable', 'on');
    eval(['TextHandle = handles.' DataType 'Text;']);
    set(TextHandle, 'String', ' ');
end
%}


function clear_data2Button_Callback(hObject, eventdata, handles)
% hObject    handle to clear_data2Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%{
for ii = 4:5
    button = handles.buttons(ii);
    Indices = get(button, 'UserData');
    DataType = handles.DataStruct(Indices(1)).Data(Indices(2)).Type;
    set(button, 'String', ['Find ' DataType ' File'], 'BackgroundColor', 'white', 'Enable', 'on');
    eval(['TextHandle = handles.' DataType 'Text;']);
    set(TextHandle, 'String', ' ');
end

set(handles.load_dataInstr,'Visible','on');
set(handles.load_dataError, 'Visible','off');
%}

d2 = struct('Type', {'Activity', 'Consistency', 'Blank'}, 'Data', {[] [] cell(1,2)}, 'Path', [], 'Loaded', false);
fields = {'Rows' 'Cols' 'RowSelect', 'ColSelect'};
handles.DataStruct(2).Data = d2;
for ff = 1:length(fields)
    handles.DataStruct = setfield(handles.DataStruct, {2}, fields{ff}, []);
end
handles.cleared=1;
handles.ButtonEnable(8:9) = {'on' 'on'};
guidata(handles.figure1, handles);
uiresume

function chipMenu_Callback(hObject, eventdata, handles)
% hObject    handle to chipMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns chipMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from chipMenu
ChipList = get(handles.chipMenu, 'String');
ChipPick = get(handles.chipMenu, 'Value');
if ChipPick == length(ChipList)
    NewChipName = inputdlg('Enter a name for the new platform');
    [filename filedir] = uigetfile('*.txt', 'Select a probe-geneID linking file for the new platform');
    NewChipList = [ChipList(1:end-1);NewChipName; ChipList(end)];
    NewChipFiles = [handles.ChipData.Files;{[filedir '/' filename]}]; 
    set(handles.chipMenu, 'String', NewChipList);
    handles.ChipData.Names = NewChipList;
    handles.ChipData.Files = NewChipFiles;
    guidata(handles.figure1, handles);
end

function backButton_Callback(hObject, eventdata, handles)
% hObject    handle to backButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for ii = 1:length(handles.Buttons)
    handles.ButtonEnable{ii} = get(handles.Buttons(ii),'Enable');
end
guidata(handles.figure1, handles);
if handles.DataStruct(1).Data(2).Loaded || handles.DataStruct(1).Data(1).Loaded
    chip = get(handles.chipMenu, 'Value');
    handles.ChipData.Selection=chip;
    if chip > 1
        fid = fopen(handles.ChipData.Files{chip});
        pl = textscan(fid, '%s%s', 'delimiter', '\t');
        fclose(fid);
        handles.ChipData.Links = horzcat(pl{:});
        guidata(handles.figure1, handles);
        uiresume
    else
        if handles.DataStruct(1).Data(2).Loaded
        set(handles.choose_chipText, 'Visible', 'on')
        else
          uiresume  
        end
    end
else
    uiresume
end


% --- Executes on button press in rma_calcButton.
function rma_calcButton_Callback(hObject, eventdata, handles)
% hObject    handle to rma_calcButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    set(handles.RMA_progressText, 'String', 'Reading Files...');
    probeData = celintensityread(' ', ' ','pmonly', true);
    chipName = probeData.CDFName(1:strfind(probeData.CDFName,'.')-1);
    if strcmp(chipName(1:3), 'HT_')
        chipName = chipName(4:end);
    end
    chipNum = strmatch(chipName, handles.ChipData.Names, 'exact');
    if ~isempty(chipNum)
    set(handles.chipMenu, 'Value', chipNum, 'Enable', 'on');
    end
    set(handles.RMA_progressText, 'String', 'Performing Background Adjustment...');
    drawnow expose
    pms_bg = rmabackadj(probeData.PMIntensities);

set(handles.RMA_progressText, 'String', 'Performing Quantile Normalization...');
    drawnow expose
    pms_bgnorm = quantilenorm(pms_bg);
   
    set(handles.RMA_progressText, 'String', 'Summarizing Results...');
    %drawnow expose
    rma_exp = rmasummary(probeData.ProbeIndices, pms_bgnorm, 'Output', 'log2');
    
   
    
    [filename pathname] = uiputfile('*.txt', 'Save RMA as:');
    
    RMA_data = [[{chipName} probeData.CELNames'];[probeData.ProbeSetIDs num2cell(rma_exp)]];
    fid = fopen([pathname '/' filename], 'wt');
    fprintf(fid, '%s\t', RMA_data{1,1:end-1});
    fprintf(fid, '%s', RMA_data{1,end});
    for ii = 2:size(RMA_data,1)
       fprintf(fid, '\n%s', RMA_data{ii,1});
       fprintf(fid, '\t%10.6f', RMA_data{ii,2:end});
    end
    fclose(fid)
    
    set(handles.RMA_progressText, 'String', '');
    set(handles.RMAText, 'String',[pathname filename]);
    drawnow expose
    
    handles.DataStruct(1).Data(1).Loaded = true;
    handles.DataStruct(1).Rows= probeData.ProbeSetIDs;
    handles.DataStruct(1).Data(1).Data = rma_exp;
    handles.DataStruct(1).Data(1).Path = [pathname filename];
    handles.DataStruct(1).Cols = probeData.CELNames';
    guidata(hObject, handles);
    update_loaddata(handles, handles.find_rmaButton);
catch ME
    Error = ME.identifier;
    switch Error
        case 'Bioinfo:celintensityread:RowsOrColsNotMatch'
            set(handles.RMA_progressText, 'String', 'Error: Library (.cdf) file does not match .cel files.');
        otherwise
            set(handles.RMA_progressText, 'String', 'Unknown Error.');
    end
end

% --- Executes on button press in udp_calcButton.
function udp_calcButton_Callback(hObject, eventdata, handles)
% hObject    handle to udp_calcButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% try
    RMA = handles.DataStruct(1).Data(1).Data;
    [rows cols] = size(RMA);
    invalid = zeros(1,rows);
    %{
        for i = 1:rows
            for j = 1:cols
                if RMA(i,j)<0
                invalid(i) = 1;
            end
        end
    %}
    invalid = logical(invalid);
    sample_names = handles.DataStruct(1).Cols;
    probe_names = handles.DataStruct(1).Rows;
    UDP = UDP_calc(RMA);
 
    if ~isempty(UDP)
        UDP(invalid,:) = NaN;
        UDP_data = [[{'-'} sample_names]; [probe_names num2cell(UDP)]];
        handles.DataStruct(1).Data(2).Loaded = true;
        handles.DataStruct(1).Data(2).Data = UDP;

        Directory = cd;
        [FileName FileDir] = uiputfile('*.txt', 'Save UDP as: ');
        handles.DataStruct(1).Data(2).Path = [FileName '\' FileDir];

        cd(FileDir);
        fid = fopen(FileName, 'wt');
        fprintf(fid, '%s\t', UDP_data{1,1:end-1});
        fprintf(fid, '%s', UDP_data{1,end});
        for rr = 2:rows+1
            fprintf(fid,'\n%s', UDP_data{rr,1});
            fprintf(fid, '\t%1.6f', UDP_data{rr,2:end});
        end
        fclose(fid);
        cd(Directory);

        set(handles.RMA_progressText, 'String', '');
        set(handles.UDPText, 'String',[FileDir FileName]);
        drawnow expose
    
        guidata(hObject, handles);  
        update_loaddata(handles, handles.find_udpButton);
    end
% catch ME
%     w=waitbar(.5);
%     delete(w)
%     LastStack = ME.stack(1).name;
%     switch LastStack
%         case 'fzero'
%             set(handles.RMA_progressText, 'String', 'Error: Not enough samples to perform normalization.');
%         otherwise
%             set(handles.RMA_progressText, 'String', 'Unknown Error.');
%     end
% end



function UDP_matrix = UDP_calc(RMA)
%takes a MxN matrix 'RMA' of rma-adjusted expression values (with M probes, N samples for each) and creates an MxN matrix of 
%probabilities that each value is in the 'up' state
RMA = double(RMA);
min_exprs = min(min(RMA));
[num_probes num_samples] = size(RMA);
UDP_matrix = zeros(size(RMA));
W = waitbar(0, 'Calculating... 0% complete', 'CreateCancelBtn', 'set(gcf, ''Clipping'', ''off'')', 'Clipping', 'on');
PctComplete=0;
for i = 1:num_probes
    if strcmp(get(W, 'Clipping'), 'on')
        if round(i*100/num_probes)>PctComplete
            PctComplete= round(i*100/num_probes);
            waitbar(i/num_probes, W, ['Calculating... ' num2str(PctComplete) '% complete']);
        end
        if all(RMA(i,:)==min_exprs)~=1 % smallest exp value indicate the raw expr are zero (2019-05-06)
            if min(RMA(i,:))<0
                rma = RMA(i,:) - min(RMA(i,:)) + 0.1; % 2019_05_08 add this syntax to avoid any data point less than 0 within the data
            else
                rma = RMA(i,:);
            end
            [AB E M adj] = mixgam(rma);
            % prob = zeros(1,num_samples);
            if M == 1 %for a mixed distr 
                % for j = 1:num_samples
                    prob = getP_Present_Absent_x(AB(1,1), AB(1,2), AB(2,1), AB(2,2), E, adj, rma);
                % end

            else %if unmixed, decide whether all up or down based on gaussian alpha parameter
                if normfit(rma) > 1
                    prob(1:num_samples) = 1;
                else
                    prob(1:num_samples) = 0;
                end
            end
        else
            prob(1:num_samples) = 0;
        end
        UDP_matrix(i,:) = prob(:);
        drawnow
    else
        UDP_matrix = [];
        delete(W)
        break
    end 
    
end
try
        delete(W)
    end

function [alphabeta eta mixture adj]= mixgam(data)

%data should be a row vector
if (size(data,1)>size(data,2))
    data = data';
end
%data values smaller than eps should be moved to be eps
% data(data<eps)=eps;
C=2;
loop = 0;
while 1
    % try
        idx = kmeans(data',C);
%         data1 = data(find(idx==1));
%         data2 = data(find(idx==2));
% 
%         phat1 = gamfit(data1);
%         phat2 = gamfit(data2);
% 
%         alphabeta(1,1) = phat1(1);
%         alphabeta(1,2) = phat1(2);
%         alphabeta(2,1) = phat2(1);
%         alphabeta(2,2) = phat2(2);
%     catch
%         alphabeta(1,1) = 0.9;
%         alphabeta(1,2) = 200;
%         alphabeta(2,1) = 3;
%         alphabeta(2,2) = 400;
%     end
% 
%     aa = rand(1);
%     eta = [aa 1-aa];
%   generate the initial gamma coefficients and weight (2019-05-08)
    for conta=1:C
        eta(conta) = sum(idx==conta)/length(data);        
        alphabeta(conta,2) = var(data(idx==conta))/mean(data(idx==conta));
        alphabeta(conta,1) = mean(data(idx==conta))/alphabeta(conta,2);
    end
    iter = 0;
    while ( ( ( iter < 3 ) || ( abs((loglike(end) - loglike(end-1))/loglike(end)) > 1e-5 ) ) && (iter < 1000) )
        for ii=1:C
            w(ii,:) = eta(ii) * gampdf(data,alphabeta(ii,1),alphabeta(ii,2));
        end
        w = w ./ (ones(2,1)*sum(w,1));

        eta = sum(w,2)/sum(sum(w));
        A = sum(w,2);
        B = w*data';
        D = w*(log(data)');

        for ii=1:C
            f = @(x) ( (log(x) + log(A(ii)) - log(B(ii)) ) * A(ii) + D(ii) - A(ii)*psi(x) );
            success=0;
            maxrange=1e3;
            minrange=10*eps;
            while ( (success==0) && (maxrange < 1e20) &&  (alphabeta(ii,1)<1e6))
                try
                    alphabeta(ii,1) = fzero(@(x) f(x), [minrange maxrange]);
                    success=1;
                catch
                    maxrange = maxrange*10;
                    minrange = minrange*1000;
    %                disp(sprintf('Increasing maxrange to %d and minrange to %d\n', maxrange, minrange));
                end
            end
            alphabeta(ii,2) = 1/(alphabeta(ii,1) * A(ii)/ B(ii));
        end
        iter = iter + 1;
        loglike(iter) = sum(log(eta(1)*gampdf(data,alphabeta(1,1), alphabeta(1,2)) + eta(2)*gampdf(data, alphabeta(2,1), alphabeta(2,2))));
    end
    if alphabeta(1,1)*alphabeta(1,2) < alphabeta(2,1)*alphabeta(2,2)
            alphabeta=flip(alphabeta);
            eta = flip(eta);
    end
    % for the single case
    onealphabeta = gamfit(data);
    oneloglike = sum(log(gampdf(data, onealphabeta(1), onealphabeta(2))));

    % Make a decision if this is one or two
    dm = C*2 +1;
    % BIC = loglike(end) - log(length(data))*dm/2;
    AIC = loglike(end) - dm;
    if AIC > oneloglike % which means that this is a mixture
        mixture=1;
    else
        mixture=0;
        alphabeta = onealphabeta;
        eta = 1;
    end
    
    adj = 0;
    P1 = 0;
    P2 = 0;
    if sum([sum(isinf(alphabeta)) sum(isnan(alphabeta))])<1 % avoid non-numeric parameters 2019-05-06
        if mixture == 0
            break
        elseif ( sum(idx==1)==1 || sum(idx==2)==1 )
            break
        else % evaluate the precision of the estimated parameters of the gamma mixture model (2019-05-06)
            % exprs = sort(data);
            [freq,exprs] = ecdf(sort(data));
            exprs = exprs';
            pxp = gampdf(exprs, alphabeta(1,1), alphabeta(1,2));
            % pd1=makedist('Gamma','a',alphabeta(2,1),'b',alphabeta(2,2));
            pd1 = cdf(makedist('Gamma','a',alphabeta(1,1),'b',alphabeta(1,2)),exprs);
            pp = eta(1);
            pa = 1-pp;
            pxa = gampdf(exprs, alphabeta(2,1), alphabeta(2,2));
            % pd2=makedist('Gamma','a',alphabeta(1,1),'b',alphabeta(1,2));
            pd2 = cdf(makedist('Gamma','a',alphabeta(2,1),'b',alphabeta(2,2)),exprs);
            P = pp*pxp + pa*pxa;
            PD = pp*pd1 + pa*pd2;
            ppx = (pxp.*pp)./P;
            if ( find(ppx==max(ppx),1)~=1 && (length(find(single(ppx)==min(single(ppx))))+length(find(single(ppx)==max(single(ppx)))))~=length(ppx) ) % avoid single data point to fit in the generalized regression model for logistic regression (2019-05-08)
                % goodness-of-fit test for gamma mixture model (2019-05-06)
                % [freq,observed] = ecdf(exprs);
                bin = 1:size(freq,1);
                [~,~,st] = chi2gof(bin,'Ctrs',bin,'Frequency',freq','Expected',PD,'NParams',1);
                P1 = 1-chi2cdf(st.chi2stat,1);
                % goodness-of-fit test for probability of the x presenting in logitic regression model (2019-05-08)
                [Coef,~] = glmfit(exprs(1:find(ppx==max(ppx)))',ppx(1:find(ppx==max(ppx)))','binomial','link','logit');
                warning off;
                Fit=glmval(Coef,exprs,'logit');
                chi2stat = sum((ppx'-Fit).^2 ./ Fit);
                P2 = 1-chi2cdf(chi2stat,1);
            end
        end
    end
    
    if ((P1 > 0.01) && (P2 > 0.05)) % fit both gamma mixture model and ligistic regression (2019-05-08)
        break
	elseif ((P1 > 0.01) && (loop > 5)) % fit gamma mixture model only (2019-05-08)
        adj = 1;
        break
    elseif loop > 10
        break
    end
    loop=loop+1;
end

function p_present_x = getP_Present_Absent_x(alphaP, betaP, alphaA, betaA, eta, adj, rma)

for i = 1:size(rma,2)
    x = rma(i);
    P_x_present(i) = gampdf(x, alphaP, betaP);
    P_present(i) = eta(1);
    P_absent(i) = 1-P_present(i);
    P_x_absent(i) = gampdf(x, alphaA, betaA);
    P_x(i) = P_present(i)*P_x_present(i) + P_absent(i)*P_x_absent(i);
    p_present_x(i) = P_x_present(i)*P_present(i)/P_x(i);
    if isnan(p_present_x(i))
        p_present_x(i) = 0.5;
    end
end
% adjusting the up/down score curve (2019-05-08)
tmp = sortrows([rma' p_present_x' linspace(1,size(rma,2),size(rma,2))']);
maxp = find(tmp(:,2)==max(tmp(:,2)));
if adj == 1
	[logitCoef,~] = glmfit(tmp(1:maxp,1),tmp(1:maxp,2),'binomial','logit');
	logitFit=glmval(logitCoef,tmp(:,1),'logit');
	warning off;
    tmp(:,2) = logitFit;
    tmp = sortrows(tmp,3);
    p_present_x = tmp(:,2)';
elseif maxp<length(tmp(:,2))
    m = max(tmp(:,2));
	for ii = 1:maxp
        tmp(ii,2)=m;
    end
    tmp = sortrows(tmp,3);
    p_present_x = tmp(:,2)';
end

function handles = button_enable(ButtonHandle, handles)

try
    ButtonNum = find(handles.Buttons==ButtonHandle);
catch
    ButtonNum = 0;
end

switch ButtonNum         
    case {1,2}
        OnOff = [0 0 1 1 1 1 0 1 1];
    case {3,4}
        OnOff = [0 0 0 0 1 1 1 1 1];
    case 6
        OnOff = [2 2 2 2 2 0 2 2 2];
    case 8
        OnOff = [2 2 2 2 2 2 2 0 2];
    case 9
        OnOff = [2 2 2 2 2 2 2 2 0];
    otherwise
        OnOff = [1 1 0 1 0 1 0 1 1]; 
end

for ii = 1:length(OnOff)
    if OnOff(ii)==0
        set(handles.Buttons(ii),'Enable','off');
    elseif OnOff(ii) ==1 
        set(handles.Buttons(ii), 'Enable', 'on');
    end
end



% --- Executes on button press in sample_labelsButton.
function sample_labelsButton_Callback(hObject, eventdata, handles)
% hObject    handle to sample_labelsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Samples = handles.DataStruct(1).Cols;
data = [];
[H data] = entersampledata_gui('UserData', {Samples,data});
close(H)
if ~isempty(data)
    handles.DataStruct(1).Cols = data';
end
guidata(hObject, handles);
