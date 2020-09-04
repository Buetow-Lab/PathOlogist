function varargout = UDP_gui(varargin)
%UDP_GUI M-file for UDP_gui.fig
%      UDP_GUI, by itself, creates a new UDP_GUI or raises the existing
%      singleton*.
%
%      H = UDP_GUI returns the handle to a new UDP_GUI or the handle to
%      the existing singleton*.
%
%      UDP_GUI('Property','Value',...) creates a new UDP_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to UDP_gui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      UDP_GUI('CALLBACK') and UDP_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in UDP_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UDP_gui

% Last Modified by GUIDE v2.5 15-Aug-2008 16:51:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UDP_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @UDP_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
end
% End initialization code - DO NOT EDIT

function UDP_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for UDP_gui
handles.output = hObject;
u = get(hObject, 'UserData');
[handles.DataStruct handles.ChipData] = deal(u{:}); 
handles.next = 0;
% Update handles structure
guidata(hObject, handles);
if handles.DataStruct(1).Data(1).Loaded
    set(handles.calculateButton, 'Enable', 'on');
end    
if handles.DataStruct(1).Data(1).Loaded || handles.DataStruct(1).Data(2).Loaded
    set(handles.chipMenu, 'Enable', 'on', 'String', handles.ChipData.Names, 'Value', handles.ChipData.Selection);
end

 uiwait
end

function varargout = UDP_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.DataStruct;
varargout{3} = handles.ChipData;
varargout{4} = handles.next;
end

function chipMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chipMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sample_labelButton_Callback(hObject, eventdata, handles)
% hObject    handle to sample_labelButton (see GCBO)
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
end

function calculateButton_Callback(hObject, eventdata, handles)
% hObject    handle to calculateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
type = get(handles.udpRadio, 'Value');

RMA = handles.DataStruct(1).Data(1).Data;
[rows cols] = size(RMA);
invalid = zeros(1,rows);
for i = 1:rows
    for j = 1:cols
        if RMA(i,j)<0
            invalid(i) = 1;
        end
    end
end
invalid = logical(invalid);
sample_names = handles.DataStruct(1).Cols;
probe_names = handles.DataStruct(1).Rows;
if type 
    UDP = UDP_calc(RMA);
else
    UDP = PUD_calc(RMA);
end
if ~isempty(UDP)
UDP(invalid,:) = NaN;
UDP_data = [[{'-'} sample_names]; [probe_names num2cell(UDP)]];
handles.DataStruct(1).Data(2).Loaded = true;
handles.DataStruct(1).Data(2).Data = UDP;

guidata(handles.figure1, handles);
Directory = cd;
[FileName FileDir] = uiputfile('*.txt', 'Save UDP as: ');
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
handles.DataStruct(1).Data(2).Path = [FileName '\' FileDir];
set(hObject, 'String', 'Calculated');
set(handles.calculateText, 'Visible', 'on');
set(handles.returnButton, 'Visible', 'on');
guidata(handles.figure1, handles);
end
end
    function PUD_matrix = PUD_calc(RMA)
[num_probes num_samples] = size(RMA);
PUD_matrix = zeros(size(RMA));
for ii = 1:num_probes
    r = RMA(ii,:);
    PUD_matrix(ii,:) = (r-min(r))/(max(r)-min(r));
end
    end
    function UDP_matrix = UDP_calc(RMA)
%takes a MxN matrix 'RMA' of rma-adjusted expression values (with M probes, N samples for each) and creates an MxN matrix of 
%probabilities that each value is in the 'up' state

[num_probes num_samples] = size(RMA);
UDP_matrix = zeros(size(RMA));
W = waitbar(0, 'Calculating... 0% complete', 'CreateCancelBtn', 'set(gcf, ''UserData'', false)', 'UserData', true);
PctComplete=0;
for i = 1:num_probes
    if get(W, 'UserData')
        if round(i*100/num_probes)>PctComplete
            PctComplete= round(i*100/num_probes);
            waitbar(i/num_probes, W, ['Calculating... ' num2str(PctComplete) '% complete']);
            drawnow
        end
    [AB E M] = mixgam(RMA(i,:));
    prob = zeros(1,num_samples);
   if M == 1 %for a mixed distr 
       for j = 1:num_samples
           prob(j) = getP_Present_Absent_x(AB(1,1), AB(1,2), M, AB(2,1), AB(2,2), E, RMA(i,j));
       end
   else %if unmixed, decide whether all up or down based on gaussian alpha parameter
            if AB(1,1)>1
                prob(1:num_samples) = 1;
            else
                prob(1:num_samples) = 0;
            end
   end
    
    UDP_matrix(i,:) = prob(:);
    else
        i = num_probes;
        break
    end 
end
aborted = ~get(W,'UserData');
delete(W)
if  aborted
   UDP_matrix = [];
end
function [alphabeta eta mixture]= mixgam(data)

%data should be a row vector
if (size(data,1)>size(data,2))
    data = data';
end
%data values smaller than eps should be moved to be eps
data(data<eps)=eps;

C=2;


try
    rng(1) % for reproductivity
    idx = kmeans(data',2);
    data1 = data(find(idx==1));
    data2 = data(find(idx==2));

    phat1 = gamfit(data1);
    phat2 = gamfit(data2);

    alphabeta(1,1) = phat1(1);
    alphabeta(1,2) = phat1(2);
    alphabeta(2,1) = phat2(1);
    alphabeta(2,2) = phat2(2);
catch
    alphabeta(1,1) = 0.9;
    alphabeta(1,2) = 200;
    alphabeta(2,1) = 3;
    alphabeta(2,2) = 400;
end


aa = rand(1);
eta = [aa 1-aa];
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

% for the single case

onealphabeta = gamfit(data);
oneloglike = sum(log(gampdf(data, onealphabeta(1), onealphabeta(2))));


% Make a decision if this is one or two
dm = C*2 +1;
BIC = loglike(end) - log(length(data))*dm/2;
AIC = loglike(end) - dm;
if AIC > oneloglike % which means that this is a mixture
    mixture=1;
    avrg1 = alphabeta(1,1)*alphabeta(1,2);
    avrg2 = alphabeta(2,1)*alphabeta(2,2);
    if avrg2>avrg1
        newalphabeta(1,1)=alphabeta(2,1);
        newalphabeta(2,1)=alphabeta(1,1);
        newalphabeta(1,2)=alphabeta(2,2);
        newalphabeta(2,2)=alphabeta(1,2);
        alphabeta=newalphabeta;
    end
else
    mixture=0;
    alphabeta = onealphabeta;
    eta = 1;
end
end
function p_present_x = getP_Present_Absent_x(alphaP, betaP, mixture, alphaA, betaA, eta, x)

P_x_present = gampdf(x, alphaP, betaP);
P_present = eta(1);
P_absent = 1-P_present;
P_x_absent = gampdf(x, alphaA, betaA);
P_x = P_present*P_x_present + P_absent*P_x_absent;
p_present_x = P_x_present*P_present/P_x;
if isnan(p_present_x)
    p_present_x = 0.5;
end
end
    end

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
    NewChipFiles = [handles.ChipData.Files;{[filedir '\' filename]}]; 
    set(handles.chipMenu, 'String', NewChipList);
    handles.ChipData.Names = NewChipList;
    handles.ChipData.Files = NewChipFiles;
    guidata(handles.figure1, handles);
end
end

function nextButton_Callback(hObject, eventdata, handles)
% hObject    handle to nextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.next = 1;
guidata(handles.figure1, handles);
returnButton_Callback(hObject, eventdata, handles)
end

function returnButton_Callback(hObject, eventdata, handles)
% hObject    handle to returnButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.DataStruct(1).Data(2).Loaded
    chip = get(handles.chipMenu, 'Value');
    handles.ChipData.Selection = chip;
    if chip > 1
        fid = fopen(handles.ChipData.Files{chip});
        pl = textscan(fid, '%s%s', 'delimiter', '\t');
        fclose(fid);
        handles.ChipData.Links = horzcat(pl{:});
        guidata(handles.figure1, handles);
        uiresume
    else
        set(handles.chipError, 'Visible', 'on')
    end
else
    uiresume
end

end