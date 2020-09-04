function varargout = PathOlogist(varargin)
% PATHOLOGIST M-file for PathOlogist.fig
%      PATHOLOGIST, by itself, creates a new PATHOLOGIST or raises the existing
%      singleton*.
%
%      H = PATHOLOGIST returns the handle to a new PATHOLOGIST or the handle to
%      the existing singleton*.
%
%      PATHOLOGIST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PATHOLOGIST.M with the given input arguments.
%
%      PATHOLOGIST('Property','Value',...) creates a new PATHOLOGIST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PathOlogist_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PathOlogist_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PathOlogist

% Last Modified by GUIDE v2.5 26-Jun-2010 16:54:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PathOlogist_OpeningFcn, ...
                   'gui_OutputFcn',  @PathOlogist_OutputFcn, ...
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


% --- Executes just before PathOlogist is made visible.
function PathOlogist_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PathOlogist (see VARARGIN)
%additions to handles structure: 
%DataStruct  
%SampleData 
%PwData

%%set up variables to hold data
handles.output = hObject;
handles.version = 1;
d1 = struct('Type', {'RMA', 'UDP', 'Hits'}, 'Data', {[] [] cell(1,2)}, 'Path', [], 'Loaded', false);
d2 = struct('Type', {'Activity', 'Consistency', 'Blank'}, 'Data', {[] [] cell(1,2)}, 'Path', [], 'Loaded', false);
handles.DataStruct = struct('Type', {'Pre-Metrics', 'Metrics'}, 'Data', {d1 d2}, 'Rows', [], 'Cols', [], 'RowSelect', [], 'ColSelect', []);
handles.SampleData = struct('IDs', [], 'Labels', [], 'Class', [], 'Value', [], 'Survival', [], 'Censor', []);
handles.ButtonEnable = {'on' 'on' 'off' 'on' 'off' 'on' 'off' 'on' 'on'};
handles.MainHandle = handles.figure1;
guidata(hObject, handles);

%read in PID database
fid = fopen('pathologist.db.txt');
fd = textscan(fid, '%s%f%s%s%f%s%s%s%s%s%s%s', 'delimiter', '\t');
fclose(fid);
handles.PwData.DB = horzcat(fd{1},num2cell(fd{2}),fd{3:4},num2cell(fd{5}),fd{6:end});

%read in PID complex data
fid = fopen('pathologist.complexes.txt');
fc = textscan(fid, '%f%s%s%f%s%s%s%s%s%f', 'delimiter', '\t');
fclose(fid);
handles.PwData.Complexes = horzcat(num2cell(fc{1}),fc{2:3},num2cell(fc{4}),fc{5:7});

%set up list of pathways
fid = fopen('pathologist.pathways.txt');
% pws = textscan(fid, '%10.f%s', 'delimiter', '\t', 'BufSize', 9000);
pws = textscan(fid, '%f%s', 'delimiter', '\t');
fclose(fid);
[pnames ix]= sort(lower(pws{2}));
pnums = pws{1}(ix);
handles.PwData.Names = pnames;
handles.PwData.Nums = pnums;

%read in uniprot-locuslink mapping
up2ll = importdata('UP2LL.txt', '\t');
handles.PwData.UP2LL = [up2ll.textdata num2cell(up2ll.data)];

clear fid fd fc pws pnames pnums up2ll

%set up platform data
handles.ChipData.Names = {
'Choose a microarray chip';
'HG-U133A';
'HG-U133B';
'HG-U133_Plus_2';
'HG-U95';
'HG-U133AAofAv2';
'MG-430A';
'MG-430B';
'Mouse430_2';
'Add new platform'};
handles.ChipData.Files = strcat('probelinks/', ...
    {'-';
    'U133a.txt';
    'U133b.txt';
    'U133plus2.txt';
    'U95.txt';
    'U133AA_of_Av2.txt';
    'Mouse430A.txt';
    'Mouse430B.txt';
    'Mouse4302.txt'});
handles.ChipData.Selection = 1;
handles.ChipData.Links = [];

guidata(hObject, handles);

function varargout = PathOlogist_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% Button Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function load_dataButton_Callback(hObject, eventdata, handles)
% hObject    handle to load_dataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cleared = 1;
while cleared

%open Load Data window
[H handles.DataStruct handles.SampleData handles.ChipData handles.ButtonEnable cleared next]= loaddata_gui('UserData', {handles.DataStruct, handles.SampleData, handles.ChipData, handles.ButtonEnable});
close(H);
if cleared
    handles.DataStruct(1).ColSelect=[];
    handles.SampleData = struct('IDs', [], 'Labels', [], 'Class', [], 'Value', [], 'Survival', [], 'Censor', []);
    guidata(handles.figure1, handles);
end
guidata(handles.figure1, handles);
end



%open Calculate Probabilities window
if next
    calc_udpButton_Callback(hObject, eventdata, handles)
end

 
function calc_pwButton_Callback(hObject, eventdata, handles)
% hObject    handle to calc_pwButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[H handles.DataStruct handles.SampleData next] = calcpw_gui('UserData', {handles.version handles.DataStruct , handles.SampleData, handles.PwData, handles.ChipData});
close(H)
guidata(handles.figure1, handles);
if sum([handles.DataStruct(2).Data(:).Loaded])
    set(handles.calculateCheck, 'Value', 1);
end
if next 
    choose_pwButton_Callback(hObject, eventdata, handles)
end

function draw_pwButton_Callback(hObject, eventdata, handles)
% hObject    handle to draw_pwButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[H drawn] = drawpw_gui('UserData', {handles.version handles.PwData, handles.DataStruct handles.ChipData});
close(H)
if drawn
    set(handles.drawCheck, 'Value', 1);
end

function sample_dataButton_Callback(hObject, eventdata, handles)
% hObject    handle to sample_dataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.SampleData.IDs) 
    handles.SampleData.IDs = handles.DataStruct(2).Cols';
end
[H handles.SampleData next] = sampledata_gui('UserData', handles.SampleData);
guidata(handles.figure1, handles);
close(H)
if next
    heatmapButton_Callback(hObject, eventdata, handles)
end

function heatmapButton_Callback(hObject, eventdata, handles)
% hObject    handle to heatmapButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[H next] = view_heatmaps('UserData', {handles.DataStruct; handles.SampleData});
close(H);
if next
    analyzeButton_Callback(hObject, eventdata, handles)
end

function analyzeButton_Callback(hObject, eventdata, handles)
% hObject    handle to analyzeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[H handles.DataStruct next] = analyzepws_gui('UserData', {handles.DataStruct; handles.SampleData; handles.PwData});
guidata(handles.figure1, handles);
close(H);
if next
    heatmapButton_Callback(hObject, eventdata, handles)
end
    
% Empty Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function checkbox7_Callback(hObject, eventdata, handles)
function checkbox8_Callback(hObject, eventdata, handles)
function udpCheck_Callback(hObject, eventdata, handles)
function loadCheck_Callback(hObject, eventdata, handles)
function drawCheck_Callback(hObject, eventdata, handles)
function samplesCheck_Callback(hObject, eventdata, handles)
function pwCheck_Callback(hObject, eventdata, handles)
function calculateCheck_Callback(hObject, eventdata, handles)
