function varargout = RMA_gui(varargin)
% RMA_GUI M-file for RMA_gui.fig
%      RMA_GUI, by itself, creates a new RMA_GUI or raises the existing
%      singleton*.
%
%      H = RMA_GUI returns the handle to a new RMA_GUI or the handle to
%      the existing singleton*.
%
%      RMA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RMA_GUI.M with the given input arguments.
%
%      RMA_GUI('Property','Value',...) creates a new RMA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RMA_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RMA_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RMA_gui

% Last Modified by GUIDE v2.5 27-May-2008 17:32:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RMA_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @RMA_gui_OutputFcn, ...
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


% --- Executes just before RMA_gui is made visible.
function RMA_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RMA_gui (see VARARGIN)

% Choose default command line output for RMA_gui
handles.output = hObject;
handles.loaded = [0 0];
handles.DataStruct = get(handles.figure1, 'UserData');
handles.samples = [];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RMA_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
uiwait

% --- Outputs from this function are returned to the command line.
function varargout = RMA_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.DataStruct;

% --- Executes on button press in celButton.
function celButton_Callback(hObject, eventdata, handles)
% hObject    handle to celButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cel_pathname = uigetdir('','Select the .cel folder');
if handles.cel_pathname
    set(handles.celText,'String', handles.cel_pathname, 'BackgroundColor', 'yellow' );
end
handles.loaded(1) = 1;
D = dir(handles.cel_pathname);
files = {D(3:end).name};
for ii = 1:length(files)
cel_files(ii) = strcmp(files{ii}(end-3:end),'.CEL') || strcmp(files{ii}(end-3:end),'.cel');
end
handles.cels = files(cel_files); 
guidata(hObject, handles);

% --- Executes on button press in cdfButton.
function cdfButton_Callback(hObject, eventdata, handles)
% hObject    handle to cdfButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.cdf_filename,handles.cdf_pathname] = uigetfile('*.cdf','Select the .cdf file');
if handles.cdf_filename
    set(handles.cdfText,'String', [handles.cdf_pathname '\' handles.cdf_filename], 'BackgroundColor', 'yellow' );
end
handles.loaded(2) = 1;
guidata(hObject, handles);

% --- Executes on button press in calculateButton.
function calculateButton_Callback(hObject, eventdata, handles)
% hObject    handle to calculateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if sum(handles.loaded) == 2
    set(handles.progressText, 'String', '');
try
    W = waitbar(0, 'Reading Files...', 'CreateCancelBtn', 'set(gcf, ''UserData'', false);', 'UserData', true);
    %set(handles.progressText, 'String', );
    drawnow expose
    probeData = celintensityread(handles.cels, handles.cdf_filename,'celpath', handles.cel_pathname, 'cdfpath', handles.cdf_pathname, 'pmonly', true);
    chipName = ProbeStructure.CDFName(1:strfind(ProbeStructure.CDFName,'.')-1);
    if get(W, 'UserData')
        waitbar(.25, W, 'Performing Background Adjustment...');
    %set(handles.progressText, 'String', 'Performing Background Adjustment...');
    drawnow expose
    pms_bg = rmabackadj(probeData.PMIntensities);
    end
    if get(W, 'UserData')
            waitbar(.5, W, 'Performing Quantile Normalization...');
    %set(handles.progressText, 'String', 'Performing Quantile Normalization...');
    drawnow expose
    pms_bgnorm = quantilenorm(pms_bg);
    end
    if get(W, 'UserData')
       waitbar(.75, W, 'Summarizing Results...');
    %set(handles.progressText, 'String', 'Summarizing Results...');
    %drawnow expose
    rma_exp = rmasummary(probeData.ProbeIndices, pms_bgnorm, 'Output', 'natural');
    end
    aborted = ~get(W, 'UserData');
    delete(W)
    if~aborted
    [filename pathname] = uiputfile('*.txt', 'Save RMA as:');
    
    RMA_data = [[{chipName} probeData.CELNames'];[probeData.ProbeSetIDs num2cell(rma_exp)]];
    fid = fopen([pathname '\' filename], 'wt');
    fprintf(fid, '%s\t', RMA_data{1,1:end-1});
    fprintf(fid, '%s', RMA_data{1,end});
    for ii = 2:size(rma_exp,1)
       fprintf(fid, '\n%s', RMA_data{ii,1});
       fprintf(fid, '\t%10.6f', RMA_data{ii,2:end});
    end
    fclose(fid)
    
    set(handles.progressText, 'String', 'RMA values have been calculated and written to');
    drawnow expose
    
    handles.DataStruct(1).Data(1).Loaded = true;
    handles.DataStruct(1).Rows= probeData.ProbeSetIDs;
    handles.DataStruct(1).Data(1).Data = rma_exp;
    handles.DataStruct(1).Data(1).Path = [pathname '\' filename];
    handles.DataStruct(1).Cols = probeData.CELNames';
    guidata(hObject, handles);

    uiresume
    end
catch
    delete(W)
    set(handles.progressText, 'String', 'Error: Ensure that chosen folder contains .cel files corresponding to chosen .cdf file');
end
    else
    set(handles.progressText, 'String', 'Please enter all required information.');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume
