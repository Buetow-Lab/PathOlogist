function varargout = entersampledata_gui(varargin)
% ENTERSAMPLEDATA_GUI M-file for entersampledata_gui.fig
%      ENTERSAMPLEDATA_GUI, by itself, creates a new ENTERSAMPLEDATA_GUI or raises the existing
%      singleton*.
%
%      H = ENTERSAMPLEDATA_GUI returns the handle to a new ENTERSAMPLEDATA_GUI or the handle to
%      the existing singleton*.
%
%      ENTERSAMPLEDATA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ENTERSAMPLEDATA_GUI.M with the given input arguments.
%
%      ENTERSAMPLEDATA_GUI('Property','Value',...) creates a new ENTERSAMPLEDATA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before entersampledata_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to entersampledata_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help entersampledata_gui

% Last Modified by GUIDE v2.5 16-Sep-2008 16:40:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @entersampledata_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @entersampledata_gui_OutputFcn, ...
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


% --- Executes just before entersampledata_gui is made visible.
function entersampledata_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to entersampledata_gui (see VARARGIN)

% Choose default command line output for entersampledata_gui
handles.output = hObject;
u = get(hObject, 'UserData');
[handles.Samples, handles.data] = deal(u{:});
set(handles.samplesListbox, 'String', handles.Samples);
if ~isempty(handles.data)
    set(handles.dataText, 'String', handles.data);
    set(handles.entrycount, 'String', size(handles.data,1));
end
set(handles.samplecount, 'String', length(handles.Samples));
% Update handles structure
guidata(hObject, handles);
uiwait
% UIWAIT makes entersampledata_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = entersampledata_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.data;



function dataText_Callback(hObject, eventdata, handles)
% hObject    handle to dataText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataText as text
%        str2double(get(hObject,'String')) returns contents of dataText as a double
set(handles.entrycount, 'String', size(get(hObject, 'String'),1));

% --- Executes during object creation, after setting all properties.
function dataText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if size(get(handles.dataText, 'String'),1)==length(handles.Samples)
handles.data = cellstr(get(handles.dataText, 'String'));
guidata(hObject, handles);
uiresume

else
    set(handles.error, 'Visible', 'on');
end




% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume


% --- Executes on button press in clearButton.
function clearButton_Callback(hObject, eventdata, handles)
% hObject    handle to clearButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data = [];
guidata(hObject, handles);
uiresume


% --- Executes on selection change in samplesListbox.
function samplesListbox_Callback(hObject, eventdata, handles)
% hObject    handle to samplesListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns samplesListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from samplesListbox
try
set(handles.dataText, 'Value', get(hObject, 'Value'));
end

% --- Executes during object creation, after setting all properties.
function samplesListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to samplesListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


