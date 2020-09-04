function varargout = sampledata_gui(varargin)
% SAMPLEDATA_GUI M-file for sampledata_gui.fig
%      SAMPLEDATA_GUI, by itself, creates a new SAMPLEDATA_GUI or raises the existing
%      singleton*.
%
%      H = SAMPLEDATA_GUI returns the handle to a new SAMPLEDATA_GUI or the handle to
%      the existing singleton*.
%
%      SAMPLEDATA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SAMPLEDATA_GUI.M with the given input arguments.
%
%      SAMPLEDATA_GUI('Property','Value',...) creates a new SAMPLEDATA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sampledata_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sampledata_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sampledata_gui

% Last Modified by GUIDE v2.5 05-Sep-2008 14:03:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sampledata_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @sampledata_gui_OutputFcn, ...
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


% --- Executes just before sampledata_gui is made visible.
function sampledata_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sampledata_gui (see VARARGIN)

% Choose default command line output for sampledata_gui
handles.output = hObject;
handles.SampleData = get(hObject, 'UserData');
if~isempty(handles.SampleData.IDs)
set(handles.IDsListbox, 'String', handles.SampleData.IDs');
else
    set(handles.pushbutton11, 'Enable', 'off')
end
set(handles.uipanel3, 'UserData', handles.labelsRadio);
handles.data_types = fieldnames(handles.SampleData);
for ii = 2:length(handles.data_types)
    eval(['set(handles.' handles.data_types{ii} 'Listbox, ''String'', handles.SampleData.' handles.data_types{ii} ');']);
end
handles.next = 0;
% Update handles structure
guidata(hObject, handles);
uiwait
% UIWAIT makes sampledata_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sampledata_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.SampleData;
varargout{3} = handles.next;

% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for ii = 1:length(handles.data_types)
    eval(['flipped=size(handles.SampleData.' handles.data_types{ii} ',2)>size(handles.SampleData.' handles.data_types{ii} ',1);']);
    if flipped
        eval(['handles.SampleData.' handles.data_types{ii} '= handles.SampleData.' handles.data_types{ii} ''';']);
    end
end
guidata(hObject, handles);

uiresume

% --- Executes on selection change in IDsListbox.
function IDsListbox_Callback(hObject, eventdata, handles)
% hObject    handle to IDsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns IDsListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from IDsListbox
ix = get(hObject,'Value');
for ii = 1:length(handles.data_types)
    try
    eval(['set(handles.' handles.data_types{ii} 'Listbox, ''Value'', ix);']);
    end
end

% --- Executes during object creation, after setting all properties.
function IDsListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IDsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in LabelsListbox.
function LabelsListbox_Callback(hObject, eventdata, handles)
% hObject    handle to LabelsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns LabelsListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LabelsListbox


% --- Executes during object creation, after setting all properties.
function LabelsListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LabelsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ClassListbox.
function ClassListbox_Callback(hObject, eventdata, handles)
% hObject    handle to ClassListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ClassListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ClassListbox


% --- Executes during object creation, after setting all properties.
function ClassListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ClassListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ValueListbox.
function ValueListbox_Callback(hObject, eventdata, handles)
% hObject    handle to ValueListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ValueListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ValueListbox


% --- Executes during object creation, after setting all properties.
function ValueListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ValueListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TimeListbox.
function TimeListbox_Callback(hObject, eventdata, handles)
% hObject    handle to TimeListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns TimeListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TimeListbox


% --- Executes during object creation, after setting all properties.
function TimeListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in SurvivalListbox.
function SurvivalListbox_Callback(hObject, eventdata, handles)
% hObject    handle to SurvivalListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns SurvivalListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SurvivalListbox


% --- Executes during object creation, after setting all properties.
function SurvivalListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SurvivalListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CensorListbox.
function CensorListbox_Callback(hObject, eventdata, handles)
% hObject    handle to CensorListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns CensorListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CensorListbox


% --- Executes during object creation, after setting all properties.
function CensorListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CensorListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%data = cell(length(get(handles.IDsListbox, 'String')),1);
DataType = get(get(handles.uipanel3, 'UserData'), 'UserData');
%DataTypeIX = strmatch(DataType,handles.data_types);
Samples = get(handles.IDsListbox, 'String');
eval(['data = handles.SampleData.' DataType ';']); 
[H data] = entersampledata_gui('UserData', {Samples,data});
close(H)

if~isempty(data)
if ~sum(strcmp(DataType, {'Labels' 'Class'}))
    for ii = 1:length(data)
        try
        data2(ii,1) = str2double(data(ii,:));
        catch
            data2(ii,1) = NaN;
        end
    end
    data = data2;
else
    if strcmp(DataType, 'Labels')
        handles.SampleData.IDs = data;
    end
end
end
eval(['set(handles.' DataType 'Listbox, ''String'', data);']);
eval(['handles.SampleData.' DataType '= data;']);

guidata(hObject, handles);


% --- Executes when selected object is changed in uipanel3.
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel3 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanel3, 'UserData', hObject);


    


% --- Executes on button press in nextButton.
function nextButton_Callback(hObject, eventdata, handles)
% hObject    handle to nextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.next=1;
guidata(hObject, handles);

pushbutton16_Callback(hObject, eventdata, handles)
