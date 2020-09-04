function varargout = selectsample_gui(varargin)
% SELECTSAMPLE_GUI M-file for selectsample_gui.fig
%      SELECTSAMPLE_GUI, by itself, creates a new SELECTSAMPLE_GUI or raises the existing
%      singleton*.
%
%      H = SELECTSAMPLE_GUI returns the handle to a new SELECTSAMPLE_GUI or the handle to
%      the existing singleton*.
%
%      SELECTSAMPLE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTSAMPLE_GUI.M with the given input arguments.
%
%      SELECTSAMPLE_GUI('Property','Value',...) creates a new SELECTSAMPLE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before selectsample_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to selectsample_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help selectsample_gui

% Last Modified by GUIDE v2.5 21-Aug-2008 15:27:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @selectsample_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @selectsample_gui_OutputFcn, ...
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


% --- Executes just before selectsample_gui is made visible.
function selectsample_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to selectsample_gui (see VARARGIN)

% Choose default command line output for selectsample_gui
handles.output = hObject;
handles.next = 0;
handles.DataStruct = get(hObject, 'UserData');

if isempty(handles.DataStruct(1).Cols)
    handles.samples_selected = 1;
    set(handles.samplesListbox,'String','Load data to see samples');
else
    %if~isempty(handles.SampleData.Class)
     %   set(handles.classRadio,'Enable','on');
    %end
        set(handles.samplesListbox,'String',handles.DataStruct(1).Cols, 'Value', handles.DataStruct(1).ColSelect);
end
% Update handles structure
guidata(hObject, handles);
%set(handles.samplesListbox, 'Value', handles.samples_selected); %set initial selected samples

uiwait
% UIWAIT makes selectsample_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = selectsample_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.DataStruct;
varargout{3} = handles.next;


% --- Executes on selection change in samplesListbox.
function samplesListbox_Callback(hObject, eventdata, handles)
% hObject    handle to samplesListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns samplesListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from samplesListbox


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



% --- Executes on button press in returnButton.
function returnButton_Callback(hObject, eventdata, handles)
% hObject    handle to returnButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if ~isempty(handles.DataStruct(1).Cols)
       handles.DataStruct(1).ColSelect = get(handles.samplesListbox, 'Value');
    end

guidata(hObject, handles);

uiresume


% --- Executes on button press in allCheck.
function allCheck_Callback(hObject, eventdata, handles)
% hObject    handle to allCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allCheck
if get(hObject, 'Value')
    set(handles.samplesListbox, 'Value', 1:length(get(handles.samplesListbox, 'String')));
else
    set(handles.samplesListbox, 'Value', 1);
end

% --- Executes when selected object is changed in uipanel3.
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel3 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject == handles.nameRadio
    set(handles.samplesListbox, 'String', handles.SampleData.IDs)
elseif hObject == handles.labelRadio
    set(handles.samplesListbox, 'String', handles.SampleData.Labels)
elseif hObject == handles.classRadio
    [c ia ib] = unique(handles.SampleData.Class);
    set(handles.samplesListbox, 'String', c, 'Value', 1);
    set(handles.classRadio, 'UserData', ib);
end


% --- Executes on button press in nextButton.
function nextButton_Callback(hObject, eventdata, handles)
% hObject    handle to nextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.next=1;
returnButton_Callback(handles.returnButton, eventdata, handles);
