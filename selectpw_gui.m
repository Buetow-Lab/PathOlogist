function varargout = selectpw_gui(varargin)
% SELECTPW_GUI M-file for selectpw_gui.fig
%      SELECTPW_GUI, by itself, creates a new SELECTPW_GUI or raises the existing
%      singleton*.
%
%      H = SELECTPW_GUI returns the handle to a new SELECTPW_GUI or the handle to
%      the existing singleton*.
%
%      SELECTPW_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTPW_GUI.M with the given input arguments.
%
%      SELECTPW_GUI('Property','Value',...) creates a new SELECTPW_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before selectpw_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to selectpw_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help selectpw_gui

% Last Modified by GUIDE v2.5 31-Jul-2008 12:43:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @selectpw_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @selectpw_gui_OutputFcn, ...
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


% --- Executes just before selectpw_gui is made visible.
function selectpw_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to selectpw_gui (see VARARGIN)

% Choose default command line output for selectpw_gui
u  = get(hObject, 'UserData');
[handles.DataStruct handles.PwData] = deal(u{:});
handles.output = hObject;
if isempty(handles.DataStruct(1).RowSelect)
    handles.DataStruct(1).RowSelect = 1;
end
handles.next = 0;
set(handles.pw_Listbox, 'String', handles.PwData.Names, 'Value', handles.DataStruct(1).RowSelect, 'Visible', 'on');


% Update handles structure
guidata(hObject, handles);
uiwait
% UIWAIT makes selectpw_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = selectpw_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.DataStruct;
varargout{3} = handles.next;


% --- Executes on selection change in pw_Listbox.
function pw_Listbox_Callback(hObject, eventdata, handles)
% hObject    handle to pw_Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pw_Listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pw_Listbox


% --- Executes during object creation, after setting all properties.
function pw_Listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pw_Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in backButton.
function backButton_Callback(hObject, eventdata, handles)
% hObject    handle to backButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataStruct(1).RowSelect= get(handles.pw_Listbox, 'Value');
guidata(hObject, handles);
uiresume

% --- Executes on button press in nextButton.
function nextButton_Callback(hObject, eventdata, handles)
% hObject    handle to nextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.next = 1;
guidata(handles.figure1, handles);
backButton_Callback(hObject, eventdata, handles)


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

if get(hObject, 'Value')
    set(handles.pw_Listbox, 'Value', 1:length(handles.PwData.Names));
else
    set(handles.pw_Listbox, 'Value', 1);
end
