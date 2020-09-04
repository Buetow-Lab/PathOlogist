function varargout = view_pws(varargin)
% VIEW_PWS M-file for view_pws.fig
%      VIEW_PWS, by itself, creates a new VIEW_PWS or raises the existing
%      singleton*.
%
%      H = VIEW_PWS returns the handle to a new VIEW_PWS or the handle to
%      the existing singleton*.
%
%      VIEW_PWS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW_PWS.M with the given input arguments.
%
%      VIEW_PWS('Property','Value',...) creates a new VIEW_PWS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before view_pws_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to view_pws_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help view_pws

% Last Modified by GUIDE v2.5 13-May-2019 16:06:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @view_pws_OpeningFcn, ...
                   'gui_OutputFcn',  @view_pws_OutputFcn, ...
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


% --- Executes just before view_pws is made visible.
function view_pws_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to view_pws (see VARARGIN)

% Choose default command line output for view_pws
handles.output = hObject;

u = get(hObject, 'UserData');
[handles.choices handles.selected] = deal(u{:});
if isempty(handles.selected)
    handles.selected = 1;
end
set(handles.listbox1, 'String', handles.choices, 'Value', handles.selected);
% Update handles structure
guidata(hObject, handles);
uiwait
% UIWAIT makes view_pws wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = view_pws_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.selected;

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.selected=get(handles.listbox1, 'Value');
guidata(hObject, handles);
uiresume
