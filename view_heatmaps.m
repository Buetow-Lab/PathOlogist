function varargout = view_heatmaps(varargin)
set(0,'RecursionLimit',1000)
% VIEW_HEATMAPS M-filsete for view_heatmaps.fig
%      VIEW_HEATMAPS, by itself, creates a new VIEW_HEATMAPS or raises the existing
%      singleton*.
%
%      H = VIEW_HEATMAPS returns the handle to a new VIEW_HEATMAPS or the handle to
%      the existing singleton*.
%
%      VIEW_HEATMAPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW_HEATMAPS.M with the given input arguments.
%
%      VIEW_HEATMAPS('Property','Value',...) creates a new VIEW_HEATMAPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before view_heatmaps_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to view_heatmaps_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help view_heatmaps

% Last Modified by GUIDE v2.5 05-Nov-2019 12:00:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @view_heatmaps_OpeningFcn, ...
                   'gui_OutputFcn',  @view_heatmaps_OutputFcn, ...
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


% --- Executes just before view_heatmaps is made visible.
function view_heatmaps_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to view_heatmaps (see VARARGIN)

% Choose default command line output for view_heatmaps
handles.output = hObject;
u = get(hObject, 'UserData');
[handles.DataStruct handles.SampleData] = deal(u{:});
if ~isempty(handles.SampleData.Class)
    [classes ia ib] = unique(handles.SampleData.Class);
    set(handles.classListbox, 'String', classes, 'UserData', ib);
else
    set(handles.sample_classRadio, 'Enable', 'off')
    set(handles.selected_samplesRadio, 'Value', 1)
end

if handles.DataStruct(2).Data(1).Loaded
    set(handles.actCheck, 'Enable', 'on', 'Value', 1)
end
if handles.DataStruct(2).Data(2).Loaded
    set(handles.conCheck, 'Enable', 'on', 'Value', 1)
end
handles.next=0;

% Update handles structure
guidata(hObject, handles);
uiwait
% UIWAIT makes view_heatmaps wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = view_heatmaps_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.next;


% --- Executes on button press in actCheck.
function actCheck_Callback(hObject, eventdata, handles)
% hObject    handle to actCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of actCheck


% --- Executes on button press in conCheck.
function conCheck_Callback(hObject, eventdata, handles)
% hObject    handle to conCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of conCheck


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in select_samplesButton.
function select_samplesButton_Callback(hObject, eventdata, handles)
% hObject    handle to select_samplesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.SampleData.Labels)
    choices = handles.SampleData.Labels;
else
    choices = handles.DataStruct(2).Cols;
end
[H handles.DataStruct(2).ColSelect] = view_pws('UserData', {choices handles.DataStruct(2).ColSelect});
close(H)
guidata(hObject, handles);


% --- Executes on selection change in classListbox.
function classListbox_Callback(hObject, eventdata, handles)
% hObject    handle to classListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns classListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from classListbox


% --- Executes during object creation, after setting all properties.
function classListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function variable_pwsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to variable_pwsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of variable_pwsEdit as text
%        str2double(get(hObject,'String')) returns contents of variable_pwsEdit as a double


% --- Executes during object creation, after setting all properties.
function variable_pwsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to variable_pwsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in select_pwsButton.
function select_pwsButton_Callback(hObject, eventdata, handles)
% hObject    handle to select_pwsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[H handles.DataStruct(2).RowSelect] = view_pws('UserData', {handles.DataStruct(2).Rows handles.DataStruct(2).RowSelect});
close(H)
guidata(hObject, handles);

% --- Executes on button press in returnButton.
function returnButton_Callback(hObject, eventdata, handles)
% hObject    handle to returnButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume

% --- Executes on button press in nextButton.
function nextButton_Callback(hObject, eventdata, handles)
% hObject    handle to nextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.next=1;
guidata(hObject, handles);

uiresume

% --- Executes on button press in viewButton.
function viewButton_Callback(hObject, eventdata, handles)
% hObject    handle to viewButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.error, 'Visible', 'off');
PwCount = length(handles.DataStruct(2).Rows);
row_labels = [strcat('act:', handles.DataStruct(2).Rows); strcat('con:', handles.DataStruct(2).Rows)];
col_labels = handles.SampleData.IDs;

if handles.DataStruct(2).Data(1).Loaded
data(1:PwCount,:) = handles.DataStruct(2).Data(1).Data;
end
if handles.DataStruct(2).Data(2).Loaded
data(1+PwCount:2*PwCount,:) = handles.DataStruct(2).Data(2).Data;
else
    data(1+PwCount:2*PwCount,:) = 0;
end

if get(handles.sample_classRadio, 'Value')
    ClassSelected = get(handles.classListbox, 'Value');
    IX =  get(handles.classListbox, 'UserData');
    s_ix = ismember(IX, ClassSelected);
    col_labels = strcat(handles.SampleData.Class, '(', col_labels, ')');
elseif get(handles.selected_samplesRadio, 'Value')
    s_ix = handles.DataStruct(2).ColSelect;
else
    s_ix = 1:size(data,2);
end

col_labels = col_labels(s_ix);
data = data(:,s_ix);

if get(handles.variable_pwsRadio, 'Value')
    try
        pct = str2double(get(handles.variable_pwsEdit, 'String'));
        if pct>100 || pct<0
            set(pw_error, 'Visible', 'on')
        else
            pw_ix = [];
            if get(handles.actCheck, 'Value')
                pw_ix = 1:PwCount;
                V = var(data(pw_ix,:)');
                [v_sort v_ix] = sort(V,'descend');
                pw_ix = v_ix(find(~isnan(v_sort),1):end);
                pw_ix = pw_ix(1:round(pct*size(data(pw_ix,:),1)/100));
            end
            if get(handles.conCheck, 'Value')
                pw_ix = [pw_ix (1:PwCount)+PwCount];
                V = var(data(pw_ix,:)');
                [v_sort v_ix] = sort(V,'descend');
                pw_ix = v_ix(find(~isnan(v_sort),1):end);
                pw_ix = pw_ix(1:round(pct*size(data(pw_ix,:),1)/100));
                pw_ix = pw_ix+PwCount;
            end
        end
    catch
        set(handles.error, 'Visible', 'on')
    end
elseif get(handles.selected_pwsRadio, 'Value') && ~isempty(handles.DataStruct(2).RowSelect)
    pw_ix = [];
    if get(handles.actCheck, 'Value')
        pw_ix = [pw_ix handles.DataStruct(2).RowSelect];
    end
    if get(handles.conCheck, 'Value')
        pw_ix = [pw_ix handles.DataStruct(2).RowSelect+PwCount];
    end
else
    pw_ix = [];
    if get(handles.actCheck, 'Value')
        pw_ix = 1:PwCount;
    end
    if get(handles.conCheck, 'Value')
        pw_ix = [pw_ix (1:PwCount)+PwCount];
    end
end

nan_data = isnan(data);
invalid = find(sum(nan_data')>0);
pw_ix(ismember(pw_ix,invalid)) = [];
row_labels = row_labels(pw_ix);
data = data(pw_ix,:);



if ~isempty(data)
    try
        % dlmwrite('HCC_matrix.txt', data, 'delimiter','\t','precision',5);
        % clustergram(data, 'RowLabels', row_labels, 'Dendrogram', {'ColorThreshold', 'default'}, 'Symmetric', false);
        clustergram(data, 'RowLabels', row_labels, 'ColumnLabels', col_labels, 'Linkage','ward','Colormap', jet(256), 'Symmetric', false);
        % clustergram(data, 'RowLabels', row_labels, 'Colormap', redbluecmap, 'Symmetric', false);
        if sum(invalid)>0
            % set(handles.error, 'Visible', 'on', 'String', 'Note: Some selected pathways had invalid values and were not included');
        end
    catch
        % set(handles.error, 'Visible', 'on', 'String', 'Not enough variation to perform clustering.');
    end
else
    set(handles.error, 'Visible', 'on');
end


% --- Executes on button press in select_pwsButton.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to select_pwsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
