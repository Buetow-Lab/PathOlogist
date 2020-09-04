function varargout = calcpw_gui(varargin)
% CALCPW_GUI M-file for calcpw_gui.fig
%      CALCPW_GUI, by itself, creates a new CALCPW_GUI or raises the existing
%      singleton*.
%
%      H = CALCPW_GUI returns the handle to a new CALCPW_GUI or the handle to
%      the existing singleton*.
%
%      CALCPW_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALCPW_GUI.M with the given input arguments.
%
%      CALCPW_GUI('Property','Value',...) creates a new CALCPW_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calcpw_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to calcpw_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help calcpw_gui

% Last Modified by GUIDE v2.5 17-Dec-2008 15:21:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calcpw_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @calcpw_gui_OutputFcn, ...
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
end

% --- Executes just before calcpw_gui is made visible.
function calcpw_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to calcpw_gui (see VARARGIN)
warning('off', 'all');
% Choose default command line output for drawpw_gui
handles.output = hObject;
handles.calculated = false;
handles.next=0;
u = get(hObject, 'UserData');
[handles.version handles.DataStruct handles.SampleData handles.PwData handles.ChipData] = deal(u{:});

if isempty(handles.DataStruct(1).RowSelect)
    handles.DataStruct(1).RowSelect = 1;
end
set(handles.pw_Listbox, 'String', handles.PwData.Names, 'Value', handles.DataStruct(1).RowSelect);

if isempty(handles.DataStruct(1).Cols)
    handles.samples_selected = 1;
    set(handles.samplesListbox,'String','Load data to see samples');
else
    set(handles.samplesListbox,'String',handles.DataStruct(1).Cols);
    if ~isempty(handles.DataStruct(1).ColSelect)
        set(handles.samplesListbox, 'Value', handles.DataStruct(1).ColSelect);
    end
end

if  handles.DataStruct(1).Data(2).Loaded
    set(handles.calculateButton, 'Enable', 'on');
end
guidata(hObject, handles);
uiwait
end

function varargout = calcpw_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.DataStruct;
varargout{3} = handles.SampleData;
varargout{4} = handles.next;
end

% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calculateButton_Callback(hObject, eventdata, handles)
% hObject    handle to calculateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DataStruct(1).RowSelect= get(handles.pw_Listbox, 'Value');
handles.DataStruct(1).ColSelect = get(handles.samplesListbox, 'Value');

set(handles.errorText, 'Visible', 'off')
set(hObject, 'BackgroundColor', 'white');
drawnow expose
try
%calculate activity and consistency for chosen pathway and uploaded data

pwnames = handles.PwData.Names(handles.DataStruct(1).RowSelect);
pwnums = handles.PwData.Nums(handles.DataStruct(1).RowSelect);

udp_s = handles.DataStruct(1).ColSelect;

p_activity = zeros(length(pwnums),length(udp_s));
p_consistency = zeros(length(pwnums),length(udp_s));

W = waitbar(0, 'Calculating... 0% complete', 'CreateCancelBtn', 'set(gcf, ''Clipping'', ''off'')');
PctComplete=0;
for pw = 1:length(pwnums)
    if strcmp(get(W, 'Clipping'), 'on')
        if round(pw*100/length(pwnums))>PctComplete
            PctComplete= round(pw*100/length(pwnums));
            waitbar(pw/length(pwnums), W, ['Calculating... ' num2str(PctComplete) '% complete']);
        end
    drawnow
    [molList molmap A] = get_molList(pwnums(pw), udp_s, handles);

%perform calculations
    [i_activity i_consistency omit_i omit_i_c] = pw_calc(A, molList, molmap, udp_s);
    group = ones(size(i_activity,1),1);
    group(omit_i) = 0;
    group = logical(group);

    p_activity(pw,:) = mean(i_activity(group,:),1);
    group(omit_i_c) = false;
    p_consistency(pw,:) = mean(i_consistency(group,:),1);
    else
        pw = length(pwnums);
        break;
    end
end

aborted = ~strcmp(get(W, 'Clipping'), 'on');
delete(W)

if  ~aborted
sample_labels = handles.DataStruct(1).Cols(udp_s);
row_labels = ['Pathway Name'; pwnames];

a = [row_labels [sample_labels; num2cell(p_activity)]];
c = [row_labels [sample_labels; num2cell(p_consistency)]];

[afile_name,afile_dir] = uiputfile('*.txt','Save activity as: ');
Directory= cd;
cd(afile_dir);
fid = fopen(afile_name, 'wt');
fprintf(fid, '%s\t', a{1,1:end-1});
fprintf(fid, '%s', a{1,end});
for rr = 2:size(a,1)
    fprintf(fid,'\n%s', a{rr,1});
    fprintf(fid, '\t%1.4f', a{rr,2:end});
end
fclose(fid);
cd(Directory);
            
[cfile_name,cfile_dir] = uiputfile('*.txt','Save consistency as: ');
Directory = cd;
cd(cfile_dir);
fid = fopen(cfile_name, 'wt');
fprintf(fid, '%s\t', c{1,1:end-1});
fprintf(fid, '%s', c{1,end});
for rr = 2:size(c,1)
    fprintf(fid,'\n%s', c{rr,1});
    fprintf(fid, '\t%1.4f', c{rr,2:end});
end
fclose(fid);
cd(Directory);
            
f_text = ['Metrics written to ' afile_name ' and ' cfile_name]; 
   
set(handles.errorText, 'Visible', 'on', 'String', f_text);
data_types = {'Activity' 'Consistency'; p_activity p_consistency; afile_dir cfile_dir; afile_name cfile_name};
for dt = 1:2
handles.DataStruct(2).Data(dt).Loaded = true;
handles.DataStruct(2).Data(dt).Data = data_types{2,dt};
handles.DataStruct(2).Rows= row_labels(2:end);
handles.DataStruct(2).Cols= sample_labels;
handles.DataStruct(2).Data(dt).Path = [data_types{3,dt} '\' data_types{4,dt}];
handles.SampleData.IDs = sample_labels;
end           
guidata(hObject, handles);
end
catch
    delete(W)
    set(handles.errorText, 'Visible', 'on')
end

set(hObject, 'String', 'Calculate');
end

function [molList molmap A] = get_molList(pwnum, udp_s, handles)

A = handles.PwData.DB([handles.PwData.DB{:,2}] == pwnum,:);
A(strcmp('null',A(:,4)),4)= A(strcmp('null',A(:,4)),6); 

D = dataset({A, 'pwName', 'pwNum', 'molType', 'molName', 'molNum', 'molLink', 'c7', 'c8', 'c9', 'molRole', 'c11', 'intxnType'});

%create list of mols involved in pathway, and a map to their rows in A
[N iaa molmap]= unique([A{:,5}]);
number_of_mols = length(N);
N = A(iaa,4);
Arows = iaa;
molList = D(Arows, [4 5 3 6]);

clear D
%for each mol in pathway, use LL (or UP)# to find up/down prob for each sample using
%probe_links table   
if ~isempty(udp_s)
    probe_links = handles.ChipData.Links;
for j = 1:number_of_mols
    class = get_class(molList.molType{j}, molList.molLink{j});  %classes: 1:protein, 2:complex, 3:compound, 4:other
    % j
    [molProbs(j,1:length(udp_s)) compProbs{j}]= get_probs(class, molList.molLink{j}, molList.molNum{j}, udp_s, []);
end
else
    molProbs(1:number_of_mols,1:length(udp_s)) = NaN;
    compProbs = cell(1,number_of_mols);
end
%add probabilities to list of mols
molList = [molList dataset({molProbs, 'molProbs'}, {compProbs', 'compProbs'}) ];

function class = get_class(type,link) %determine class of a molecule - 1:protein, 2:complex, 3:compound, 4:other
    % (code changed because both protein and complex can be calculated
    % using the same routine 2019-06-20)
        if strcmp(type, 'protein')
            class = 1;
        elseif strcmp(type, 'complex')
            class = 1;
        elseif strcmp(type, 'other')
            class = 1;
        elseif strcmp(type, 'compound')
            class = 3;
        else
            class = 4;
        end
end
function [probs cProbs]= get_probs(class, link, num, snums, cProbs)
    %use molecule link to find associated probes, and compile upd/down probabilities
    %for those probes into a probability for the molecule
   
        if class == 1
            cProbs = [];
            parsedlink = get_parsedlink(link);
            if ~isempty(parsedlink)
                proberows = get_proberows(parsedlink);
                if proberows
                    Probs = get_pr(proberows, snums);
                    probs(1,1:length(snums)) = 1;
                        for cc = 1:size(Probs,1)
                            probs = probs .* Probs(cc,:);
                        end
                else
                    probs = NaN;    %no matching proberows
                end
            else
                probs = NaN;    %invalid link
            end
            
        elseif class == 2 % this sub-routine is abandoned for complex probability calculation; complex probability be used as the same routine as protein 2019-06-20
 
            parsedlink = get_parsedlink(link);
            % complex_data =
            % handles.PwData.Complexes([handles.PwData.Complexes{:,1}] == num,:);   Don't know why use 1st column to find complex content (ID at 1st column isn't the mol ID 2019-06-20)
            % complex_data = handles.PwData.Complexes([handles.PwData.Complexes{:,4}] == num,:);
            % found_prob = 0;
            % probs(1,1:length(snums)) = 1;
            if ~isempty(parsedlink)
                proberows = get_proberows(parsedlink);
                if proberows
            % try
                % for cc = 1:size(complex_data,1)
                    % class = get_class(complex_data{cc,2}, complex_data{cc,5});
                    % cProbs = get_probs(class, complex_data{cc,5}, complex_data{cc,4}, snums, cProbs);
                    cProbs = get_pr(proberows, snums);
                    probs(1,1:length(snums)) = 1;
                    % if ~isnan(cProbs(1))
                        % found_prob = 1;
                        for cc = 1:size(cProbs,1)
                            probs = probs .* cProbs(cc,:);
                        end
                    % end
                % end
                % cProbs = [unique(complex_data(:,3)) num2cell(probs)];
                % if ~found_prob %if none of the components has a probability
                else
                    probs = NaN;    %no components probs found
                end
            % catch
            else
                probs = NaN;    %no component data
            end

        elseif class == 3
             cProbs = [];
            probs = 1;	%compounds are assumed to always be present
        elseif class == 4
             cProbs = [];
            probs = 1;    %others are assumed to always be present
        end
end   
function parsedlink = get_parsedlink(str)
    %turn linking string for molecule into useful search term in linking table
        L = findstr('LL:', str);
        U = findstr('UP:', str);
        C = findstr(',', str);

        if L
            for i = 1:length(L)
                start = L(i) + 3;
                if length(C)>=i
                    parsedlink{i,1} = (str(start:C(i) - 1));
                    parsedlink{i,2} = 2;
                else
                    parsedlink{i,1} = (str(start:end)); 
                    parsedlink{i,2} = 2;
                end
            end
            
        elseif U
            for i = 1:length(U)
                start = U(i) + 3;
                if length(C)>=i
                    parsedlink{i,1} = str(start:C(i) - 1);
                    parsedlink{i,2} = 3;
                else
                    parsedlink{i,1} = str(start:end);
                    parsedlink{i,2} = 3;
                end
            try
            parsedlink{i,2} = UP2LL{strmatch(parsedlink{i,1}, UP2LL(:,1)), 2};
            parsedlink{i,1} = num2str(parsedlink{i,2});
            catch
                parsedlink = [];
            end
            end
        else 
            parsedlink = [];
        end
end
function proberows = get_proberows(parsedlink)
    %search linking table for probes associated with LL# or UP#
    proberows = [];

    for pl = 1:size(parsedlink,1)
        try
            if parsedlink{pl,2} ==3
                parsedlink{pl,1} = num2str(handles.PwData.UP2LL{strmatch(parsedlink{pl,1}, handles.PwData.UP2LL(:,1)), 2});
            end
            pnames = probe_links(strmatch(parsedlink{pl,1}, probe_links(:,2), 'exact'),1);
            for pn = 1:length(pnames)
                proberows = [proberows; strmatch(pnames{pn}, handles.DataStruct(1).Rows, 'exact');];
            end
        catch
        end
    end
        
end
function pr = get_pr(proberows, snums)
        for n = 1:length(proberows) %get probabilities from UDP for each associated probe
            p(n,:) = handles.DataStruct(1).Data(2).Data(proberows(n),snums);
        end
        % pr = max(p, [], 1); %for each sample, take the max probability in
        % the probeset (code changed because using RNAseq data not microarray data 2019-06-20)
        pr = p;
end
end

function [i_activity i_consistency omit_i omit_i_c mult_out] = pw_calc(A, molList, molmap, udp_s)
[intxn_levels ia intxn_map]= unique(A(:,11));
number_of_intxns = length(intxn_levels);

%set up variables
i_activity = ones(number_of_intxns,length(udp_s)); 
i_consistency = ones(number_of_intxns,length(udp_s)); 
omit_i = [];
omit_i_c = [];
% skip_rest = 0;
mult_out = [];
%for each interaction, use mol roles (eg. input/output) and probabilities
%to calculate activity and consistency
for in = 1:number_of_intxns 
    skip_rest = 0;
    intxn_group = (intxn_map == in);
    m_roles = A(intxn_group,10);
    m_probs = molList.molProbs(molmap(intxn_group),:);
    m_index = 1:size(A,1);
    m_index = m_index(intxn_group);
    outs = strcmp(m_roles, 'output');
    if sum(outs)>1
        mult_out_i = [in*ones(sum(outs),1) m_index(outs)'];
    end    
    for jj = 1:length(m_roles)
        if isnan(m_probs(jj,1))
            skip_rest = 1; %if a mol does not have a probability, skip the interaction
            omit_i = [omit_i in];
        end
    end
            
    if skip_rest ~= 1
        for j = 1:length(m_roles)
            if ( strcmp(m_roles{j},'agent')||strcmp(m_roles{j},'input')||strcmp(m_roles{j},'regulator') )
                i_activity(in,:) = i_activity(in,:).* m_probs(j,:);
            elseif strcmp(m_roles{j}, 'inhibitor')
                i_activity(in,:) = i_activity(in,:).* (1-m_probs(j,:));
            end
        end
    
       
       clear i_con
       if sum(outs)>0
       for j = 1:length(udp_s)
          i_con(:,j) = i_activity(in,j).* m_probs(outs,j) + (1 - i_activity(in,j)).*(1-m_probs(outs,j));
       end
       else
           omit_i_c = [omit_i_c in];
           i_con(:,1:length(udp_s)) = NaN;
       end
       if size(i_con,1)>1
           for ss = 1:size(i_con,2)
                i_consistency(in,ss) = mean(i_con(~isnan(i_con(:,ss)),ss));
                
           end
           mult_out_i = [mult_out_i i_con];
           mult_out = [mult_out; mult_out_i];
       else
           i_consistency(in,:) = i_con;
           mult_out_i = [];
       end
    
    else
      if sum(outs)>1
          mult_out = [mult_out; [mult_out_i NaN*ones(size(mult_out_i,1),length(udp_s))]];
      end
    end
    
end
end

function nextButton_Callback(hObject, eventdata, handles)
% hObject    handle to nextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.next=1;
guidata(hObject, handles);
returnButton_Callback(hObject, eventdata, handles)
end

function returnButton_Callback(hObject, eventdata, handles)
% hObject    handle to returnButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.DataStruct(2).Data(1).Loaded
    handles.calculated = true;
    guidata(hObject, handles);
end

uiresume
end


% --- Executes on selection change in samplesListbox.
function samplesListbox_Callback(hObject, eventdata, handles)
% hObject    handle to samplesListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns samplesListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from samplesListbox
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

end
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
end
% --- Executes on selection change in pw_Listbox.
function pw_Listbox_Callback(hObject, eventdata, handles)
% hObject    handle to pw_Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pw_Listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pw_Listbox
end

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
end

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
end
