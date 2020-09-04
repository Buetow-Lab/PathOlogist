function varargout = analyzepws_gui(varargin)
% ANALYZEPWS_GUI M-file for analyzepws_gui.fig
%      ANALYZEPWS_GUI, by itself, creates a new ANALYZEPWS_GUI or raises the existing
%      singleton*.
%
%      H = ANALYZEPWS_GUI returns the handle to a new ANALYZEPWS_GUI or the handle to
%      the existing singleton*.
%
%      ANALYZEPWS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYZEPWS_GUI.M with the given input arguments.
%
%      ANALYZEPWS_GUI('Property','Value',...) creates a new ANALYZEPWS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before analyzepws_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to analyzepws_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help analyzepws_gui

% Last Modified by GUIDE v2.5 27-Apr-2016 12:25:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analyzepws_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @analyzepws_gui_OutputFcn, ...
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

function analyzepws_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to analyzepws_gui (see VARARGIN)

% Choose default command line output for analyzepws_gui
handles.output = hObject;
u = get(hObject, 'UserData');
[handles.DataStruct handles.SampleData handles.PwData] = deal(u{:});
handles.next=0;
if ~isempty(handles.SampleData.Class)
    set(handles.binaryRadio, 'Enable', 'on')
    [classes ia ib] = unique(handles.SampleData.Class);
    set(handles.group1Listbox, 'String', classes, 'UserData', ib);
    set(handles.corr_classListbox, 'String', classes, 'UserData', ib, 'Value', 1:length(classes));
end

if ~isempty(handles.SampleData.Value)
    set(handles.corrRadio, 'Enable', 'on');
end

if ~isempty(handles.SampleData.Survival) && ~isempty(handles.SampleData.Censor)
    set(handles.survivalRadio, 'Enable', 'on');
    set(handles.minsamplesPopup, 'String', num2cell(2:length(handles.SampleData.Survival)/2-1));
end

if handles.DataStruct(1).Data(3).Loaded
    set(handles.hitsRadio, 'Enable', 'on');
end

a_pws = []; c_pws = [];
for ii = 1:length(handles.DataStruct(2).Rows)
    if handles.DataStruct(2).Data(1).Loaded
    a_pws{ii} = strcat('act:', handles.DataStruct(2).Rows{ii});
    end
    if handles.DataStruct(2).Data(2).Loaded
    c_pws{ii} = strcat('con:', handles.DataStruct(2).Rows{ii});
    end
end

handles.pws = [a_pws';c_pws'];
handles.metrics = [handles.DataStruct(2).Data(1).Data; handles.DataStruct(2).Data(2).Data];
handles.bgcolor = get(handles.writeButton, 'BackgroundColor');           
guidata(hObject, handles);
uiwait

function varargout = analyzepws_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.DataStruct;
varargout{3} = handles.next;

% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function goButton_Callback(hObject, eventdata, handles)
% hObject    handle to goButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set(handles.goButton, 'String', 'Working...', 'FontSize', 8);
    drawnow expose
try  
    switch get(handles.test_typePanel, 'SelectedObject')
        %t test
        case handles.binaryRadio
        if strcmp(get(handles.group2Listbox, 'Enable'), 'on')
            Class2IX = get(handles.group2Listbox, 'UserData');
            ClassSelected = Class2IX(get(handles.group2Listbox, 'Value'));
            IX =  get(handles.group1Listbox, 'UserData');
            group2 = ismember(IX, ClassSelected);
            
            metrics = [handles.DataStruct(2).Data(1).Data; handles.DataStruct(2).Data(2).Data];
            data1 = metrics(:,handles.group1);
            data2 = metrics(:,group2);
            data = [data1 data2];
            
            group = [ones(sum(handles.group1),1); 2*ones(sum(group2),1)];
            group1classes=unique(handles.SampleData.Class(handles.group1));
            group2classes=unique(handles.SampleData.Class(group2));
            if length(group1classes)>1
                group1label = strcat(group1classes(1:end-1),'/');
                group1label = [group1label{:} group1classes{end}];
            else
                group1label = group1classes{1};
            end
            if length(group2classes)>1
                group2label = strcat(group2classes(1:end-1),'/');
                group2label = [group2label{:} group2classes{end}];
            else
                group2label = group2classes{1};
            end
            
            %[feature_set] = loofs(data', group, 5, handles.pws);
            for ii = 1:size(data,1)
                if var(data(ii,:))>eps
                    try
                        [h p(ii)] = ttest2(data1(ii,:), data2(ii,:));
                    catch
                        p(ii) = NaN;
                    end
                else
                    p(ii)=NaN;
                end
            end
            
            % [psort pix] = sort(p, 'ascend');
            len = length(p);
            p_adj(1:len/2) = p(1:len/2)*sum(~isnan(p(1:len/2)));
            p_adj((len/2)+1:len) = p((len/2)+1:len)*sum(~isnan(p((len/2)+1:len)));
            p_adj(find(p_adj>=1)) = 1;
            
            [p_nom pix] = sort(p, 'ascend');
            % add adj p value (bonferroni correction basis 2019-05-09)
            % p_adj = p_nom*(length(p_nom)/2);
            p_adj = p_adj(pix);
             
            
            set(handles.accuracyText, 'String', num2cell(p_nom));
            
           RowHeaders = cell(2,max(sum(handles.group1), sum(group2))+1);
           RowHeaders(1,1:sum(handles.group1)+1) = [['Class:' group1label]  handles.SampleData.IDs(handles.group1)'];
           RowHeaders(2,1:sum(group2)+1) = [['Class:' group2label]  handles.SampleData.IDs(group2)'];
           Cols = [handles.pws(pix) num2cell(p_nom)' num2cell(p_adj)' num2cell(ones(length(p_nom),1))];
           handles.group = group;
           handles.labels = {group1label;group2label};
           handles.data = data;
        else
            set(handles.errorText, 'Visible', 'on', 'String', 'Choose classes');
        end
        %correlation test
        case handles.corrRadio
            if ~isempty(handles.SampleData.Class)
                ClassSelected = get(handles.corr_classListbox, 'Value');
                IX =  get(handles.group1Listbox, 'UserData');
                SampleSelected = ismember(IX, ClassSelected);
            else
                SampleSelected = true(length(handles.SampleData.Value),1);
            end

            s_ix = SampleSelected & ~isnan(handles.SampleData.Value);
            data = [handles.DataStruct(2).Data(1).Data;handles.DataStruct(2).Data(2).Data];
            data = data(:,s_ix);
            vals = handles.SampleData.Value(s_ix);
            
            for ii = 1:length(handles.pws)
                try
                    d1 = data(ii,:)';
                    [rho(ii) p(ii)] = corr(vals, d1);
                catch
                    rho(ii)=NaN;
                    p(ii) = NaN;
                end
            end
            
            len = length(p);
            p_adj(1:len/2) = p(1:len/2)*sum(~isnan(p(1:len/2)));
            p_adj((len/2)+1:len) = p((len/2)+1:len)*sum(~isnan(p((len/2)+1:len)));
            p_adj(find(p_adj>=1)) = 1;
            
            % [psort pix] = sort(p, 'ascend');
            [p_nom pix] = sort(p, 'ascend');
            % add adj p value (bonferroni correction basis 2019-05-09)
            % p_adj = p_nom*(length(p_nom)/2);
            % p_adj(find(p_adj>=1)) = 1;
            p_adj = p_adj(pix);
            
            set(handles.accuracyText, 'String', [strcat('pval: ',num2str(p_nom', '%1.3e')) strcat(' rho: ', num2str(rho(pix)', '%1.3f'))]);
            RowHeaders = [[{'SampleNames'} handles.SampleData.IDs(s_ix)']; [{'Value'} num2cell(handles.SampleData.Value(s_ix))']];
            Cols = [handles.pws(pix) num2cell(p_nom') num2cell(p_adj') num2cell(rho(pix))' num2cell(ones(length(p_nom),1)*2)];
            handles.data = data;
            handles.vals = vals;
        % suvival test
        case handles.survivalRadio
            if ~isempty(handles.SampleData.Class)
               ClassSelected = get(handles.corr_classListbox, 'Value');
               IX =  get(handles.group1Listbox, 'UserData');
               SampleSelected = ismember(IX, ClassSelected);
           else
               SampleSelected = true(length(handles.SampleData.Survival),1);
           end
           survival = handles.SampleData.Survival;
           censor = handles.SampleData.Censor;
           min_grp = get(handles.minsamplesPopup,'Value')+1;

           invalid = isnan(survival)|isnan(censor);
           survival = survival(SampleSelected & ~invalid);
           censor = censor(SampleSelected & ~invalid);
           metrics = handles.metrics(:,SampleSelected & ~invalid);
           clear handles.group;

           metrics(any(isnan(metrics),2),:) = []; % exclude missing data (avoid warning message; 2016-04-17 updated)
           % metrics = (metrics+1)/0.1; % shift scores to avoid singluar matrics calculation (2016-04-17 updated)
           % re-write function code to survival analysis (2019_05_13)
           for ii = 1:size(metrics,1)
                try
                    if get(handles.randomCheck, 'Value') == 0
                        rng('default'); % for reproductivity
                    end;
                    
                    scores = metrics(ii,:)';
                    IDX = kmeans(metrics(ii,:)',2);
                    % uno = IDX == 1;
                    % dos = IDX == 2;
                    if (length(scores(IDX==1))>=5 && length(scores(IDX==2))>=5) % only both groups size greater than 5
                        if mean(scores(IDX==1))>=mean(scores(IDX==2))
                            unodata = cat(2, survival(IDX==1),censor(IDX==1));
                            dosdata = cat(2, survival(IDX==2),censor(IDX==2));
                        else
                            unodata = cat(2, survival(IDX==2),censor(IDX==2));
                            dosdata = cat(2, survival(IDX==1),censor(IDX==1));
                        end
                    
                        p1(ii) = logrank1(unodata,dosdata,'plotflag',0); % logrank test 2019_05_13

                        % original matlab code for survival test
                        % figure
                        % [cummuSurviva N] = kaplanmeier({survival(uno); survival(dos)}, {censor(uno); censor(dos)}, .5, 0);
                        % close
                        % p(ii) = ranksum(cummuSurviva{1}', cummuSurviva{2}'); % originally used Wilcoxon rank test

                        % Using cox proportional hazard test (code from statistics toolbox) (20156-04-17 updated)
                        % [urow ~] = size(survival(IDX==1));
                        % unodata = cat(2, survival(IDX==1),ones(urow,1),censor(IDX==1));
                        % [drow ~] = size(survival(IDX==2));
                        % dosdata = cat(2, survival(IDX==2),zeros(drow,1),censor(IDX==2));
                        % unodos = cat(1,unodata,dosdata);

                        [~,~,~,stats] = coxphfit(scores,survival,'censoring',censor); % Cox proportional test 2019_05_13
                        p2(ii) = stats.p;
                        hr(ii) = exp(stats.beta);
                        ci1(ii) = exp(stats.beta-1.96*stats.se);
                        ci2(ii) = exp(stats.beta+1.96*stats.se);

                        handles.group{ii}=IDX;
                    else
                    % if sum(uno)<min_grp || sum(dos)<min_grp
                        p1(ii) = NaN;
                        p2(ii) = NaN;
                        hr(ii) = NaN;
                        ci1(ii) = NaN;
                        ci2(ii) = NaN;
                        handles.group{ii} = NaN;
                    end
                    % end
                catch
                    p1(ii) = 1;
                    p2(ii) = 1;
                    hr(ii) = NaN;
                    ci1(ii) = NaN;
                    ci2(ii) = NaN;
                    handles.group{ii} = NaN;
                end
           end

          % [psort pix] = sort(p, 'ascend');
          len = length(p1);
          p_adj1(1:len/2) = p1(1:len/2)*sum(~isnan(p1(1:len/2)));
          p_adj1((len/2)+1:len) = p1((len/2)+1:len)*sum(~isnan(p1((len/2)+1:len)));
          p_adj1(find(p_adj1>=1)) = 1;
          
          p_adj2(1:len/2) = p2(1:len/2)*sum(~isnan(p2(1:len/2)));
          p_adj2((len/2)+1:len) = p2((len/2)+1:len)*sum(~isnan(p2((len/2)+1:len)));
          p_adj2(find(p_adj2>=1)) = 1;
          
          % p_adj1 = p1*(length(p1)/2);
          % p_adj1(find(p_adj1>=1)) = 1;
          % p_adj2 = p1*(length(p1)/2);
          % p_adj2(find(p_adj2>=1)) = 1;
          
          p = [p1,p_adj1,p2,p_adj2,hr,ci1,ci2];
          % [p_nom pix] = sort(p1, 'ascend');
          [p_nom,pix] = sortrows(p,1,'ascend');
          % add adj p value (bonferroni correction basis 2019-05-09)
          % p_adj = p_nom*(length(p_nom)/2);
          set(handles.accuracyText, 'String', num2cell(p_nom(:,1)));
          handles.data = metrics; % handles.data = metrics; need to convert metrics back to the original calculated scores (2016-11-16 updated)
          % handles.data = (metrics*0.1)-1;
          handles.survival=survival;
          handles.censor=censor;

          RowHeaders = [[{'SampleNames'} handles.SampleData.IDs(SampleSelected & ~invalid)']; [{'Survival'} num2cell(survival')]; [{'Censor'} num2cell(censor')]];
          Cols = [handles.pws(pix) num2cell(p_nom(:,1)') num2cell(p_nom(:,2)') num2cell(p_nom(:,3)') num2cell(p_nom(:,4)') num2cell(p_nom(:,5)') num2cell(p_nom(:,6)') num2cell(p_nom(:,7)') num2cell(ones(length(p_nom),1)*3)];
        % enrichment test
        case handles.hitsRadio
            if ~isempty(handles.SampleData.Class)
                ClassSelected = get(handles.corr_classListbox, 'Value');
                IX =  get(handles.group1Listbox, 'UserData');
                SampleSelected = ismember(IX, ClassSelected);
            else
                SampleSelected = true(length(handles.DataStruct(1).Cols),1);
            end
            %[PwNames PwGenes] = get_pwgenesLL(handles.PwData.DB, handles.PwData.Complexes, handles.PwData.UP2LL);
            load PwNames;
            load PwGenes;
            [A D p fisheromnibusp] = HitTest([handles.DataStruct(1).Data(3).Data{1}{:}], handles.DataStruct(1).Data(3).Data{2}(:,SampleSelected), PwNames, PwGenes);
            % [psort pix] = sort(fisheromnibusp, 'ascend');
            
            len = length(p);
            p_adj(1:len/2) = p(1:len/2)*sum(~isnan(p(1:len/2)));
            p_adj((len/2)+1:len) = p((len/2)+1:len)*sum(~isnan(p((len/2)+1:len)));
            p_adj(find(p_adj>=1)) = 1;
            
            [p_nom pix] = sort(fisheromnibusp, 'ascend');
            % add adj p value (bonferroni correction basis 2019-05-09)
            % p_adj = p_nom*(length(p_nom)/2);
            p_adj = p_adj(pix);
            set(handles.accuracyText, 'String', num2cell(p_nom));
            RowHeaders = [{'SamplesIncluded'} handles.DataStruct(1).Cols(SampleSelected)];
            Cols = [PwNames(pix) num2cell(p_nom') num2cell(p_adj)' num2cell(ones(length(p_nom),1))];
    end
    
    %{
for ii = 1:length(pix)
        try
            top_pws{ii} = PwNames{pix(ii)}(1:70);
        catch
            top_pws{ii} = PwNames{pix(ii)};
        end
    end
%}
    if get(handles.test_typePanel, 'SelectedObject') ~=handles.hitsRadio
    set(handles.plotButton, 'Enable', 'on');
    end
    set(handles.heatmapButton, 'Enable', 'on', 'UserData', pix);
    set(handles.writeButton, 'Enable', 'on', 'UserData', {RowHeaders;Cols});
    set(handles.top_pwsText, 'String', Cols(:,1), 'UserData', {pix;get(handles.test_typePanel, 'SelectedObject')});
    set(handles.goButton, 'String', 'GO');
catch
    set(handles.errorText, 'Visible', 'on', 'String', 'Error')
    set(handles.goButton, 'String', 'GO');
end
  guidata(hObject, handles);
  
function sampleclassPanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in sampleclassPanel 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if hObject == handles.listRadio
    set(handles.samplesListbox, 'Visible', 'on')
    set(handles.indicesText, 'Visible', 'off');
    set(handles.indicesEdit, 'Visible', 'off');
else
    set(handles.indicesText, 'Visible', 'on');
    set(handles.indicesEdit, 'Visible', 'on');
    set(handles.samplesListbox, 'Visible', 'off')
end

function returnButton_Callback(hObject, eventdata, handles)
% hObject    handle to returnButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume

function writeButton_Callback(hObject, eventdata, handles)
% hObject    handle to writeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   u = get(hObject, 'UserData');
   [RowHeaders Cols] = deal(u{:});
   write_data(RowHeaders, Cols);
   set(handles.writeButton, 'BackgroundColor', 'yellow')

function setgroupButton_Callback(hObject, eventdata, handles)
% hObject    handle to setgroupButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Classes = get(handles.group1Listbox, 'String');
ClassSelected = zeros(length(handles.Classes),1);
ClassSelected(get(handles.group1Listbox, 'Value'))=1;
IX =  get(handles.group1Listbox, 'UserData');
handles.group1 = ismember(IX, find(ClassSelected));
set(handles.group2Listbox, 'Enable', 'on', 'String', handles.Classes(~ClassSelected), 'UserData', find(~ClassSelected), 'Value', 1);
guidata(hObject, handles);

function heatmapButton_Callback(hObject, eventdata, handles)
% hObject    handle to heatmapButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
num_pws =get(handles.percentText, 'String');
try 
    num_pws = str2double(num_pws);
    top_pws = get(hObject, 'UserData');
    for ii = 1:num_pws
        if top_pws(ii)>length(handles.pws)/2;
            top_pws(ii) = top_pws(ii) - length(handles.pws)/2;
        end
    end
    handles.DataStruct(2).RowSelect = unique(top_pws(1:num_pws));
    set(hObject, 'BackgroundColor', 'yellow');
    guidata(hObject, handles);
catch
    set(handles.errorText, 'Visible', 'on');
end
  

function test_typePanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in test_typePanel 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.writeButton, 'BackgroundColor', handles.bgcolor, 'Enable', 'off');
set(handles.top_pwsText, 'String', '');
set(handles.accuracyText, 'String', '');
set(handles.heatmapButton, 'BackgroundColor', handles.bgcolor, 'Enable', 'off');

function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.next = 1;
guidata(hObject, handles);
uiresume

function plotButton_Callback(hObject, eventdata, handles)
% hObject    handle to plotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
u = get(handles.top_pwsText, 'UserData');
[pix type] = deal(u{:});
pw = get(handles.top_pwsText, 'Value');

% figure
if type == handles.binaryRadio
    figure
    gscatter(1:length(handles.group), handles.data(pix(pw),:), handles.group)
    xlabel('samples')
    title(handles.pws{pix(pw)});
    legend(handles.labels, 'Location', 'Best');
elseif type == handles.corrRadio
    figure
    scatter(handles.vals, handles.data(pix(pw),:), '.');
    hold on
    coeffs = polyfit(handles.vals, handles.data(pix(pw),:)',1);
    x = min(handles.vals):range(handles.vals)/100:max(handles.vals);
    y = polyval(coeffs,x);
    plot(x, y)
    title(handles.pws{pix(pw)});
    hold off
elseif type == handles.survivalRadio
    labels1 = {'high';'low'};
    labels2 = {'low';'high'};
    % be careful for the graph object setting, matlab version changes, graph position might change (2016-04-17 updated)
    f= figure('Position', [150 150 1220 500], 'Tag', 'survival');
    leftPanel = uipanel('bordertype','etchedin', 'Units','characters','Position',[1 1 100 50],'Parent',f);
    middlePanel = uipanel('bordertype','etchedin', 'Units','characters','Position',[100 1 100 50],'Parent',f);
    rightPanel = uipanel('bordertype','etchedin', 'Units','characters','Position',[145 1 100 50],'Parent',f);
    leftAxes = axes('Parent', leftPanel, 'Position', [.1 .07 .85 .65]); % 2016-04-17 updated
    middleAxes = axes('Parent', middlePanel, 'Position', [.1 .07 .85 .65]);
    rightAxes = axes('Parent', rightPanel, 'Position', [.1 .07 .85 .65]); % 2016-04-17 updated
    
    IDX = handles.group{pix(pw)};
    uno = IDX==1;
    uno_avg = sum(handles.data(pix(pw),uno))/sum(uno);
    dos = IDX==2;
    dos_avg = sum(handles.data(pix(pw),dos))/sum(dos);
     if uno_avg>dos_avg
        labels = labels1;
    else
        labels = labels2;
     end
     
    axes(leftAxes)
    gscatter(1:length(IDX), handles.data(pix(pw),:), IDX);
    xlabel('samples');
    if ~isempty(strfind(handles.pws{pix(pw)}, 'act:'))
       ylabel('activity');
    else
       ylabel('consistency');
    end
    title(handles.pws{pix(pw)}(5:end));
    legend(labels);
    
    % additional table showing the scores w/ assigned group ids by kmean
    % function (row numbers are the same order along w/ x scale on the
    % scatter plot) (2016-11-16 updated)
    % may update the code to associate the scatter plot in near future
    axes(middleAxes)
    if ~isempty(strfind(handles.pws{pix(pw)}, 'act:'))
       uitable('Data', [handles.data(pix(pw),:);IDX']', 'ColumnName',{'activity','group id'},'Position',[505 1 219 500]);
    else
       uitable('Data', [handles.data(pix(pw),:);IDX']', 'ColumnName',{'consistency','group id'},'Position',[505 1 219 500]);
    end
    
    % original code for KM plot (not used; 2016-04-17 updated)
    % [cummuSurviva N] = kaplanmeier({handles.survival(uno); handles.survival(dos)}, {handles.censor(uno); handles.censor(dos)}, .5, rightAxes); % orignally used for KM plot
    
    %Using Kaplan-Meier code (code from statistics toolbox; 2016-04-17 updated)
    axes(rightAxes)
%     unodata = cat(2, handles.survival(IDX==1), handles.censor(find(IDX==1)));
%     [uf ux] = ecdf(unodata(:,1),'censoring',unodata(:,2),'function','survivor');
%     stairs(ux, uf);
%     hold on
%     dosdata = cat(2, handles.survival(IDX==2), handles.censor(find(IDX==2)));
%     [df dx] = ecdf(dosdata(:,1),'censoring',dosdata(:,2),'function','survivor')
%     stairs(dx, df, 'color','r');
%     legend(labels);
%     xlabel('time');
%     ylabel('percent surviving');
%     title('Kaplan Meier Survival Plot');
    if mean(scores(IDX==1))>=mean(scores(IDX==2))
        unodata = cat(2, survival(IDX==1),censor(IDX==1));
        dosdata = cat(2, survival(IDX==2),censor(IDX==2));
    else
        unodata = cat(2, survival(IDX==2),censor(IDX==2));
        dosdata = cat(2, survival(IDX==1),censor(IDX==1));
    end
    
    p = logrank1(unodata,dosdata,'plotflag',1); % logrank test 2019_05_13
    hold off;
end

function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
if get(hObject, 'Value')
    set(handles.corr_classListbox, 'Value', 1:length(get(handles.corr_classListbox, 'String')));
else
    set(handles.corr_classListbox, 'Value',1);
end

% Blank Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function group1Listbox_Callback(hObject, eventdata, handles)
% hObject    handle to group1Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function group1Listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to group1Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function group2Listbox_Callback(hObject, eventdata, handles)
% hObject    handle to group2Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function group2Listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to group2Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function minsamplesPopup_Callback(hObject, eventdata, handles)
% hObject    handle to minsamplesPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function minsamplesPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minsamplesPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function corr_classListbox_Callback(hObject, eventdata, handles)
% hObject    handle to corr_classListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function top_pwsText_Callback(hObject, eventdata, handles)
% hObject    handle to top_pwsText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns top_pwsText contents as cell array
%        contents{get(hObject,'Value')} returns selected item from top_pwsText

set(handles.accuracyText, 'Value', get(hObject, 'Value'));
function percentText_Callback(hObject, eventdata, handles)
% hObject    handle to percentText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function percentText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to percentText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function corr_classListbox_CreateFcn(hObject, ~, handles)
% hObject    handle to corr_classListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function indicesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to indicesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function indicesEdit_Callback(hObject, eventdata, handles)
% hObject    handle to indicesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function samplesListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to samplesListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function samplesListbox_Callback(hObject, eventdata, handles)
% hObject    handle to samplesListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Helper Functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PwNames PwGenes] = get_pwgenesLL(db, complexes, UP2LL)

[geneIDs ia ib]  = unique([db{:,5}]);
unparsedLLs = db(ia,6);
for ii = 1:length(unparsedLLs)
db(ib==ii,6) = {get_parsedlink(unparsedLLs{ii}, UP2LL)};
end

[geneIDs ia ib]  = unique([complexes{:,4}]);
unparsedLLs = complexes(ia,5);
for ii = 1:length(unparsedLLs)
complexes(ib==ii,5) = {get_parsedlink(unparsedLLs{ii}, UP2LL)};
end

[pws ia ib] = unique([db{:,2}]);
PwNames = db(ia,1);

PwGenes = cell(length(pws),1);
for ii = 1:length(pws);
pw1 = db(ib==ii,:);
LL = vertcat(pw1{:,6});
if isempty(LL)
    LL={[]};
end
c_id = unique([pw1{strmatch('complex',pw1(:,3)), 5}]);
c_LL = vertcat(complexes{ismember([complexes{:,1}], c_id),5});
if isempty(c_LL)
    c_LL = {[]};
end

PwGenes{ii} = unique([LL{:} c_LL{:}]);

end
function parsedlink = get_parsedlink(str, UP2LL)
    %turn linking string for molecule into useful search term in linking table
        L = findstr('LL:', str);
        U = findstr('UP:', str);
        C = findstr(',', str);

        if L
            for i = 1:length(L)
                start = L(i) + 3;
                if length(C)>=i
                    parsedlink{i,1} = str2num((str(start:C(i) - 1)));
                else
                    parsedlink{i,1} = str2num((str(start:end))); 
                end
            end
            
        elseif U
            for i = 1:length(U)
                start = U(i) + 3;
                if length(C)>=i
                    pl = str(start:C(i) - 1);
                else
                    pl = str(start:end);
                end
                try
                parsedlink{i,1} = UP2LL{strmatch(pl, UP2LL(:,1)), 2};
                catch
                    parsedlink{i,1} = [];
                end
            end
            
        else 
            parsedlink = [];
        end
function [A D p fisheromnibusp] = HitTest(GeneHitsLabels, HitsMatrix, PwNames, PwGenes)
HitsMatrix = abs(HitsMatrix);
B = length(GeneHitsLabels);
%hits{1} = abs(ld-nd)>1;
%hits{2} = abs(hd-nd)>1;
%for h = 1:2
C = sum(HitsMatrix);
for ii = 1:length(PwNames)
    [GeneMatch ia ib] = intersect(PwGenes{ii}, GeneHitsLabels);
    D(ii) = length(GeneMatch);
    %if D(ii) ==1
        %pname = pws{ii};
    %end
    PwHits = HitsMatrix(ib,:);
    A(ii,:) = sum(PwHits,1);

    for jj = 1:size(A,2)
        ptmp = hygecdf(A(ii,jj)-1, B, C(jj), D(ii));
        
        p(ii,jj) = 1 - ptmp;  
    end
fisheromnibusp(ii) = testfisheromnibus(p(ii,~isnan(p(ii,:))));
end
function pval = testfisheromnibus(val)

% TESTFISHEROMNIBUS computes p-value for the Fisher Omnibus test.
% -----------------------------
% pval = testfisheromnibus(val)
% -----------------------------
% Description: computes p-value for the Fisher Omnibus test.
% Input:       {val} list of p-values derived from a series of one-sided
%                   tests.
% Output:      {pval} p-value of Fisher Omnibus.

% © Liran Carmel
% Classification: Hypothesis testing
% Last revision date: 30-Apr-2006

% replace zeros by small numbers
val(val==0) = eps;

% test statistics
stat = -2*sum(log(val));

% p-value
pval = 1 - chi2cdf(stat,2*length(val));

function write_data(HeaderRows, Cols)
numCols = size(Cols,2);
[file_name,file_dir] = uiputfile('*.txt','Save pathway p-vals as: ');
Directory= cd;
cd(file_dir);
fid = fopen(file_name, 'wt');
for hr = 1:size(HeaderRows,1)
fprintf(fid, '%s', HeaderRows{hr,1});
if ischar(HeaderRows{hr,2})
    fprintf(fid, '\t%s', HeaderRows{hr,2:end});
else
    fprintf(fid, '\t%f', HeaderRows{hr,2:end});
end
fprintf(fid, '\n');
end
fprintf(fid, '%s', 'Pathways');
if Cols{1,4} == 1
    fprintf(fid, '\t%s', 'P-Vals');
    fprintf(fid, '\t%s', 'adjP-Vals');
elseif Cols{1,5} == 2
    fprintf(fid, '\t%s', 'P-Vals');
    fprintf(fid, '\t%s', 'adjP-Vals');
    fprintf(fid, '\t%s', 'Rho');
elseif Cols{1,9} == 3
    fprintf(fid, '\t%s', 'logrank P-Vals');
    fprintf(fid, '\t%s', 'logrank adjP-Vals');
    fprintf(fid, '\t%s', 'Cox prop P-Vals');
    fprintf(fid, '\t%s', 'Cox prop adjP-Vals');
    fprintf(fid, '\t%s', 'Hazard Ratio');
    fprintf(fid, '\t%s', 'upper 95% CI');
    fprintf(fid, '\t%s', 'lower 95% CI');
end

for cc = 1:size(Cols,1)
    fprintf(fid,'\n%s', Cols{cc,1});
    if numCols == 4
        fprintf(fid, '\t%1.3e', Cols{cc,2});
        fprintf(fid, '\t%1.4e', Cols{cc,3});
    elseif numCols == 5
        fprintf(fid, '\t%1.3e', Cols{cc,2});
        fprintf(fid, '\t%1.4e', Cols{cc,3});
        fprintf(fid, '\t%1.4f', Cols{cc,4});
    elseif numCols == 8
        fprintf(fid, '\t%1.3e', Cols{cc,3});
        fprintf(fid, '\t%1.4e', Cols{cc,4});
        fprintf(fid, '\t%1.3e', Cols{cc,5});
        fprintf(fid, '\t%1.4e', Cols{cc,6});
        fprintf(fid, '\t%1.4f', Cols{cc,7});
        fprintf(fid, '\t%1.4f', Cols{cc,7});
        fprintf(fid, '\t%1.4f', Cols{cc,7});
    end
end
fclose(fid);
cd(Directory);

% original kaplan-meiyer code (not used; 2016-04-17 updated)
% function [cummuSurviva N] = kaplanmeier(datacollection, censordatacollection, binsize, drawkm)
% 
% colorchoser = {'r', 'b', 'g', 'c', 'm', 'y', 'k', [0.5 0.5 0.5], [0.1 0.5 0.9]};
% span=0;
% for ii=1:length(datacollection)
%     if max(datacollection{ii}) > span
%         span = max(datacollection{ii});
%     end
% end
%     
% strlegend = {};
% for idx = 1:length(datacollection)
%     data = datacollection{idx};
%     strlegend(length(strlegend)+1) = {num2str(length(data))};
%     span = max(data);
%     N(idx) = length(data);
%     censordata = censordatacollection{idx};
%     cummuSurvivalAtEndOfInterval=[];
% 
%     rmvdata = find(data==-999);
%     data(rmvdata) = [];
%     censordata(rmvdata) = [];
% 
%     intervals = 0:binsize:span;
% 
%     riskAtIntervalStart = length(find(data>0));
%     censoredDuringInterval = length(find(censordata(find(data<=binsize))==0));
%     atRiskAtEndOfInterval = riskAtIntervalStart - censoredDuringInterval;
%     diedAtEndOfInterval = length(find(censordata(find(data<=binsize))==1));
%     proportionSurvivingThisInterval = (atRiskAtEndOfInterval - diedAtEndOfInterval)/ atRiskAtEndOfInterval ;
%     cummuSurvivalAtEndOfInterval(1) = proportionSurvivingThisInterval;
% 
% 
%     for ii=2:length(intervals)
%         riskAtIntervalStart = length(find(data>binsize*(ii-1)));
%         censoredDuringInterval = length(find(censordata(find((data<=binsize*ii) & (data>binsize*(ii-1))))==0));
%         atRiskAtEndOfInterval = riskAtIntervalStart - censoredDuringInterval;
%         diedAtEndOfInterval = length(find(censordata(find((data<=binsize*ii) & (data>binsize*(ii-1))))==1));
%         if atRiskAtEndOfInterval>0
%             proportionSurvivingThisInterval = (atRiskAtEndOfInterval - diedAtEndOfInterval)/ atRiskAtEndOfInterval ;
%         else
%             proportionSurvivingThisInterval = 1;
%         end
%         cummuSurvivalAtEndOfInterval(ii) = proportionSurvivingThisInterval*cummuSurvivalAtEndOfInterval(ii-1);
%     end
% % if drawkm
% if ~isempty(findobj('Tag', 'survival'))  % verify the graph obejct tag (survival plot)
%     try
%         axes(drawkm)
%     end
%     stairs(intervals, cummuSurvivalAtEndOfInterval, '-', 'Color', colorchoser{idx}); hold on;
% end
%     cummuSurviva{idx} = cummuSurvivalAtEndOfInterval;
% end
% 
% % legend(strlegend);
% 
% hold off;

function randomCheck_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject, 'Value')
    set(handles.randomCheck, 'Value', 1);
else
    set(handles.randomCheck, 'Value', 0);
end;
