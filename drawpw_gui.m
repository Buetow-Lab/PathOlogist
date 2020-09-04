function varargout = drawpw_gui(varargin)
% DRAWPW_GUI M-file for drawpw_gui.fig
%      DRAWPW_GUI, by itself, creates a new DRAWPW_GUI or raises the existing
%      singleton*.
%
%      H = DRAWPW_GUI returns the handle to a new DRAWPW_GUI or the handle to
%      the existing singleton*.
%
%      DRAWPW_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DRAWPW_GUI.M with the given input arguments.
%
%      DRAWPW_GUI('Property','Value',...) creates a new DRAWPW_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before drawpw_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to drawpw_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help drawpw_gui

% Last Modified by GUIDE v2.5 01-Mar-2016 21:51:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @drawpw_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @drawpw_gui_OutputFcn, ...
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

% --- Executes just before drawpw_gui is made visible.
function drawpw_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to drawpw_gui (see VARARGIN)
warning('off', 'all');
% Choose default command line output for drawpw_gui
handles.output = hObject;
u = get(hObject, 'UserData');
[version handles.PwData handles.DataStruct handles.ChipData] = deal(u{:});
handles.drawn = 0;
if version ~=1
    set(handles.bgCheck, 'Value', 0, 'Enable', 'off');
end
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
    if handles.DataStruct(1).Data(2).Loaded
        set(handles.acRadio, 'Enable', 'on', 'Value', 1);
    end   
    if handles.DataStruct(1).Data(3).Loaded
        set(handles.cnaRadio, 'Enable', 'on', 'Value', 1);
    end


guidata(hObject, handles);
uiwait

end

function varargout = drawpw_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.drawn;
end

% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_pwButton_Callback(hObject, eventdata, handles)
% hObject    handle to draw_pwButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.one_pwText, 'Visible', 'off');
set(handles.errorText,'Visible','off');
set(handles.finishedText, 'Visible', 'off');

handles.DataStruct(1).RowSelect= get(handles.pw_Listbox, 'Value');
handles.DataStruct(1).ColSelect = get(handles.samplesListbox, 'Value');

try
    set(handles.draw_pwButton, 'String', 'Drawing...', 'BackgroundColor', 'white');
    drawnow expose
    
    if length(handles.DataStruct(1).RowSelect) ~= 1;
        set(handles.one_pwText, 'Visible', 'on');
    else
        pwname = handles.PwData.Names{handles.DataStruct(1).RowSelect};
        pwnum = handles.PwData.Nums(handles.DataStruct(1).RowSelect);
        

%         if get(handles.PIDRadio, 'Value')|| ~(handles.DataStruct(1).Data(2).Loaded || handles.DataStruct(1).Data(3).Loaded)
%             %if no data has been loaded, link to online pathway drawing
%             pwname_trim = pwname(1:strfind(pwname, '(')-1);
%             url = ['http://pid.nci.nih.gov/pathway_names?what=graphic&svg=&gif=true&xml=&biopax=&format=html&word=', pwname_trim, '&Submit=Go'];
%             web(url);
%             set(handles.draw_pwButton, 'String', 'Draw Pathway');
%         else
            tc = get(handles.textCheck, 'Value');
            bc = get(handles.bgCheck, 'Value');
            draw = get(handles.acRadio, 'Value') + 2*get(handles.cnaRadio, 'Value');
            %draw = 0:neither, 1:ac, 2:cna, 3:both
            udp_s = handles.DataStruct(1).ColSelect;
            
            if get(handles.addlevelCheck, 'Value')
                levels=2;
            else
                levels=1;
            end
            [molList molList_map A] = get_molList(pwnum, levels, udp_s, handles);       
 
            if mod(draw,2)==1
                [i_activity i_consistency omit_i mult_out] = pw_calc(A, molList, molList_map, udp_s);
            else 
                i_activity=[]; i_consistency=[]; mult_out=[];
            end
            
            paren = strfind(pwname, '(');
            if paren(end)>35
               pwname = [pwname(1:35) pwname(paren(end):end)];
            end
            pwname(isstrprop(pwname,'punct')) = '_';
            if (levels>1)
                pwname = [pwname '_plus'];
            end
           
             c_p = molList.compProbs;
            [i_list bgwrite] = pw_draw();
            
            if tc
                Directory = cd;
                dir_name = uigetdir(Directory,'Select an output directory');
                mkdir(strcat(dir_name, '\',pwname));
                
                sample_labels = handles.DataStruct(1).Cols(handles.DataStruct(1).ColSelect);              
                
                comp_write = [{'complexes'} num2cell(NaN*ones(1,length(udp_s)))];        
                for ii = 1:size(c_p,1)
                    if ~isempty(c_p{ii})
                        comp_write = [comp_write; strcat('*',molList.molName{ii}) num2cell(NaN*ones(1,length(udp_s))); c_p{ii}];
                    end
                end

                cd(dir_name);
                cd(pwname);
            %{
                for ss = 1:length(sample_labels)
                    fid = fopen(['PWstructure_' sample_labels{ss} '.txt'], 'wt');
                    for ii = 1:length(bgwrite{ss})
                        fprintf(fid, bgwrite{ss}{ii});
                    end
                end
                   fclose(fid);
              %}  
                
                if mod(draw,2) ==1
                    mol_write = [[{'molecules'} sample_labels]; [molList.molName num2cell(molList.molProbs)]];
                    if size(comp_write,1)>1
                        mol_write = [mol_write; num2cell(NaN*ones(1,length(udp_s)+1)); comp_write];
                    end
                    act_write = [[{'intxn num' 'inputs' 'outputs'} sample_labels]; [i_list' num2cell(i_activity)]];
                    con_write = [[{'intxn num' 'inputs' 'outputs'} sample_labels]; [i_list' num2cell(i_consistency)]];
                    
                    fid = fopen(strcat(pwname,'_molprobs.txt'), 'wt');
                    fprintf(fid, '%s\t', mol_write{1,1:end-1});
                    fprintf(fid, '%s', mol_write{1,end});
                    for rr = 2:size(mol_write,1)
                        fprintf(fid,'\n%s', mol_write{rr,1});
                        fprintf(fid, '\t%1.4f', mol_write{rr,2:end});
                    end
                    fclose(fid);
                
                    fid = fopen(strcat(pwname,'_act.txt'), 'wt');
                    fprintf(fid, '%s\t', act_write{1,1:end-1});
                    fprintf(fid, '%s', act_write{1,end});
                    for rr = 2:size(act_write,1)
                        fprintf(fid,'\n%s\t%s\t%s', act_write{rr,1:3});
                        fprintf(fid, '\t%1.4f', act_write{rr,4:end});
                    end
                    fclose(fid);
            
                    fid = fopen(strcat(pwname,'_con.txt'), 'wt');
                    fprintf(fid, '%s\t', con_write{1,1:end-1});
                    fprintf(fid, '%s\t', con_write{1,end});
                    for rr = 2:size(con_write,1)
                        fprintf(fid,'\n%s\t%s\t%s', con_write{rr,1:3});
                        fprintf(fid, '\t%1.4f', con_write{rr,4:end});
                    end
                    fclose(fid);
                end
           
                if draw >1
                    cna_write = [[{'molecules'} sample_labels]; [molList.molName num2cell(molList.cnaNum)]];
                    fid = fopen(strcat(pwname,'_cna.txt'), 'wt');
                    fprintf(fid, '%s\t', cna_write{1,:});
                    for rr = 2:size(cna_write,1)
                        fprintf(fid,'\n%s', cna_write{rr,1});
                        fprintf(fid, '\t%1.0f', cna_write{rr,2:end});
                    end
                    fclose(fid);
                end
                cd(Directory);
                set(handles.finishedText, 'Visible', 'on', 'String', ['Pathway data has been written to the folder ''' pwname ''' in the chosen directory.']);
            end
                handles.drawn = 1;
                guidata(hObject, handles);
%             end
        end
 catch
    set(handles.errorText, 'Visible','on');
end
set(handles.draw_pwButton, 'String', 'Draw Pathway');

function [i_list bgwrite] = pw_draw()

%make interaction node labels based on type of intxn
[intxn_ids ia intxn_map]= unique(A(:,11));
number_of_intxns = length(intxn_ids);

for kk = 1:number_of_intxns
    type = A(intxn_map==kk,12);
    intxnLabels{kk,1} = strcat(num2str(kk), ':', type{1});
end

%make molecule node labels, with active/inactive form of molecule a separate node from regular form
adj_mol_ids = A(:,4);
for i = 1:length(adj_mol_ids)
    if strcmp(A{i,7}, 'active')
        adj_mol_ids{i} = strcat(adj_mol_ids{i},'+');
    elseif strcmp(A{i,7}, 'inactive')
        adj_mol_ids{i} = strcat(adj_mol_ids{i},'-');
    else
        adj_mol_ids{i} = adj_mol_ids{i};
    end          
end
[molLabels ia molLabels_map]= unique(adj_mol_ids);
number_of_mols = length(molLabels);
molmap(molLabels_map) = molList_map;

%set up space to store info
agent_edges = [];
inhibitor_edges = [];
loc_matrix = [];

%create empty biograph data matrix
node_names = [molLabels; intxnLabels];
number_of_nodes = length(node_names);
bg_data = zeros(number_of_nodes, number_of_nodes);
edge_color = cell(number_of_nodes, number_of_nodes);
edge_width = zeros(number_of_nodes, number_of_nodes);
con_nodes = zeros(number_of_intxns,2);

%fill in biograph data matrix and set edge colors
for i = 1:number_of_intxns; 
    input_string = [];
    output_string = [];
	intxn_group = intxn_map == i;
    mols = molLabels_map(intxn_group); 
    roles = A(intxn_group, 10);
    locs = A(intxn_group, 9);
    intxn_ind = number_of_mols +i;
    m_names = molLabels(mols);
  
    for j = 1:length(mols);
        if strcmp(roles{j}, 'output'); 
            x = intxn_ind;
            y = mols(j);
            con_nodes(i,1:2) = [x y]; 
            output_string = [output_string, m_names{j}, ', '];
         
        else
            x = mols(j);
            y = intxn_ind;
            input_string = [input_string, m_names{j}, ', '];
            if strcmp(roles{j}, 'agent');
                agent_edges = [agent_edges; i mols(j)]; %store location of agent edges
                edge_color{x,y} = [.5 .75 0];
            elseif strcmp(roles{j}, 'inhibitor');
                inhibitor_edges = [inhibitor_edges; i mols(j)];%store location of inhibitor edges
                edge_color{x,y} = [.95 .6 .4];
            else
                edge_color{x,y} = [.5 .5 .5];
            end
        end
     
        bg_data(x,y) = 1;
        loc_matrix = [loc_matrix; {x} {y} {locs{j}}]; %store location of each step in the interaction
          
    end
i_inputs{i} = input_string;
i_outputs{i} = output_string;
end
i_list = [intxnLabels'; i_inputs; i_outputs];

%set up copy number outline colors
if draw >1
    cna_colors = {[1 .5 0]; [0 0 0]; [1 .92 0]};
    CNA = handles.DataStruct(1).Data(3).Data{2};
    cna_genes = handles.DataStruct(1).Data(3).Data{1};
    cna_genes = vertcat(cna_genes{:});
    cna_s = handles.DataStruct(1).ColSelect;
    for mm = 1:length(molList.molName)
        try
            LL = molList.molLink{mm}(:,2);
            LL = vertcat(LL{:});
            cc = mode(CNA(ismember(cna_genes,LL),cna_s),1);
            cc(isnan(cc)) = 0;
            molList.cnaList{mm} = cna_colors(cc+2)';
            %{
            if cc<8
                molList.cnaList{mm} = {[0 cc/7 0 ]};
            else
                molList.cnaList{mm} = {[(cc-7)/7 1 0]};
            end
            %}
            molList.cnaNum(mm,1:length(cna_s)) = cc;
        catch
            %molList.cnaList{mm} = {[]};
            molList.cnaNum(mm,1:length(cna_s)) = NaN;
        end
    end
end
if bc || tc
%store ids of complexes
complex_ids = [];mol_ids = A(:,5);
for k = 1:number_of_mols
    group = molLabels_map == k;
    type = A(group, 3);
    if strcmp(type{1}, 'complex') 
        is_complex = mol_ids(group);
        complex_num = is_complex{1};
        complex_ids = [complex_ids; {k} c_p(molmap(k))];
    end
end


num_samples =1;
%assign each mol a probability
if mod(draw,2)==1
mol_probs = molList.molProbs(molmap,:);
num_samples = size(mol_probs,2);
end

if draw >1
mol_cnas = molList.cnaList(molmap,:);
mol_cna = molList.cnaNum(molmap,:);
num_samples = size(mol_cna,2);
end

if mod(draw,2)==1
    i_con_plus = i_consistency;
    
    if numel(mult_out)
        
        mult_out(:,1) = round(mult_out(:,1))+number_of_mols;
        mult_out(:,2) = molLabels_map(round(mult_out(:,2)));
        con_nodes = [con_nodes; mult_out(:,1:2)];
        i_con_plus = [i_con_plus; mult_out(:,3:end)];
    end
end

for ss = 1:num_samples
    if mod(draw,2)==1
        edge_width(logical(bg_data)) = .5;
        for i = 1:size(i_con_plus,1)
            if con_nodes(i,1)~= 0
               edge_width(con_nodes(i,1),con_nodes(i,2))= 2;
                if ~isnan(i_con_plus(i,ss))
                    edge_color{con_nodes(i,1),con_nodes(i,2)}=[0 .25 i_con_plus(i,ss)];
                else    
                    edge_color{con_nodes(i,1),con_nodes(i,2)}=[1 1 0];
                end
            end
        end
    else 
        mol_probs = [];
        i_list = [];
        i_activity = [];
        i_consistency = [];
        edge_color(logical(bg_data)) = {[.75 .75 .75]};
        edge_width(logical(bg_data)) = .5;
    end


    node_color=cell(number_of_nodes,1); outline_color=cell(number_of_nodes,1); outline_width=zeros(number_of_nodes,1);shape=cell(number_of_nodes,1);node_size=cell(number_of_nodes,1);
 
    for ii = 1:number_of_nodes
        intxn_num = ii - number_of_mols;
    if intxn_num>0
        node_color{ii} = [.75 .85 1];
        outline_color{ii} = [.4 .4 .4];
        outline_width(ii) = .5;
        shape{ii} = 'circle';
        if mod(draw,2)==1 && ~isnan(i_activity(intxn_num,ss))
            node_size{ii} = [(10+20*i_activity(intxn_num,ss)) (10+20*i_activity(intxn_num,ss))];
        else
            node_size{ii} = [20 20];
        end
    else
        shape{ii} = 'box';
        node_size{ii} = [7*length(molLabels{ii}) 20];
        if mod(draw,2)==1&& ~isnan(mol_probs(ii,ss))
            node_color{ii} = [1-mol_probs(ii,ss) 1 1];
        else
            node_color{ii} = [.85 .85 .85];
        end
        if draw>1 && ~isempty(mol_cnas{ii})
            outline_color{ii} = mol_cnas{ii}{ss};
            if isnan(mol_cna(ii,ss))
                outline_width(ii) =  .5;
            %elseif mol_cna(ii,ss)==0
             %   outline_width(ii) =  .1;
            else
                outline_width(ii)=2.5;
            end
        else
            outline_color{ii} = [0 0 0];
            outline_width(ii) = .5;
        end
    end
    
 end
 sample_labels = handles.DataStruct(1).Cols(handles.DataStruct(1).ColSelect);
        
if bc
    make_BG()
end
if tc
    %bgwrite{ss} = bg2text();
    bgwrite = [];
else
    bgwrite = [];
end
end
end
function make_BG()
%make biograph
   
    BG = biograph(bg_data, node_names, 'ID',sample_labels{ss} , 'Label', strcat('Pathway ', num2str(pwnum), ': ', pwname), 'ShowWeights', 'off', 'EdgeFontSize', 7, 'LayoutType', 'hierarchical');

    %set node properties
    for ii = 1:number_of_nodes;
        set(BG.nodes(ii), 'Color', node_color{ii}, 'LineColor', outline_color{ii}, 'LineWidth', outline_width(ii), 'TextColor', [0 0 0], 'Shape', shape{ii});
    end
    for iii = 1:size(complex_ids,1)
        %set(BG.nodes(complex_ids(iii,1)), 'Description', num2str(complex_ids(iii,2)), 'UserData', c_p);
        try
        set(BG.nodes(complex_ids{iii,1}), 'Description', 'complex', 'UserData', [complex_ids{iii,2}(:,1) complex_ids{iii,2}(:,ss+1)]);
        end
    end
    

    %set edge properties
    for bgr = 1:number_of_nodes
        for bgc = 1:number_of_nodes
            if bg_data(bgr,bgc)
                edge = getedgesbynodeid(BG,BG.nodes(bgr).ID, BG.nodes(bgc).ID);
                set(edge, 'LineColor', edge_color{bgr, bgc});
                set(edge,'LineWidth', edge_width(bgr,bgc));
            end
        end
    end

    %add location labels to edges
    for kk = 1:size(loc_matrix,1);
        if ~strcmp(loc_matrix{kk,3}, 'null')
            source = BG.nodes(loc_matrix{kk,1});
            sink = BG.nodes(loc_matrix{kk,2});
            edge = getedgesbynodeid(BG, source.ID, sink.ID);
            try
            set(edge, 'Label', loc_matrix{kk,3}(1:3));
            end
        end
    end

    %set callback to open gene webpage, or show complex components
    set(BG, 'NodeCallback', {@(node)nodeclick(node, handles)});

    dolayout(BG)


    %make intxn node size separately after initial layout
    
        set(BG, 'NodeAutoSize', 'off');
        for ii = number_of_mols+1:number_of_nodes
            set(BG.nodes(ii), 'Size', node_size{ii});
            if mod(draw,2)==1
                set(BG.nodes(ii), 'Label', ['activity: ', num2str(i_activity(ii-number_of_mols,ss))]);
            end
        end
        
        for ii = 1:number_of_mols
            if mod(draw,2)==1
            try
                set(BG.nodes(ii), 'Label', ['probability: ' num2str(mol_probs(ii,ss))]);
            catch
                set(BG.nodes(ii), 'Label', 'probability unknown');
            end
            end
        end

    
    dolayout(BG)
    view(BG)
 end
function bgwrite = bg2text()
        
        bgwrite{1} = ['Pathway ' num2str(pwnum), ': ' pwname '\n' sample_labels{ss} '\n\n*Nodes \t' num2str(number_of_nodes) '\n'];
        bgwrite{2} = '*0 \t Example: molecule_name \t shape \t color \t linecolor \t linewidth \t size \n';
        op = '[';
        cl =']\t';
        for nn = 1:number_of_nodes
            bgwrite{nn+2} = [num2str(nn) '\t' node_names{nn} '\t' shape{nn} '\t' op deblank(num2str(node_color{nn},'%.2f ')) cl op deblank(num2str(outline_color{nn},'%.2f ')) cl op num2str(outline_width(nn),'%.1f') cl op deblank(num2str(node_size{nn},'%.2f ')) cl '\n'];
        end
        bgwrite{nn+3} = ['\n\n*Edges \t' num2str(sum(sum(bg_data))) '\n'];
        bgwrite{nn+4} = '*Example: source_node \t sink_node \t linecolor \t linewidth\n';
        
        bgLine = length(bgwrite)+1;
       
        for bgr = 1:number_of_nodes
            for bgc = 1:number_of_nodes
                if bg_data(bgr,bgc)
                   bgwrite{bgLine} = [num2str(bgr) '\t' num2str(bgc) '\t' op deblank(num2str(edge_color{bgr,bgc},'%.2f ')) cl op num2str(edge_width(bgr,bgc),'%.2f') cl '\n'];
                    bgLine = bgLine+1;
                   %fprintf(fid, [op b.LineColor cl op b.LineWidth cl]);
                end
            end
        end
      
end
end
end

function [molList molmap A] = get_molList(pwnum, levels, udp_s, handles)

A = addpw_levels(pwnum,levels, handles.PwData.DB);
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
    %probe_links = handles.DataStruct(1).Data(3).Data{1};
for j = 1:number_of_mols
    class = get_class(molList.molType{j}, molList.molLink{j});%classes: 1:protein, 2:complex, 3:compound, 4:other
    [molProbs(j,1:length(udp_s)) compProbs{j} mLink]= get_probs(handles, class, molList.molLink{j}, molList.molNum{j}, udp_s, []);
    molList.molLink{j}= mLink;
end
else
    molProbs(1:number_of_mols,1:length(udp_s)) = NaN;
    compProbs = cell(1,number_of_mols);
end
%add probabilities to list of mols
molList = [molList dataset({molProbs, 'molProbs'}, {compProbs', 'compProbs'}) ];
end

function [probs cProbs parsedlink]= get_probs(handles, class, link, num, snums, cProbs)
    %use molecule link to find associated probes, and compile upd/down probabilities
    %for those probes into a probability for the molecule
   
        if class == 1
            cProbs = [];
            parsedlink = get_parsedlink(link, handles.PwData.UP2LL);
            if ~isempty(parsedlink)
                proberows = get_proberows(handles, parsedlink);
                if proberows
                    Probs = get_pr(handles, proberows, snums);
                    probs(1,1:length(snums)) = 1;
                        for cc = 1:size(Probs,1)
                            probs = probs .* Probs(cc,:);
                        end
                else
                    probs = NaN;
                    %s = [s; {sprintf('\t Molecule: %d, link: %s, no matching proberows', num, link)}];
                end
            else
                probs = NaN;
                %s = [s; {sprintf('\t Molecule: %d, invalid link', num)}];
            end
            
        elseif class == 2
            parsedlink = get_parsedlink(link, handles.PwData.UP2LL);
            % complex_data = getcomplex(num2str(num));
            % complex_data = handles.PwData.Complexes([handles.PwData.Complexes{:,1}] == num,:);
            % found_prob = 0;
            % probs(1,1:length(snums)) = 1;
            % parsedlink = [];
            if ~isempty(parsedlink)
                proberows = get_proberows(parsedlink);
                if proberows
            % try
                % for cc = 1:size(complex_data,1)
                    % class = get_class(complex_data{cc,2}, complex_data{cc,5});
                    % [cProbs(cc,1:length(snums)) a cLink] = get_probs(handles, class, complex_data{cc,5}, complex_data{cc,4}, snums, cProbs);
                    cProbs = get_pr(handles, proberows, snums);
                    probs(1,1:length(snums)) = 1;
                    % if ~isnan(cProbs(1))
                        % found_prob = 1;
                        % probs = probs .* cProbs(cc,:);
                        % parsedlink= [parsedlink;cLink];
                        for cc = 1:size(cProbs,1)
                            probs = probs .* cProbs(cc,:);
                        end
                    % end
                    
                % end
                % if isempty(complex_data)
                %    cProbs = [{'No complex data'} num2cell(NaN*ones(1,length(snums)))];
                % else
                %     cProbs = [complex_data(:,3) num2cell(cProbs)];
                % end
                % if ~found_prob %if none of the components has a probability
                else
                    probs = NaN;
                    %s= [s; {sprintf('\t Complex: %d, no components probs found', num)}];
                end
            % catch
            else
                probs = NaN;
                % cProbs = [{'No complex data'} num2cell(NaN*ones(1,length(snums)))];
                 %s = [s; {sprintf('\t Complex: %d, no component data', num)}];
            end

        elseif class == 3
             cProbs = [];
             parsedlink = [];
            probs = 1; %compounds are assumed to always be present
        elseif class == 4
             cProbs = [];
             parsedlink = [];
            probs = 1;
             %s = [s; {sprintf('\t Molecule: %d, no identifier', num)}];
        end
end

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

function parsedlink = get_parsedlink(str, UP2LL)
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
                parsedlink{i,2} = str2double(parsedlink{i,1});
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

function proberows = get_proberows(handles, parsedlink)
    %search linking table for probes associated with LL# or UP#
    proberows = [];
    probe_links = handles.ChipData.Links;

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

function pr = get_pr(handles, proberows, snums)
        for n = 1:length(proberows) %get probabilities from UDP for each associated probe
            p(n,:) = handles.DataStruct(1).Data(2).Data(proberows(n),snums);
        end
        % pr = max(p, [], 1); %for each sample, take the max probability in the probeset
        pr = p;
end

function [i_activity i_consistency omit_i mult_out] = pw_calc(A, molList, molmap, udp_s)
[intxn_levels ia intxn_map]= unique(A(:,11));
number_of_intxns = length(intxn_levels);

%set up variables
i_activity = ones(number_of_intxns,length(udp_s)); 
i_consistency = ones(number_of_intxns,length(udp_s)); 
omit_i = [];
% skip_rest = 0;
mult_out = [];
%for each interaction, use mol roles (eg. input/output) and probabilities
%to calculate activity and consistency
for in = 1:number_of_intxns; 
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
           %if strcmp(m_roles{j}, 'output')
               %nodes = [nodes; in m_index(j)];
               i_con(:,j) = i_activity(in,j).* m_probs(outs,j) + (1 - i_activity(in,j)).*(1-m_probs(outs,j));
       end
       else
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

function [] = nodeclick(node, handles)
%opens cgap gene webpage, or depicts components of a complex

if strcmp(node.Shape, 'box')
    %do nothing for intxn nodes
    if strcmp(node.Description, '')
        web(strcat('http://cgap.nci.nih.gov/Genes/RunUniGeneQuery?PAGE=1&ORG=Hs&SYM=&PATH=&TERM=', node.ID));
    else
         
        B = node.UserData;
        
        num_comps = size(B,1);
        for ii = 1:num_comps
               new_node_labels{ii,1} = ['comp', num2str(ii), ': ', B{ii,1}];
        end
        new_node_labels = [new_node_labels; node.ID];

new_bg_data = zeros(num_comps+1);
new_bg_data(1:size(B,1), num_comps+1) = 1;
new_BG = biograph(new_bg_data, new_node_labels, 'ShowWeights', 'off', 'ShowTextInNodes', 'id', 'NodeCallback', {@(newnode)web(strcat('http://cgap.nci.nih.gov/Genes/RunUniGeneQuery?PAGE=1&ORG=Hs&SYM=&PATH=&TERM=', newnode.ID))});

if get(handles.data_typePanel, 'SelectedObject') ~= handles.cnaRadio
    node_probs = [B{:,2}];
    for ii = 1:num_comps
        set(new_BG.nodes(ii), 'Label', strcat('probability: ', num2str(node_probs(ii))));
    end
    set(new_BG.nodes(end), 'Label', node.Label);
end
view(new_BG)
    end
end
end

function returnButton_Callback(hObject, eventdata, handles)
% hObject    handle to returnButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume
end

function [joint_pw intxn_levels]= addpw_levels(pwnum, levels, DB)

A = DB([DB{:,2}] == pwnum,:);
A(strcmp('null',A(:,4)),4)= A(strcmp('null',A(:,4)),6); 
pathway_name = A{1,1};
size_A = size(A,1);
joint_intxn_ids = getlabels(nominal(A(:,11)));
joint_mol_ids = [];
joint_pw = A;

intxn_levels = ones(length(joint_intxn_ids),1);

for L = 2:levels
    new_mols = setdiff(unique([joint_pw{:,5}]), joint_mol_ids);
    for ii = 1:length(new_mols)        
        intxn_count = length(joint_intxn_ids);
        B = DB([DB{:,5}] == new_mols(ii),:);
        B_intxn_ids = getlabels(nominal(B(:,11)));
        joint_intxn_ids = [joint_intxn_ids setdiff(B_intxn_ids, joint_intxn_ids)];
        intxn_levels = [intxn_levels; L*ones(length(joint_intxn_ids)-intxn_count,1)];
    end


    for ii = find(intxn_levels == L)'
        joint_pw = [joint_pw; DB(strmatch(joint_intxn_ids{ii}, DB(:,11)) ,:)];
    end
        joint_mol_ids = [joint_mol_ids new_mols];
end

end

% function cnaRadio_Callback(hObject, eventdata, handles)
% hObject    handle to cnaRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cnaRadio

% if get(hObject, 'Value')
%     set(handles.PIDRadio, 'Value', 0)
% end
% end

% function acRadio_Callback(hObject, eventdata, handles)
% hObject    handle to acRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of acRadio

% if get(hObject, 'Value')
%     set(handles.PIDRadio, 'Value', 0)
% end
% end
% Blank Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function textCheck_Callback(hObject, eventdata, handles)
end
function excelRadio_Callback(hObject, eventdata, handles)
end
function textRadio_Callback(hObject, eventdata, handles)
end
function bgCheck_Callback(hObject, eventdata, handles)
% hObject    handle to bgCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bgCheck

end
function addlevelCheck_Callback(hObject, eventdata, handles)
% hObject    handle to addlevelCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of addlevelCheck
end
function samplesListbox_Callback(hObject, eventdata, handles)
% hObject    handle to samplesListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns samplesListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from samplesListbox
end
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
function pw_Listbox_Callback(hObject, eventdata, handles)
% hObject    handle to pw_Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pw_Listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pw_Listbox
end
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
