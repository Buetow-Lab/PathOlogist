function [NodeLabels IntxnMatrix MolVals MolSelect] = create_pw(GeneList,Metrics, MetricsLabels)
%GeneList: list of Gene Symbols
%Metrics: value for many genes
%MetricsLabels: rowheaders for Metrics


%read in PID database
fid = fopen('db.txt');
fd = textscan(fid, '%s%9.0f%s%s%9.0f%s%s%s%s%s%s%s', 'delimiter', '\t');
fclose(fid);
handles.PwData.DB = horzcat(fd{1},num2cell(fd{2}),fd{3:4},num2cell(fd{5}),fd{6:end});

%read in PID complex data
fid = fopen('complexes.txt');
fc = textscan(fid, '%9.0f%s%s%9.0f%s%s%s', 'delimiter', '\t');
fclose(fid);
handles.PwData.Complexes = horzcat(num2cell(fc{1}),fc{2:3},num2cell(fc{4}),fc{5:end});

%read in LocusLink-UniProt linking data
up2ll = importdata('UP2LL.txt', '\t');
handles.PwData.UP2LL = [up2ll.textdata num2cell(up2ll.data)];

A = handles.PwData.DB(ismember(handles.PwData.DB(:,4),GeneList),:);
A = handles.PwData.DB(ismember(handles.PwData.DB(:,11),A(:,11)),:);

[molList molList_map] = get_molList(A);
[NodeLabels IntxnMatrix] = pw_draw();

function [molList molList_map A] = get_molList(A)

D = dataset({A, 'pwName', 'pwNum', 'molType', 'molName', 'molNum', 'molLink', 'c7', 'c8', 'c9', 'molRole', 'c11', 'intxnType'});

%create list of mols involved in pathway, and a map to their rows in A
[N iaa molList_map]= unique([A{:,5}]);
number_of_mols = length(N);
N = A(iaa,4);
Arows = iaa;
molList = D(Arows, [4 5 3 6]);

clear D

for j = 1:number_of_mols
    try
        molVal(j) = Metrics(strmatch(molList.molName{j}, MetricsLabels));
    catch
        molVal(j) = NaN;
    end
end

%add probabilities to list of mols
molList = [molList dataset({molVal', 'molVals'}) ];
end
function [NodeLabels IntxnMatrix molprobs mol_cnas] = pw_draw()

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

%set up outline colors

    molList.outline = cell(size(molList,1),1);
    molList.outline(ismember(molList.molName, GeneList)) = {[1 .92 0]};
            
 



%assign each mol a probability

mol_probs = molList.molVals(molmap,:);
mol_probs = mol_probs./max(mol_probs);
num_samples = size(mol_probs,2);



mol_cnas = molList.outline(molmap,:);



        edge_color(logical(bg_data)) = {[.75 .75 .75]};
        edge_width(logical(bg_data)) = .5;
  


    node_color=cell(number_of_nodes,1); outline_color=cell(number_of_nodes,1); outline_width=zeros(number_of_nodes,1);shape=cell(number_of_nodes,1);node_size=cell(number_of_nodes,1);
 
    for ii = 1:number_of_nodes
        intxn_num = ii - number_of_mols;
    if intxn_num>0
        node_color{ii} = [.75 .85 1];
        outline_color{ii} = [.4 .4 .4];
        outline_width(ii) = .5;
        shape{ii} = 'circle';
        node_size{ii} = [20 20];
       
    else
        shape{ii} = 'box';
        node_size{ii} = [7*length(molLabels{ii}) 20];
        if  ~isnan(mol_probs(ii))
            node_color{ii} = [1-mol_probs(ii) 1 1];
        else
            node_color{ii} = [.85 .85 .85];
        end
        if ~isempty(mol_cnas{ii})
            outline_color{ii} = mol_cnas{ii};
            outline_width(ii)=2.5;
        else
            outline_color{ii} = [0 0 0];
            outline_width(ii) = .5;
        end
    end
    
 end
       

    make_BG()

%if tc
    %bgwrite{ss} = bg2text();
   % bgwrite = [];
%else
 %   bgwrite = [];

%end
IntxnMatrix = bg_data;
NodeLabels = node_names;

function make_BG()
%make biograph
   
    BG = biograph(bg_data, node_names, 'ShowWeights', 'off', 'EdgeFontSize', 7, 'LayoutType', 'hierarchical');

    %set node properties
    for ii = 1:number_of_nodes;
        set(BG.nodes(ii), 'Color', node_color{ii}, 'LineColor', outline_color{ii}, 'LineWidth', outline_width(ii), 'TextColor', [0 0 0], 'Shape', shape{ii});
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
            %set(edge, 'Label', loc_matrix{kk,3}(1:3));
        end
    end

    %set callback to open gene webpage, or show complex components
    set(BG, 'NodeCallback', {@(node)nodeclick(node, handles)});

    dolayout(BG)


    %make intxn node size separately after initial layout
    
        set(BG, 'NodeAutoSize', 'off');
        for ii = number_of_mols+1:number_of_nodes
            set(BG.nodes(ii), 'Size', node_size{ii});
            
        end
        
        for ii = 1:number_of_mols
            try
                set(BG.nodes(ii), 'Label', ['value: ' num2str(mol_probs(ii))]);
            catch
                set(BG.nodes(ii), 'Label', 'unknown');
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