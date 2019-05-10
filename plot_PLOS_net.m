function [DNO] = plot_PLOS_net(X,Y,state,neighbours,edge_weight,node_weight)

% [DNO] = plot_PLOS_net(X,Y,state,neighbours,edge_weight,node_weight)
%
% Function plots the PLOS-net for the given parameters, outputted from the 
% exaustive_generate_lon function. See e.g. generate_GECCO_2019_plots for 
% example usage.
%
% Jonathan Fieldsend, University of Exeter, 2019
% See license information in package, available at 
% https://github.com/fieldsend/mo_lons


[~,k] = size(Y);

% function plots the PLOS-net for the given inputs

DNO = state{1};
for i=2:length(state)
   DNO = [DNO state{i}]; 
end
DNO = unique(DNO);

% DNO now contains indices of dominance neutral optima


fprintf('Plot graph\n');
s = zeros(1,length(DNO)*length(DNO));
t = zeros(1,length(DNO)*length(DNO));
weights = zeros(1,length(DNO)*length(DNO));
index =1;
for i=1:length(DNO)
    for j=1:length(DNO)
        if j>=i
            s(index) = i;
            t(index) = j;
            if sum(neighbours(DNO(i),:)==DNO(j))>0
                weights(index) = 1;
            end
            index = index+1;
        end
    end
end

% unhelpful error message from matlab (negative weights not allowed in graph)
% but actually 0 valued weights need removing
I = find(weights==0);
s(I)=[];
t(I)=[];
weights(I)=[];


%basin_size = [3+2/3, 3+2/3, 3+2/3, 5];
G = graph(s,t,weights);
LWidths = edge_weight*G.Edges.Weight/max(G.Edges.Weight);

% L = cell(1,length(Q));
% for i=1:length(L)
%     L{i} = int2str(Q(i));
% end
H=plot(G,'Layout','force','LineWidth',LWidths,'NodeLabel',{});%,'NodeLabel',names);


max_y = max(Y(DNO,k));
min_y = min(Y(DNO,k));

for i=1:length(DNO)
    highlight(H,i,'MarkerSize',node_weight,'NodeColor',0.1+abs(ones(1,3)*(Y(DNO(i),k)-min_y)/(max_y-min_y)-1)*0.5);
end

set(gca,'YTick',[])
set(gca,'XTick',[])

end
