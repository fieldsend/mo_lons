function [s, t, weights] = generate_lon(basin_size, Q,TP, mode_indices,edge_weight,node_weight,EC)

%
% basin_size = 1 by number_of_vertices array of basin sizes
% Q = 1 by number_f_vertices array of fitnesses of basin modes
% TP = number_of_vertices by number_of_vertices matrix of weights between
%      vertices
% mode_indices -- NOT USED CURRENTLY
% edge_weight = scalar, rescales all edge widths
% node_weight = scalar, rescales all vertex radii
%
%
%
%
fprintf('Plot graph\n');
%E = [0 0; 0 1; 0 2; 0 3; 1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
%s = E(:,1)'+1; t=E(:,2)'+1;
%weights = [50 20 20 5 50 20 5 50 5 50];
s = zeros(1,length(basin_size)*length(basin_size));
t = zeros(1,length(basin_size)*length(basin_size));
%ec = zeros(1,length(basin_size)*length(basin_size));
weights = zeros(1,length(basin_size)*length(basin_size));
index =1;
for i=1:length(basin_size)
    for j=1:length(basin_size)
        s(index) = i;
        t(index) = j;
        weights(index) = TP(i,j);
        index = index+1;
    end
end

% unhelpful error message from matlab (negative weights not allowed in graph)
% but actually 0 valued weights need removing
I = find(weights==0);
s(I)=[];
t(I)=[];
weights(I)=[];


%basin_size = [3+2/3, 3+2/3, 3+2/3, 5];
G = digraph(s,t,weights);
LWidths = edge_weight*G.Edges.Weight/max(G.Edges.Weight);

L = cell(1,length(Q));
for i=1:length(L)
    L{i} = int2str(Q(i));
end

H=plot(G,'Layout','force','LineWidth',LWidths,'NodeLabel',{});%,'NodeLabel',names);


max_y = max(Q);
min_y = min(Q);

if min_y==max_y
    for i=1:length(basin_size)
            highlight(H,i,'MarkerSize',node_weight*basin_size(i)/sum(basin_size),'NodeColor',[0 0 0]);
    end
else
    for i=1:length(basin_size)
        highlight(H,i,'MarkerSize',node_weight*basin_size(i)/sum(basin_size),'NodeColor',0.1+abs(ones(1,3)*(Q(i)-min_y)/(max_y-min_y)-1)*0.5);
    end
end
if (exist('EC','var')==1)
    for i=1:length(basin_size)
        for j=1:length(basin_size)
            if (EC(i,j)>0)
                i
                highlight(H,[j i], 'EdgeColor','r');
                
            end
        end
    end
end


set(gca,'YTick',[])
set(gca,'XTick',[])
