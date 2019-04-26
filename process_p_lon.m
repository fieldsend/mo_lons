function [V,B,Adj,EE,C] = process_p_lon(X,Y,state,neighbours)

%  [V,B,EE,C,PO] = process_p_lon(X,Y,w,YY,state, neighbours)
%
% Processes data outputted from exaustive_generate_lon to create 
% the matrices needed to plot the PLON
%
% see generate_GECCO_2019_plots for example usage
%
% state holds destination of the walk from the index
%
% Jonathan Fieldsend, University of Exeter, 2019
% See license information in package, available at 
% https://github.com/fieldsend/mo_lons

min_state = cell(0);
    for j=1:length(state)
        contained = false;
        for k=1:length(min_state)
            if (length(min_state{k})==length(state{j}))
               if sum(min_state{k}==state{j})==length(state{j})
                  contained = true;
                  break;
               end
            end
        end
        if (contained==false)
           min_state{length(min_state)+1} = state{j}; 
        end
    end
    unique_state = min_state;
    % min_state now contains each unique PLO returned by PLS, not need to
    % pull out additional unique intersections
    for k=1:length(min_state)
        for j=1:length(min_state)
            if (j~=k)
                U = intersect(min_state{k},min_state{j});
                U = sort(U); % might not be required?
                already_seen = false;
                if isempty(U)==false
                    for m=1:length(unique_state)
                        if length(unique_state{m})==length(U) && (sum(U==unique_state{m})==length(U))
                            already_seen = true;
                            break;
                        end
                    end
                else
                    already_seen = true;
                end
                if already_seen == false
                    unique_state{length(unique_state)+1} = U;
                end
            end
        end
    end
min_state    
V =   min_state;  
    %V = unique_state; % vertices
B = zeros(length(V),1);
Adj = zeros(length(V),length(V));
for i=1:length(state)
    fprintf('processing state %d of %d \n',i,length(state));
    Z = zeros(1,length(V));
    index = 0;
    for j=1:length(V)
        K = intersect(state{i},V{j});
        if (length(K)==length(V{j}))
            if (length(state{i})==length(V{j})) %same
                B(j) = B(j)+1;
                index = j;
            else %proper subset
                Z(j) = 1;
            end
        end
        
    end
    % now increase weight between vertices which all have path from ith
    Adj(Z==1,index) = Adj(Z==1,index)+1;

end


% now colour each vertex

C = zeros(1,length(V));
PO = zeros(length(X),1);
[n,m] = size(Y);
fprintf('Determining Pareto optimal');
for i=1:length(V) % for each vertex
    for j=1:length(V{i}) % for each vertex member
        wdom = zeros(n,1);
        eq = zeros(n,1);
        for k=1:m
           wdom = wdom + (Y(V{i}(j),k)>=Y(:,k));
           eq = eq + (Y(V{i}(j),k)==Y(:,k));
        end
        if (sum(wdom==m)==sum(eq==m))
           % Pareto optimal
           PO(V{i}(j)) = 1;
        end
    end
end

for i=1:length(V) % for each vertex
    C(i) = sum(PO(V{i})); % sum number of vertex members who are pareto optimal
end



% add esacpe edges to links
EE = zeros(length(B)); % matrix of edges between basins, filled with zeros initially
for i=1:length(B) % process each vertex
    for j=1:length(V{i}) % process each vertex member
        K = neighbours(V{i}(j),:); % get neighbours of each vertex element
        for k = K
           if ((k>0) && (k<=length(X))) % if a valid neighbour index
               % now need to find which vertices each element of state{k}
               % is in
               for m=1:length(B)
                    U = intersect(state{k},V{m}); 
                    if length(U)==length(V{m}) % if leads to superset of V{m}
                        EE(i,m) = EE(i,m)+1;
                    end
               end
           end
        end
    end
end

EE=EE+Adj;

