function [V,B,EE,C,PO] = process_d_lon(X,Y,w,YY,state, neighbours)


% state holds destination of all walks from the index


V=cell(1,1);
v_index=1;
processed = zeros(1,length(w));
for i=1:length(w)
    K=[];
    if (w(i)>0) % if attractor at index i
        if (processed(i)==0) % if not represented in any vertices
            K = [K i];
            processed(i) = 1;
            [N, processed] = recursively_get_neighbours(i,w,neighbours,processed);
            [K] = [K N];
            K = unique(K);
            V{v_index} = K;
            v_index = v_index+1;
        end
    end
end
% V now contains each distinct DNO set
B = zeros(1,length(V));
for i=1:length(V)
   B(i) = sum(w(V{i})); 
end
% M now contains the basin size for each vertex

EE = zeros(length(B)); % matrix of edges between basins, filled with zeros initially
for i=1:length(B) % process each vertex
    for j=1:length(V{i})
        K = neighbours(V{i}(j),:); % get neighbours of each vertex element
        for k = K
           if ((k>0) && (k<=length(X))) % if a valid neighbour index
               % now need to find which vertices each element of state{k}
               % is in
               for n=state{k}
                   for m=1:length(B)
                       if sum(V{m}==n)>0
                            EE(i,m) = EE(i,m)+1;
                       end
                   end
               end
           end
        end
    end
end

% EE now contains edge escape numbers


C = zeros(1,length(B));
PO = zeros(length(X),1);
[n,m] = size(Y);
for i=1:length(B) % for each vertex
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
PO
Y
for i=1:length(B) % for each vertex
    C(i) = sum(PO(V{i})); % sum number of vertex members who are pareto optimal
end

%C now contains numbers of global pareto optima in each basin


end

function [K, processed] = recursively_get_neighbours(index,w,neighbours, processed)

I = neighbours(index,:);
I = I((I>0) & (I<=length(w)));
K = [];
for j=I
    if (w(j)>0)
        if (processed(j)==0)
            processed(j)=1;
            K = [K j];
            [N, processed] = recursively_get_neighbours(j,w,neighbours,processed);
            [K] = [K N];
        end
    end
end
        

end
