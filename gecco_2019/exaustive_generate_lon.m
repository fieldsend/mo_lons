function [X,Y,state,neighbours,B,YY] = exaustive_generate_lon(cost_function, args, p_lon, upper_bound, lower_bound, resolution, num_of_objectives)

% exaustive_generate_lon(cost_function, type, max, min, resolution)
%
% INPUTS
%
% cost_function = handle to cost function to be minimise d, should expect a 
%                 vector of length(max_values) and an optional args argument
% args = any additional structural arguments needed by the cost function
% p_lon = LON pair type, true means PLON data generated, false means DNON
%                 data generated
% upper_bound = array of maximum input values
% lower_bound = array of minimum input values
% resolution = number of slices to grid each dimension
% num_of_objectives = number of objectives in problem
%
% OUTPUTS
%
% X = design locations queried
% Y = fitness at design locations queried
% state = cell array containing attractor(s) for each design location
%         queries 
% neighbours = neighbour indices for each design location (an Inf entered 
%         for when there is not a fixed number of neighbours for all designs  
% B = DLON basin weights, weighetd number of hillclimbs leading to location 
%         (DNON only, empty for PLON)
% YY = matrix of B values (DNON only, empty for PLON)
%
%
% Function generates data for a PLON or DNON based on gridding the continuous
% space cost_function at the given resolution, in a box constraint space
%
%
% Jonathan Fieldsend, University of Exeter, 2019
% See license information in package, available at 
% https://github.com/fieldsend/mo_lons


d = length (upper_bound);
max_integer = resolution-1;
step_size = (upper_bound-lower_bound)/max_integer;

number_of_samples = resolution^d;
X = zeros(number_of_samples,d);
Y = zeros(number_of_samples,num_of_objectives);
coordinate_multiplier = get_coordinate_multiplier(d,resolution); 


%evaluate all points on the grid
for i=1:number_of_samples
    [values_as_integers] = convert_index_to_decision_vector(i,d,resolution,coordinate_multiplier);
    X(i,:) = values_as_integers;
    decision_vector = lower_bound+values_as_integers.*step_size;
    Y(i,:) = feval(cost_function,decision_vector,args);
end
fprintf('evaluated all fitnesses\n');
if (p_lon)
    [state, neighbours] = generate_p_lon(X,Y,d,num_of_objectives,number_of_samples,resolution,coordinate_multiplier);
    B=[];
    YY=[]; 
else
    [B, state, neighbours] = generate_d_lon(X,Y,d,num_of_objectives,number_of_samples,resolution,coordinate_multiplier);
    % now convert to matrix
    YY = zeros(resolution,resolution);
    for i=1:number_of_samples
        [values_as_integers] = convert_index_to_decision_vector(i,d,resolution,coordinate_multiplier);
        YY(values_as_integers(1)+1,values_as_integers(2)+1) = B(i);
    end
end





end

function [P, state, neighbours] = generate_d_lon(X,Y,d,m,number_of_samples,resolution,coordinate_multiplier)

neighbours = zeros(number_of_samples,2*d); % get the neighbours for each location
%support = zeros(number_of_samples,1); % weighted number of elements who hill-climbs to this location
P = zeros(number_of_samples,1);

for i=1:number_of_samples
    % get the neighbours of a sample
    neighbours(i,:) =  get_integer_neighbour_indices(X,i,d,resolution,coordinate_multiplier);
end
fprintf('calculated all neighbour locations\n');
state = cell(number_of_samples,1);

for i=1:number_of_samples
    [S, state] = pareto_hill_climb(X,Y,i,neighbours,d,m,coordinate_multiplier,state);
    [U,~,N2] = unique(S);
    for j=1:length(U)
       P(U(j)) = P(U(j)) + sum(U(N2)==U(j))/length(S); % weight for multiple paths to same destination 
    end
    state{i} = U;%S;
    if (rem(i,1)==0)
       fprintf('processed %d of %d \n',i,number_of_samples); 
    end
end
fprintf('concluded all exhaustive hillclimbs\n');

end

function [state, neighbours] = generate_p_lon(X,Y,d,m,number_of_samples,resolution,coordinate_multiplier)

neighbours = zeros(number_of_samples,2*d); % get the neighbours for each location
%P = zeros(number_of_samples,1);

for i=1:number_of_samples
    % get the neighbours of a sample
    neighbours(i,:) =  get_integer_neighbour_indices(X,i,d,resolution,coordinate_multiplier);
end
fprintf('calculated all neighbour locations\n');
state = cell(number_of_samples,1);

for i=1:number_of_samples
    already_processed = zeros(1,number_of_samples);
    already_processed(i) = 1;
    PLO = i;
    [PLO] = pareto_local_search(X,Y,neighbours,d,m,coordinate_multiplier,already_processed,PLO);
    U = unique(PLO); % check if actually needed?
    state{i} = sort(U); % PLO set resulting from ith solution
    if (rem(i,1)==0)
       fprintf('processed %d of %d \n',i,number_of_samples); 
    end
end
fprintf('concluded all exhaustive set searches\n');

end

function [PLO,already_processed] = pareto_local_search(X,Y,neighbours,d,m,coordinate_multiplier,already_processed,PLO)
%PLO
[n,m] = size(Y);
indices =[];
for i=PLO
    indices =  union(indices,neighbours(i,:));
end
indices = indices(indices>0 & indices<(length(X)+1))'; % strip out illegal
%indices
for j=indices % for each neighbour of current set
        if (already_processed(j)==0) % if not already processed previously
            dominated = false;
            for k=1:length(PLO)
                if (sum(Y(PLO(k),:) <= Y(j,:))==m) && (sum(Y(PLO(k),:) == Y(j,:))~=m)
                    dominated = true;
                    break;
                end
            end
            if (dominated==false)
                PLO = [PLO j]; % not dominated by PLO, so add 
            end
            already_processed(j) = 1;
        end
end
% now strip out dominated members of PLO
to_remove = [];
for k=1:length(PLO)
    dominated = false;
    for j=1:length(PLO)
        if (sum(Y(PLO(j),:) <= Y(PLO(k),:))==m) && (sum(Y(PLO(j),:) == Y(PLO(k),:))~=m)
            dominated = true;
            break;
        end
    end
    if (dominated==true)
        to_remove = [to_remove k];
    end
end
%'to_remove'
%PLO(to_remove)
PLO(to_remove) = [];

% now look at neighbours of PLO entries
indices =[];
for i=PLO
    indices =  union(indices,neighbours(i,:));
end
indices = indices(indices>0 & indices<(length(X)+1))'; % strip out illegal
if sum(already_processed(indices))<length(indices) % if any neighbours not yet processed
   [PLO,already_processed] = pareto_local_search(X,Y,neighbours,d,m,coordinate_multiplier,already_processed,PLO);
end

%PLO

end



function [S, state] = pareto_hill_climb(X,Y,index,neighbours,d,m,coordinate_multiplier,state)

% get the neighbours of a sample
indices =  neighbours(index,:);
K =[];
%fprintf('location %d \n', index);
for j=indices
    if j>0 && j<(length(X)+1)
        % add indices of dominating neighbours to K
        if (sum(Y(index,:) >= Y(j,:))==m) && (sum(Y(index,:) == Y(j,:))~=m)
            K = [K j];
        end
    end
end
% K contains indices of all locations which dominate the location at index
S=[];
if (isempty(K)==false)
    for j=K
        if isempty(state{j})
            S = [S pareto_hill_climb(X,Y,j,neighbours,d,m,coordinate_multiplier,state)];
        else
            S = [S state{j}]; % already processed destinations from jth location 
        end
    end
else
    S = [S index]; % refer to itself as no location dominates it in neighbourhood
end
state{index} = S;
end



function [values_as_integers] = convert_index_to_decision_vector(index,d,resolution,coordinate_multiplier)

% function converts an index into a decision vector of integers (with min value = 0 and 
% max value = resolution-1)


if exist('coordinate_multiplier','var')==0
    coordinate_multiplier = get_coordinate_multiplier(d,resolution);
end
% divide through by multiplier from highest to lowest to get vector
index = index-1; % first decrement back by one as values down to zero

values_as_integers = zeros(1,d);
for j=d:-1:1
    values_as_integers(j) = floor(index/coordinate_multiplier(j));
    index = rem(index,coordinate_multiplier(j));
end


end

function [index] = convert_decision_vector_to_index(values_as_integers,resolution,coordinate_multiplier)

% function converts a decision vector of integers (with min value = 0 and 
% max value = resolution-1) into a single index
%
% Exploits gridding index process detailed in e.g. PAES algorithm by
% Knowles and Corne

% Knowles, J.D., Corne, D.W. (2000) Approximating the nondominated front
% using the Pareto Archived Evolution Strategy. Evolutionary Computation, 
% 8(2), pp. 149-172

d = length(values_as_integers);
if exist('coordinate_multiplier','var')==0
    coordinate_multiplier = get_coordinate_multiplier(d,resolution);
end
values_as_integers=values_as_integers.*coordinate_multiplier;
index=sum(values_as_integers)+1;

end    
    
function indices = get_integer_neighbour_indices(X,i,d,resolution,coordinate_multiplier)
    indices = zeros(1,2*d);
    v = [1 -1];
    index = 1;
    for j=1:d
        for k = 1:2
            x = X(i,:);
            x(j) = x(j)+v(k);
            if (x(j)>=resolution)
               x(j) = inf; 
            end
            if (x(j)<0)
               x(j) = inf; 
            end
            indices(index) = convert_decision_vector_to_index(x,resolution,coordinate_multiplier);
            index = index+1;
        end
    end
end


function coordinate_multiplier = get_coordinate_multiplier(d,resolution)


% function gets coordinate multipler for converting vector (with min value = 0 and 
% max value = resolution-1) into a single index

coordinate_multiplier = ones(1,d);
for i=2:d
    coordinate_multiplier(i)=resolution*(i-1);
end

end
