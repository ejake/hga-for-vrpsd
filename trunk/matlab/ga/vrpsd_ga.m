function [opt_pol min_exp] = vrpsd_ga(instance,pop_size,num_iter,epsilon,show_prog,show_res)
%VRPSD_GA Vehilce Routing Problem Whit Stochastic Demands (VRPSD) Genetic Algorithm (GA)
%   Finds a (near) optimal solution to the VRPSD by setting up a GA to search
%   for the shortest route (least expected distance for the vehicle to travel to
%   each customer and return to the depot)
%
% Input:
%     INSTANCE (object instance)     
%     POP_SIZE (scalar integer) is the size of the population (should be divisible by 4)
%     NUM_ITER (scalar integer) is the number of desired iterations for the algorithm to run
%     EPSILON
%     SHOW_PROG (scalar logical) shows the GA progress if true
%     SHOW_RES (scalar logical) shows the GA results if true
%
% Output:
%     OPT_POL (integer array) is the best policy found by the algorithm
%     MIN_EXP (scalar float) is the cost of the best policy
%
% Example:
%     n = 50;
%     xy = 10*rand(n,2);
%     a = meshgrid(1:n);
%     dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),n,n);
%     pop_size = 60;
%     num_iter = 1e4;
%     show_prog = 1;
%     show_res = 1;
%     [opt_rte,min_dist] = tsp_ga(xy,dmat,pop_size,num_iter,show_prog,show_res);
%
% Author: Andres Jaque
% Email: rajaquep@gmail.com
% Release: 1.0
% Release Date: 3/6/08


% Verify Inputs
n = instance.n;
[nr,nc] = size(instance.d);
if n ~= (nr-1) || n ~= (nc-1)
    error('Invalid instance inputs (distance)!')
end

% Sanity Checks
pop_size = 4*ceil(pop_size/4);
num_iter = round(real(num_iter(1)));
show_prog = logical(show_prog(1));
show_res = logical(show_res(1));

% Initialize the Population
pop = initPopulationCyclic(pop_size, n); %if pop_size >= n, pop_size - n individuals would be repeated

% Run the GA
if epsilon <= 0
    epsilon = 0.0000001;
end
cnEpsilon = 0;
gain = Inf;
offspring_pop = Individual.empty(pop_size,0); % the size of offspring is almost the size of population
offspring_counter = 0;
new_pop = Individual.empty(pop_size,0);
global_min = Inf;
total_dist = zeros(1,pop_size);
dist_history = zeros(1,num_iter);

if show_prog
    pfig = figure('Name','Current Best Solution','Numbertitle','off');
    count_subplot =1;
end

%Apply RA to a member of population
[pi cyEd]  = rollout (instance, State(instance.n, instance.Q), pop(1).tour);
%if cyclic heuristic is used  in RA expected distance is already computed
for i=1:n
    pop(i).expected_distance = cyEd(i);
    pop(i).rolledout = true;
end
offspring_pop(1) = Individual();
offspring_pop(1).policy = pi;
offspring_counter = offspring_counter+1;

iter = 0;
while ((iter < num_iter) && gain > epsilon) % stopping criterion
    % Evaluate Each Population Member (Calculate Expected Distance)
    for p = 1:pop_size
        if pop(p).expected_distance == Inf
            total_dist(p).expected_distance = backwardExpectedDistance([0 pop(p,:)], instance);            
        end
        total_dist(p) = pop(p).expected_distance;
    end

    % Find the Best Route in the Population
    [min_dist,index] = min(total_dist);
    dist_history(iter) = min_dist;
    gain = (global_min-min_dist)/global_min; % Inf/Inf
    if min_dist < global_min
        global_min = min_dist;
        opt_rte = pop(index,:);
        if show_prog && (iter == 1 || iter == floor(iter/4) || iter == floor(iter/2) || floor(3*iter/4) )
            % Plot the Best Route
            figure(pfig);
            subplot(1,4,count_subplot);
            % imagesc(dmat(opt_rte,opt_rte))
            imagesc(pop);            
            title(sprintf('Total Distance = %1.4f, Iteration = %d',min_dist,iter));
            count_subplot = count_subplot+1;
        end
    end

    % Genetic Algorithm Operators - Mutation
    rand_pair = randperm(pop_size);%random selection of individuals to mutate
    %define number of individuals to mutate (<= pop_size)
    %probability of mutation p_m
    p_m = 0.5;
    if (rand() <= p_m)
        for p = 1:ceil(randi(pop_size-1,1)*rand())% # individuals to mutate
            rtes = pop(p).tour;
            dists = total_dist(rand_pair(p-3:p));
            [ignore,idx] = min(dists);
            best_of_4_rte = rtes(idx,:);

            for k = 1:4 % Mutate the Best to get Three New Routes
                tmp_pop(k,:) = best_of_4_rte;

            end
            new_pop(p-3:p,:) = tmp_pop;
        end        
    end
    
    % Genetic Algorithm Operators - Crossover
    
    % Local search
    
    % selection
    pop = new_pop;
    if(gain < epsilon)
        cnEpsilon = cnEpsilon + 1;
        if(cnEpsilon > num_iter*0.1)
            showResults(show_res,instance, opt_rte, min_dist, num_iter, pop, dist_history);
            return
        end
    end 
end
    showResults(show_res,instance, opt_rte, min_dist, num_iter, pop, dist_history);
% Return Outputs
%if nargout
%    varargout{1} = opt_rte;
%    varargout{2} = min_dist;
end

function showResults(show_res,instance, opt_rte, min_dist, num_iter,pop, dist_history)
    if show_res
        xy = zeros(instance.n+1,2);
        for i=1:instance.n
            xy(i+1,:) = instance.Cust(i).location;
        end
        % Plots the GA Results
        figure('Name','TSPGA','Numbertitle','off');
        subplot(2,2,1);
        plot(xy(:,1),xy(:,2),'k.');
        title('City Locations');
        subplot(2,2,2);
        %imagesc(instance.d);
        %title('Distance Matrix');
        imagesc(pop);
        title('Population');
        subplot(2,2,3);        
        rte = [1 opt_rte+1 1];
        plot(xy(rte,1),xy(rte,2),'r.-');
        title(sprintf('Total Distance = %1.4f',min_dist));
        subplot(2,2,4);
        plot(dist_history,'b','LineWidth',2);
        title('Best Solution History');
        set(gca,'XLim',[0 num_iter+1],'YLim',[0 1.1*max([1 dist_history])]);
    end

end
