function [opt_rte min_dist] = vrpsd_ga(instance,pop_size,num_iter,show_prog,show_res)
%TSP_GA Traveling Salesman Problem (TSP) Genetic Algorithm (GA)
%   Finds a (near) optimal solution to the TSP by setting up a GA to search
%   for the shortest route (least distance for the salesman to travel to
%   each city exactly once and return to the starting city)
%
% Summary:
%     1. A single salesman travels to each of the cities and completes the
%        route by returning to the city he started from
%     2. Each city is visited by the salesman exactly once
%
% Input:
%     XY (float) is an Nx2 (or Nx3) matrix of cities
%     DMAT (float) is an NxN matrix of point to point distances/costs
%     POP_SIZE (scalar integer) is the size of the population (should be divisible by 4)
%     NUM_ITER (scalar integer) is the number of desired iterations for the algorithm to run
%     SHOW_PROG (scalar logical) shows the GA progress if true
%     SHOW_RES (scalar logical) shows the GA results if true
%
% Output:
%     OPT_RTE (integer array) is the best route found by the algorithm
%     MIN_DIST (scalar float) is the cost of the best route
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
% Author: Joseph Kirk
% Email: jdkirk630@gmail.com
% Release: 2.0
% Release Date: 8/23/08


% Verify Inputs
n = instance.n;
[nr,nc] = size(instance.d);
if n ~= (nr-1) || n ~= (nc-1)
    error('Invalid XY or DMAT inputs!')
end

% Sanity Checks
pop_size = 4*ceil(pop_size/4);
num_iter = max(2,round(real(num_iter(1))));
show_prog = logical(show_prog(1));
show_res = logical(show_res(1));

% Initialize the Population
pop = initPopulation(pop_size,n);

% Run the GA
epsilon = 0.0000001;
cnEpsilon = 0;
gain = Inf;
tmp_pop = zeros(4,n);
new_pop = zeros(pop_size,n);
global_min = Inf;
total_dist = zeros(1,pop_size);
dist_history = zeros(1,num_iter);
if show_prog
    pfig = figure('Name','Current Best Solution','Numbertitle','off');
    count_subplot =1;
end

for iter = 1:num_iter
    % Evaluate Each Population Member (Calculate Total Distance)
    for p = 1:pop_size        
        tau = [0 pop(p,:) 0];
        total_dist(p) = backwardExpectedDistance(instance, tau, 0, instance.Q);
    end

    % Find the Best Route in the Population
    [min_dist,index] = min(total_dist);
    dist_history(iter) = min_dist;
    gain = (global_min-min_dist)/global_min;
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

    % Genetic Algorithm Operators
    rand_pair = randperm(pop_size);
    for p = 4:4:pop_size
        rtes = pop(rand_pair(p-3:p),:);
        dists = total_dist(rand_pair(p-3:p));
        [ignore,idx] = min(dists);
        best_of_4_rte = rtes(idx,:);
        ins_pts = sort(ceil(n*rand(1,2)));
        I = ins_pts(1);
        J = ins_pts(2);
        for k = 1:4 % Mutate the Best to get Three New Routes
            tmp_pop(k,:) = best_of_4_rte;
            switch k
                case 2 % Flip
                    tmp_pop(k,I:J) = fliplr(tmp_pop(k,I:J));
                case 3 % Swap
                    tmp_pop(k,[I J]) = tmp_pop(k,[J I]);
                case 4 % Slide
                    tmp_pop(k,I:J) = tmp_pop(k,[I+1:J I]);
                otherwise % Do Nothing
            end
        end
        new_pop(p-3:p,:) = tmp_pop;
    end
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