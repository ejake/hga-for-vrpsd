function varargout = vrpsd_ga(instance, pop_size, num_iter, epsilon, m, p_m, alpha, show_prog,show_res)
%VRPSD_GA Vehilce Routing Problem Whit Stochastic Demands (VRPSD) Genetic Algorithm (GA)
%   Finds a (near) optimal solution to the VRPSD by setting up a GA to search
%   for the shortest route (least expected distance for the vehicle to travel to
%   each customer and return to the depot)
%
% Input:
%     INSTANCE (object instance)     
%     POP_SIZE (scalar integer) is the size of the population
%     NUM_ITER (scalar integer) is the number of desired iterations for the algorithm to run
%     epsilon (fraction)
%     m (fraction) percentage (or number) of iterations without change
%     alpha (fraction)
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
    show_prog = logical(show_prog(1));
    show_res = logical(show_res(1));

    % Initialize the Population
    pop = initPopulationCyclic(pop_size, n); %if pop_size >= n, pop_size - n individuals would be repeated

    if epsilon <= 0
        epsilon = 0.001;
    end
    cnEpsilon = 0;

    local_search = false;
    
    gain = Inf;
    gain_relative = 0; %|f(k)-f(k-1)\f(k-1)|
    rate_change = 0;
    m_change = 0; %number of generation whiout significative change in fitness function
    offspring_pop = Individual.empty(pop_size,0); % the size of offspring is almost the size of population
    
    
    global_min = Inf;
    total_dist = zeros(1,pop_size);
    dist_history = zeros(1,num_iter);

    if show_prog
        pfig = figure('Name','Current Best Solution','Numbertitle','off');
        count_subplot =1;
    end

    if local_search
        %Apply RA to a member of population
        [pi cyEd]  = rollout (instance, State(instance.n, instance.Q), pop(1).tour);
        %if cyclic heuristic is used  in RA expected distance is already computed
        for i=1:n
            pop(i).expected_distance = cyEd(i);
            pop(i).rolledout = true;
        end
        %Add to offspring policy rolled out
        offspring_counter = offspring_counter+1;
        offspring_pop(offspring_counter) = Individual();
        offspring_pop(offspring_counter).policy = pi;
        offspring_pop(offspring_counter) = offspring_pop(offspring_counter).setTourOfPolicy();        
    end
       

    iter = 0;
    while ((iter < num_iter) && m >= m_change) % stopping criterion
        offspring_counter = 0;
        offspring_dist = zeros(1,pop_size);% expected distance of offspring
        iter = iter + 1;
        fprintf('Iteracion %i\n', iter);
        % Evaluate Each Population Member (Calculate Expected Distance)
        for p = 1:pop_size
            if pop(p).expected_distance == Inf
                pop(p).expected_distance = backwardExpectedDistance([0 pop(p).tour], instance);            
            end
            total_dist(p) = pop(p).expected_distance;
        end

        % Find the Best Route in the Population
        [min_dist,index] = min(total_dist);
        dist_history(iter) = min_dist;
        
        %fitness gain regarded to global min
        if global_min == Inf
            gain = min_dist;
        else
            gain = max((global_min-min_dist)/global_min, 0);
        end
        if iter > 1
            gain_relative = max((dist_history(iter-1)-min_dist)/dist_history(iter-1),0);
            rate_change = abs((dist_history(iter-1)-min_dist)/dist_history(iter-1));
            %no significative change
            if rate_change < epsilon
                m_change = m_change + 1;
            end
        end
        
        %upgrade global min
        if min_dist < global_min
            global_min = min_dist;
            opt_rte = pop(index);
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
        if (rand() <= p_m)
            for p = 1:ceil(randi(pop_size-1,1)*rand())% # individuals to mutate                        
                offspring_counter = offspring_counter+1;
                offspring_pop(offspring_counter) = Individual();
                offspring_pop(offspring_counter).tour = mutation(randi([1 4],1), ...
                    pop(rand_pair(p)).tour);            
            end        
        end    
        % Genetic Algorithm Operators - Crossover
        rand_pair = randperm(pop_size);%random selection of individuals to crossover
        %# of crossovers
        if local_search
            n_c = ceil(pop_size/4 * rand() ); %maximun 1/4 population
        else
            n_c = floor(pop_size * rand() ); %maximun 1/2 population
        end
        
        for p = 1: n_c
            rnd_limit = randi(pop_size-1);
            dists = total_dist( rand_pair( 1:rnd_limit ) );
            %parent (a)
            [ignore,idx] = min(dists);
            idx_pa = rand_pair(idx); %index in population, i.e. pop(idx_pa)

            rnd_limit = 1 + randi(pop_size - 1);
            dists = total_dist( rand_pair( rnd_limit : pop_size ) );
            %parent (b)
            [ignore,idx] = min(dists);
            idx_pb = rand_pair(rnd_limit - 1 + idx); %index in population, i.e. pop(idx_pb)
            offspring_counter = offspring_counter+1;
            offspring_pop(offspring_counter) = Individual();
            offspring_pop(offspring_counter).tour = crossover(pop(idx_pa).tour, pop(idx_pb).tour, n);
        end
        
        %asses fitness of offspring
        for p = 1:offspring_counter
            if offspring_pop(p).expected_distance == Inf
                offspring_pop(p).expected_distance = backwardExpectedDistance([0 offspring_pop(p).tour], instance);
            end
            offspring_dist(p) = offspring_pop(p).expected_distance;
        end

        % Local search
        if local_search
            % Rollout best tour in offspring
            [ignore,idx] = min(offspring_dist(1:offspring_counter));
            if offspring_pop(idx).rolledout == false
                [pi cyEd]  = rollout (instance, State(instance.n, instance.Q), offspring_pop(idx).tour);
                offspring_counter = offspring_counter + 1;
                offspring_pop(offspring_counter) = Individual();
                offspring_pop(offspring_counter).policy = pi;
                offspring_pop(offspring_counter) = offspring_pop(offspring_counter).setTourOfPolicy();
            end
        end

        %selection
        %Building new population

        %compute size of new population
        new_pop_size = pop_size;%pending
        if new_pop_size > pop_size
            %resize pop and offspring
        end
        pop_size = new_pop_size;
        %new population:
        for p = 1: offspring_counter
            pop(p) = offspring_pop(p);
        end
        %complete new population with cyEd of offspring
        %Review: offspring_pop(idx).tour is reapeated in this process
        if local_search
            i = 1;
            tau = offspring_pop(idx).tour;
            for p = offspring_counter+1: min(pop_size,instance.n)
                pop(p) = Individual();
                pop(p).expected_distance = cyEd(i);
                pop(p).tour = tau;
                tau = circshift(tau, [1,1]);
                i = i+1;
            end
        else
            %best tour in offspring
            [ignore,idx] = min(offspring_dist(1:offspring_counter));
            i = 1;
            tau = offspring_pop(idx).tour;
            for p = offspring_counter+1: min(pop_size,instance.n)
                pop(p) = Individual();                
                pop(p).tour = tau;
                tau = circshift(tau, [1,1]);
                i = i+1;
            end
        end

        if(gain < epsilon)
            cnEpsilon = cnEpsilon + 1;
            if(cnEpsilon > num_iter*0.1)
                showResults(show_res,instance, opt_rte, min_dist, num_iter, pop, dist_history);
                
            end
        end 
    end
        showResults(show_res,instance, opt_rte, min_dist, num_iter, pop, dist_history);
    % Return Outputs
    if nargout
        varargout{1} = global_min;
        varargout{2} = min_dist;
    end
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
