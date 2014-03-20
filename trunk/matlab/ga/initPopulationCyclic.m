function [ pop ] = initPopulationCyclic( pop_size, num_nodes )
%initPopulation Initialize the population
%   pop_size: Population size
%   num_nodes: Number of nodes
%   pop: population (pop_size X num_nodes)

    pop = zeros(pop_size,num_nodes);    
    for k = 1:pop_size
        pop(k,:) = [k:num_nodes 1:k-1];
    end

end

