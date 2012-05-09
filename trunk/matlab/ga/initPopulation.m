function [ pop ] = initPopulation( pop_size, num_nodes )
%initPopulation Initialize the population
%   pop_size: Population size
%   num_nodes: Number of nodes
%   pop: population (pop_size X num_nodes)

    pop = zeros(pop_size,num_nodes);
    for k = 1:pop_size
        pop(k,:) = randperm(num_nodes);
    end

end

