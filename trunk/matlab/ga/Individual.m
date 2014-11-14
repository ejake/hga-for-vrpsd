classdef Individual
    %Customer Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        tour = []; % no include depot
        policy = [];
        rolledout = logical(0);
        expected_distance = Inf;
        selection_probability = 1;
        operator = 0; %0 undifined, 1 mutation, 2 crossover, 3 cyclic heuristic, 4 rolledout
    end
    
    methods
        function obj = Individual(n)
            if nargin > 0
                obj.tour = zeros(1,n);
            end
        end
        
        function obj = setTourOfPolicy( obj )
            obj.tour = zeros(1,length(obj.policy));
            for i=1:length(obj.policy)
                obj.tour(i) = obj.policy(i).m;
            end            
        end
    end    
end

