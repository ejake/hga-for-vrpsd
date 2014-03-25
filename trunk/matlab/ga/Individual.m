classdef Individual
    %Customer Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        tour = []; 
        policy = [];
        rolledout = logical(0);
        expected_distance = Inf;
    end
    
    methods
        function obj = Instance(n)
           obj.tour = zeros(1,n);
        end        
    end    
end

