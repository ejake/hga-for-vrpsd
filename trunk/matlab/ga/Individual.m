classdef Individual
    %Customer Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        tour = []; 
        policy = [];
        rolledout = 0;
    end
    
    methods
        function obj = Instance(n)
           obj.tour = zeros(1,n);
        end        
    end    
end

