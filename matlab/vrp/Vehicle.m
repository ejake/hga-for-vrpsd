classdef Vehicle
    %Vehicle class define vehicle properties for VRPSD
    %   Detailed explanation goes here
    
    properties
        Q; %Max Capacity
        l = 0; %Current location (\in customers)
        ql; %Curren capacity in location l
    end
    
    methods
        function obj = Vehicle(maxCapacity, location, currentCapacity)
            obj.Q = maxCapacity;
            obj.l = location;
            obj.ql = currentCapacity;
        end
    end
    
end

