classdef Control
    %CONTROL Decision made in each estate
    %   Detailed explanation goes here
    
    properties
        m; %customer
        a; %0: move directly; 1: proactive replenishment
    end
    
    methods
        function obj = Control ( customer, replenishment )
            obj.m = customer;            
            if replenishment == 0 || replenishment == 1
                obj.a = replenishment;
            else
                %Error: the customer 'value' not exist in iniTau
                err = MException('Control:OutOfRange', ...
                    'replenishment decision is not a valide value, a in {0,1}');
                throw(err)
            end
        end
    end
    
end

