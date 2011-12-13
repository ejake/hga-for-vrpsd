classdef Customer
    %Customer Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        D = -99; % Demand (-99 is unknowed demand)
        location = []; %customer location [x,y]
        PD = []; %Probability distribution uniform [a,b]
    end
    
    methods
        function obj = Customer(loc, prob)
            obj.location = loc;
            obj.PD = prob;
        end        
        function obj = getDemand(obj)
            if(length(obj.PD)==2)
                obj.D = randi(obj.PD,1,1);
            else
                error('uniform probability distribution wrong')
            end
        end
    end    
end

