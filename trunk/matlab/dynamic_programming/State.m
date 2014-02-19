classdef State % < handle
    %STATE x_l = (l,q_l,r_1,...,r_n)
    %   State x_l
    
    properties
        l = 0;
        q_l = 0;
        r = [];
    end
    
    methods
        %Initial state
        function obj = State ( n, Q )
            obj.l = 0;
            obj.q_l = Q;
            obj.r = zeros(1,n)-99;
        end        
        function ifs = isFinalState( obj )
            if(obj.r == zeros(1,length(obj.r)))
                ifs = true;
            else
                ifs = false;
            end
        end
        
        %moveState:
        function obj = move2nextState( obj, instance, u)
        %   Input:
        %   obj: state before updating (this Class (State))
        %   instance: VRPSD instance before updating
        %   u: object of the class Control, decision (m,a)
        %   Output:
        %   obj: state after updating movement (this Class (State))
        %   instance: Instance after moving, demand of customer visited is updated
            %Assign demand if it is unknowed            
            
            if(instance.Cust(u.m).D == -99)
                instance.Cust(u.m) = instance.Cust(u.m).getDemand; % Review
                obj.r(u.m) = instance.Cust(u.m).D;
            end
            q = obj.q_l;
            if u.a == 0 % directly
                if(obj.r(u.m) <= q)
                    obj.q_l = q - obj.r(u.m);
                    obj.r(u.m) = 0;
                else
                    obj.q_l = max( [0, q - obj.r(u.m)]);
                    obj.r(u.m) = obj.r(u.m) - q;
                    %I'm not including replanishment here
                end
            else % replanishment
                obj.q_l = instance.Q - obj.r(u.m);
                obj.r(u.m) = 0;
            end    
            obj.l = u.m;
        end
    end    
    
end

