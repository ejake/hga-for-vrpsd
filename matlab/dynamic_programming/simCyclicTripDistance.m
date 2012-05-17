function [ cost tour ] = simCyclicTripDistance( policy, instance, typeSim )
%simTripDistance asess the total distance of a cyclic tour given for a policy
%   Input:
%       policy
%           seq: (n) Tour proposed
%           replanishments
%       instance
%           d: (n+1Xn+1)distance matrix
%           D: (n) Demand
%           Q: (Scalar) Initial capacity
%        typeSim
%               Type of simulation, i.e. 1 with earlier replanishments of the
%               policy, 0 without replanishments of policy
%   Output:
%       cost: (Scalar) total cost of the final tour
%       tour: (n+?) real tour with depot returns (? depot returns)
typeSim = logical(typeSim);
l = 0;
cost = 0;
q = instance.Q;
tour = l;
sNu = [1:instance.n]; % Set of nodes that still need to be visited (demand > 0) 
i=1;

while ~isempty(sNu)
    if instance.Cust(policy(i).m).D ~= 0
        %Earlier replanishment
        if typeSim && policy(i).a == 1
            q = instance.Q;
            if l ~= 0 % if current position l not is the depot
                cost = cost + instance.d(l+1,0+1); %return to depot
                l = 0;
            end
        end
        % Assign demand value
        if instance.Cust(policy(i).m).D == -99
            instance.Cust(policy(i).m) = instance.Cust(policy(i).m).getDemand;%Review
        end
        if q >= instance.Cust(policy(i).m).D % Capacity of vehicle >= Demand of customer
            cost = cost + instance.d(l+1, policy(i).m+1);
            q = q - instance.Cust(policy(i).m).D;
            instance.Cust(policy(i).m).D = 0;
            tour = [tour policy(i).m];
            sNu(sNu == policy(i).m) = []; %Remove node of sNu
            l = policy(i).m;
        else % Capacity of vehicle < Demand of customer (tour fail)
            cost = cost + instance.d(l+1, policy(i).m+1);
            l = policy(i).m;
            tour = [tour policy(i).m];
            auxD = instance.Cust(policy(i).m).D;

            instance.Cust(policy(i).m).D = instance.Cust(policy(i).m).D - q;

            q = q - (auxD - instance.Cust(policy(i).m).D);

            if q == 0  %return to depot
                cost = cost + instance.d(l+1,0+1);
                l=0;
                tour = [tour l];
                q = instance.Q;
            end
            l = 0;
        end        
    end
    %cyclic tour
    if i == length(policy)
        i = 1;
    else
        i = i+1;
    end
end

%Final transitio to the depot
cost = cost + instance.d(l+1,0+1);
tour = [tour 0];

end

