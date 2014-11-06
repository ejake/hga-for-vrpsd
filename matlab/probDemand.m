function [ p ] = probDemand (k, i, instance)
%Probability of customer i take the demand value k
    if( instance.Cust(i).D == k )% Demand is knowed and it's equal to k
        p = 1;
    elseif( instance.Cust(i).D ~= -99 )% Demand is knowed and it's different of k
        p = 0;
    elseif( k >= instance.Cust(i).PD(1) && k <= instance.Cust(i).PD(2) ) % Demand is unknowed and k exist in the probability distribution
        p = 1/(instance.Cust(i).PD(2)-instance.Cust(i).PD(1)+1);    
    else
        p = 0;
    end
end

