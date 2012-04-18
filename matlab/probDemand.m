function [ p ] = probDemand (k, i, instance)
%Probability of customer i take the demand value k
    if( instance.Cust(i).D == k )
        p = 1;
    elseif( instance.Cust(i).D ~= -99 )
        p = 0;
    else
        p = 1/(instance.Cust(i).PD(2)-instance.Cust(i).PD(1)+1);    
    end
end

