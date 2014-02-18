function E = backwardEd( l, ql, instance, J, tour )
%BACKWARDED Summary of this function goes here
%   e.g. tau = [0 1 2 ... n-1 n]
% (17/2/2014)

el = 0;
if(l==instance.n)
    E = instance.d(tour(l+1) + 1,1); %d(n,0)
elseif (l == 0)    
    for di = instance.Cust(tour(l+2)).PD(1):instance.Cust(tour(l+2)).PD(2)
        gamma1 = instance.d(1, tour(l+2) + 1);%d(0,l+1)
        gamma1 = (gamma1 + J(l+1, ql - di + 1))*...
            probDemand (di, tour(l+2), instance);%1/(b-a+1)
        el = el + gamma1;
    end
    E = el;
else    
    for di = instance.Cust(tour(l+2)).PD(1):instance.Cust(tour(l+2)).PD(2)
        if(ql == 0)
            gamma2 = (instance.d(tour(l+1) + 1, 1) + instance.d(tour(l+2) + 1,1));%proactive replanishment
            gamma2 = (gamma2 + J(l+1, instance.Q - di + 1))*...
                probDemand (di, tour(l+2), instance);%1/(b-a+1)
            el = el + gamma2;
        else
            if(di <= ql) % demand of customer is lesser than current vehicle capacity
                gamma1 = instance.d(tour(l+1) + 1, tour(l+2) + 1);
                gamma1 = (gamma1 + J(l+1, ql - di + 1))*...
                    probDemand (di, tour(l+2), instance);%1/(b-a+1)
                el = el + gamma1;
            else % vehicle reach its capacity without serve the customer totally
                gamma0 = (instance.d(tour(l+1) + 1, tour(l+2) +1) + 2*instance.d(tour(l+2) +1,1));
                gamma0 = (gamma0 + J(l+1, instance.Q + ql - di + 1))*...
                    probDemand (di, tour(l+2), instance);%1/(b-a+1)
                el = el + gamma0;
            end        
        end        
    end
    E = el;
end

end