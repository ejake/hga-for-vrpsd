function E = gammaBackwardEd( l, ql, instance )
%GAMMABACKWARDED Summary of this function goes here
%   Assume secuential tour, i.e. tau = [0 1 2 ... n-1 n]
el = 0;
if(l==instance.n)
    E = instance.d(l+1,1); %d(n,0)
elseif (l == 0)    
    for di = instance.Cust(l+1).PD(1):instance.Cust(l+1).PD(2)
        gamma1 = instance.d(l+1, l+2);%d(0,l+1)
        gamma1 = (gamma1 + gammaBackwardEd(l+1, ql - di, instance))*...
            probDemand (di, l+1, instance);%1/(b-a+1)
        el = el + gamma1;
    end
    E = el;
else    
    for di = instance.Cust(l+1).PD(1):instance.Cust(l+1).PD(2)
        if(ql == 0)
            gamma2 = (instance.d(l+1, 1) + instance.d(l+2,1));%proactive replanishment
            gamma2 = (gamma2 + gammaBackwardEd(l+1, instance.Q - di, instance))*...
                probDemand (di, l+1, instance);%1/(b-a+1)
            el = el + gamma2;
        else
            if(di <= ql)        
                gamma1 = instance.d(l+1, l+2);
                gamma1 = (gamma1 + gammaBackwardEd(l+1, ql - di, instance))*...
                    probDemand (di, l+1, instance);%1/(b-a+1)
                el = el + gamma1;
            else
                gamma0 = (instance.d(l+1, l+2) + 2*instance.d(l+2,1));
                gamma0 = (gamma0 + gammaBackwardEd(l+1, instance.Q + ql - di, instance))*...
                    probDemand (di, l+1, instance);%1/(b-a+1)
                el = el + gamma0;
            end        
        end        
    end
    E = el;
end

end

