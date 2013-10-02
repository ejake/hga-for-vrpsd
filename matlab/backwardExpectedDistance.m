function Exp = backwardExpectedDistance( tour, instance )
%BACKWARDEXPECTEDDISTANCE Dynamic programming algorithm to compute expected
%distance of a tour given
%   Assume secuential tour, i.e. tour = [0 1 2 ... n-1 n]

%Memorization
J = zeros(instance.n, instance.Q+1);%size: n X Q
%Base case
J(instance.n,:) = backwardEd(instance.n, 0, instance, J, tour);

for i = instance.n-1:-1:1
    for j = instance.Q+1:-1:1
        %fprintf('Computing gamma_%d(%d)\n',i,j-1);
        J(i,j) = backwardEd( i, j-1, instance, J, tour );
    end
end
Exp = backwardEd( 0, instance.Q, instance, J, tour );
end


function E = backwardEd( l, ql, instance, J, tour )
%BACKWARDED Summary of this function goes here
%   Assume secuential tour, i.e. tau = [0 1 2 ... n-1 n]
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
            if(di <= ql)        
                gamma1 = instance.d(tour(l+1) + 1, tour(l+2) + 1);
                gamma1 = (gamma1 + J(l+1, ql - di + 1))*...
                    probDemand (di, tour(l+2), instance);%1/(b-a+1)
                el = el + gamma1;
            else
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