function [ Elen ] = backwardExpectedDistance( instance, tau, l, ql )
J= zeros(instance.n + 1, instance.Q +1);
for i=l:instance.n
    for j=0:ql
        Elen = backExJ(instance, tau, instance.n-i, ql-j, J);
        J(instance.n-i+1,ql-j+1) = Elen;
    end
    %J = min(sts(n-l+1,:));
end

end

function [ E ] = backExJ( instance, tau, l, qi, J )
%BackwardExpectedDistance Expected distance of an priori solution (base
%sequence)
%   The expected distance is assesing using a backward method
%
if l == instance.n %tau(instance.n) %n is the last costumer in the tour
	E = instance.d(tau(l+1)+1,1); %d(n,0)
else
    if l == 0
        E = instance.d(tau(l+2)+1,1) + J(l+2,qi+1);        
    else
        %---
        auxE0 = instance.d(tau(l+2)+1, tau(l+1)+1);%d(i+1,i)
        for j=0:qi
            auxE0 = auxE0 + J(l+2,j+1)*...
                probDemand (j, tau(l+2), instance);%1/(b-a+1)
        end
        for j=qi+1:instance.Cust(tau(l+2)).PD(2)
            auxE0 = auxE0 + (2*instance.d(1,tau(l+2)+1)+ ...%d(0,i+1)
                J(l+2,instance.Q+qi-j))*...
                probDemand (j, tau(l+2), instance);%1/(b-a+1)
        end
        %---
        auxE1 = instance.d(1,tau(l+1)+1) + instance.d(1, tau(l+2)+1);%d(0,i)+d(0,i+1)
        for j=0:instance.Cust(tau(l+2)).PD(2)
            auxE1 = auxE1 + J(l+2,instance.Q-j)*...
                probDemand (j, tau(l+2), instance);%1/(b-a+1)
        end
        %---
        E = min(auxE1, auxE0);
    end
end

end

