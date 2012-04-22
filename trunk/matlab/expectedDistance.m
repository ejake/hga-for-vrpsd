function [ E ] = expectedDistance ( instance, tau, l, qi)
%Expected distance of an priori solution (base sequence)
%Tour tau = ( 0, l, l+1, ... , n, 0 ) is followed
%The first tau must start and finish with 0, so: (0, ... ,0)
%instance.d: [n+1,n+1] distance matrix

if l == instance.n %tau(instance.n) %n is the last costumer in the tour (instance.n assume that tau is sorted, tau(instance.n+1) is the last customer at tau
	E = instance.d(tau(l+1)+1,1); %d(n,0) %Review tau(l+1) or tau(l)
else
	if l == 0
        %E = instance.d(tau(l+1)+1,1) + expectedDistance ( instance, tau, l+1, qi);
		E = expectedDistance ( instance, tau, l+1, qi);
    else
        %---
        auxE0 = instance.d(tau(l+1)+1, tau(l)+1);%d(i+1,i)
        for j=0:min( qi, instance.Cust(tau(l+1)).PD(2) )
            auxE0 = auxE0 + expectedDistance(instance,tau,l+1,qi-j)*...
                probDemand (j, tau(l+1), instance);%1/(b-a+1)
        end
        for j=qi+1:instance.Cust(tau(l+1)).PD(2)
            auxE0 = auxE0 + (2*instance.d(1,tau(l+1)+1)+ ...%d(0,i+1)
                expectedDistance(instance,tau,l+1,instance.Q+qi-j))*...
                probDemand (j, tau(l+1), instance);%1/(b-a+1)
        end
        %---
        auxE1 = instance.d(1,tau(l)+1) + instance.d(1, tau(l+1)+1);%d(0,i)+d(0,i+1)
        for j=0:instance.Cust(tau(l+1)).PD(2)
            auxE1 = auxE1 + expectedDistance(instance,tau,l+1,instance.Q-j)*...
                probDemand (j, tau(l+1), instance);%1/(b-a+1)
        end
        %---
        E = min(auxE1, auxE0);
    end	
end

end