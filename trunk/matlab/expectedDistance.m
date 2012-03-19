function [ E ] = expectedDistance ( instance, tau, l, qi)
%Expected distance of an priori solution (base sequence)
%Tour tau = ( 0, l, l+1, ... , n, 0 ) is followed
%The first tau must start and finish with 0, so: (0, ... ,0)
%instance.d: [n+1,n+1] distance matrix

if l == instance.n %n is the last costumer in the tour
	E = instance.d(tau(l+1)+1,1); %d(n,0)
else
	if l == 0
		E = expectedDistance ( instance, tau, l+1, qi);
    else
        %---
        auxE0 = instance.d(tau(l+1)+1, tau(l)+1);%d(i+1,i)
        for j=0:qi
            auxE0 = auxE0 + expectedDistance(instance,tau,l+1,qi-j)*...
                (1/probDemand (j, tau(l+1), instance));%1/(b-a+1)
        end
        for j=qi+1:instance.Cust(tau(l+1)).PD(2)
            auxE0 = auxE0 + (2*instance.d(1,tau(l+1)+1)+ ...%d(0,i+1)
                expectedDistance(instance,tau,l+1,instance.Q+qi-j))*...
                (1/probDemand (j, tau(l+1), instance));%1/(b-a+1)
        end
        %---
        auxE1 = instance.d(1,tau(l)+1) + instance.d(1, tau(l+1)+1);%d(0,i)+d(0,i+1)
        for j=0:instance.Cust(tau(l+1)).PD(2)
            auxE1 = auxE1 + expectedDistance(instance,tau,l+1,instance.Q-j)*...
                (1/probDemand (j, tau(l+1), instance));%1/(b-a+1)
        end
        %---
        if(auxE0 >= auxE1)
            E = auxE1;
        else
            E = auxE0;
        end
    end	
end

end

%Probability of customer i take the demand value k
function [ p ] = probDemand (k, i, instance)
    if( instance.Cust(i) == k )
        p = 1;
    elseif( instance.Cust(i) ~= -99 )
        p = 0;
    else
        p = 1/(instance.Cust(i).PD(2)-instance.Cust(i).PD(1)+1);    
    end
end