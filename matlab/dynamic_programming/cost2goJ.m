function [ J ] = cost2goJ( instance, tau, l, m, q_l, replanishment )
%replanishment: 0 directly, 1 replanishment

    if(replanishment == 0)
         %Compute J0            
        J = instance.d(tau(l)+1, m+1);
        for k=0 : q_l
            p = probDemand(k, m, instance);
            J = J + p * expectedDistance(instance, tau, l-1, q_l-k); %Review: l-1 instead of l because expectedDistance(*) take first tau position in 0
        end
        for k=q_l+1 : instance.Cust(m).PD(2)
            p = probDemand(k, m, instance);
            J = J + p * 2 * instance.d(1, m+1) + ...
                expectedDistance(instance, tau, l-1, q_l + instance.Q-k); %Review if is with or without q_l
        end                                                               %Review: l-1 instead of l because expectedDistance(*) take first tau position in 0
        %end Compute J0
    else
        %Compute J1
        J = instance.d(1, tau(l)+1) + instance.d(1, m+1);
        for k=0 : instance.Cust(m).PD(2)
            p = probDemand(k, m, instance);
            J = J + p * ...
                expectedDistance(instance, tau, l-1, instance.Q-k);%Review: l-1 instead of l because expectedDistance(*) take first tau position in 0
        end        
        %end Compute J1        
    end
end