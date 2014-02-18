function [ E ] = cost2goBackwardJ( instance, tau, l, q_l, replanishment )
% replanishment: 0 directly, 1 replanishment
% instance:
% tau:
% l:
% q_l:
    replanishment = logical(replanishment);
    if(~ replanishment)
        %Compute J0
        E = instance.d(tour(l)+1, tour(l+1)+1);
        for k=0 : q_l
            p = probDemand(k, tour(l+1), instance);
            E = E + p * backwardExpectedDistancePartial([0 tau], instance, l + 1, q_l-k);
        end
        for k=q_l+1 : instance.Cust(tour(l+1)).PD(2)
            p = probDemand(k, tour(l+1), instance);
            E = E + p * (2 * instance.d(1, m+1) + ...
                backwardExpectedDistancePartial([0 tau], instance, l+1, q_l + instance.Q-k));% potential overflow in l+1
        end
    else
        %Compute J1
        E = instance.d(1, tau(l)+1) + instance.d(1, tau(l+1)+1);
        for k=0 : instance.Cust(tau(l+1)).PD(2)
            p = probDemand(k, m, instance);            
                E = E + p * ...
                    backwardExpectedDistancePartial([0 tau], instance, l+1, instance.Q-k); % Potential overflow
        end
    end    
end
