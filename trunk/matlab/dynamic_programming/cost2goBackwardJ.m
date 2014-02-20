function [ E ] = cost2goBackwardJ( instance, tau, l, q_l, replanishment )
% replanishment: 0 directly, 1 replanishment
% instance:
% tau:
% l:
% q_l:
    replanishment = logical(replanishment);
    J = [];
    if(~ replanishment)
        %Compute J0
        E = instance.d(tau(l)+1, tau(l+1)+1);
        for k=0 : q_l
            p = probDemand(k, tau(l+1), instance);
            [Ej J] = backwardExpectedDistancePartial([0 tau], instance, l + 1, q_l-k, J);
            E = E + p * Ej;
        end
        for k=q_l+1 : instance.Cust(tau(l+1)).PD(2)
            p = probDemand(k, tau(l+1), instance);
            [Ej J] = backwardExpectedDistancePartial([0 tau], instance, l+1, q_l + instance.Q-k, J);% potential overflow in l+1
            E = E + p * (2 * instance.d(1, tau(l+1)+1) + Ej);
                
        end
    else
        %Compute J1
        E = instance.d(1, tau(l)+1) + instance.d(1, tau(l+1)+1);
        for k=0 : instance.Cust(tau(l+1)).PD(2)
            p = probDemand(k, tau(l+1), instance);
            [Ej J] = backwardExpectedDistancePartial([0 tau], instance, l+1, instance.Q-k, J); % Potential overflow
            E = E + p * Ej;
                    
        end
    end    
end
