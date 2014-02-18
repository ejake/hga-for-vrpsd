function [ J ] = cost2goBackwardJ( instance, tau, l, q_l, replanishment )
% replanishment: 0 directly, 1 replanishment
% instance:
% tau:
% l:
% q_l:
    replanishment = logical(replanishment);
    
    
    if(replanishment == 0)
         %Compute J0
        E = instance.d(tour(l)+1, tour(l)+2);
        for k=0 : q_l
            p = probDemand(k, tour(l)+2, instance);
            E = E + p * backwardExpectedDistance(instance, tau, instance.n-l, q_l-k);
        end
        for k=q_l+1 : instance.Cust(m).PD(2)
            p = probDemand(k, m, instance);
            if(l~=instance.n)
                E = E + p * 2 * instance.d(1, m+1) + ...
                    backwardExpectedDistance(instance, tau, l+1, q_l + instance.Q-k); %Review if is with or without q_l
                                                                              %Review l+1 instead of l because already is taken the replanishment decision
            else
                E = E + p * 2 * instance.d(1, m+1) + ...
                    backwardExpectedDistance(instance, tau, l, q_l + instance.Q-k); %Review if is with or without q_l
            end
        end
    end
    
    %Memorization
    J = zeros(instance.n, instance.Q+1);%size: n X Q
    %Base case
    J(instance.n,:) = backwardEd(instance.n, 0, instance, J, tour);

    for i = instance.n-1:-1:l
        for j = q_l+1:-1:1
            %fprintf('Computing gamma_%d(%d)\n',i,j-1);
            J(i,j) = backwardEd( i, j-1, instance, J, tour );
        end
    end
    Exp = backwardEd( 0, instance.Q, instance, J, tour );        
end

function E = backwardEd( l, ql, instance, E, tour )

    if(replanishment == 0)
         %Compute J0            
        E = instance.d(tour(l)+1, tour(l)+2);
        for k=0 : q_l
            p = probDemand(k, m, instance);            
            E = E + p * backwardExpectedDistance(instance, tau, instance.n-l, q_l-k);
        end
        for k=q_l+1 : instance.Cust(m).PD(2)
            p = probDemand(k, m, instance);
            if(l~=instance.n)
                E = E + p * 2 * instance.d(1, m+1) + ...
                    backwardExpectedDistance(instance, tau, l+1, q_l + instance.Q-k); %Review if is with or without q_l
                                                                              %Review l+1 instead of l because already is taken the replanishment decision
            else
                E = E + p * 2 * instance.d(1, m+1) + ...
                    backwardExpectedDistance(instance, tau, l, q_l + instance.Q-k); %Review if is with or without q_l
            end
        end
        %end Compute J0
    else
        %Compute J1
        E = instance.d(1, tau(l)+1) + instance.d(1, m+1);
        for k=0 : instance.Cust(m).PD(2)
            p = probDemand(k, m, instance);
            if(l~=instance.n)
                E = E + p * ...
                    backwardExpectedDistance(instance, tau, l+1, instance.Q-k); %Review l+1 instead of l because already is taken the replanishment decision
            else
                E = E + p * ...
                    backwardExpectedDistance(instance, tau, l, instance.Q-k); %Review l+1 instead of l because already is taken the replanishment decision
            end
        end
        %end Compute J1        
    end
end