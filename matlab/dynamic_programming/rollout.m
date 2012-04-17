function [ pi ] = rollout( tau0, instance, state )
%ROLLOUT Algorithm
%   Neuro dynamic programming
%   Tour tau0 = ( 0, l, l+1, ... , n, 0 )
%   instance: data VRPSD instance
%   state: x_l = (l, q_l, j_1, ... , j_n)

    pi = [];
    x = state;
    sNu = [1:instance.n]; % Set of nodes that still need to be visited (demand > 0)    
    l=0; % Current customer
    
    %First customer
    i=1;
    minTau = tau0;
    minA = 0;
    
    while ( ~ x.isFinalState )
        tau = minTau;
        for j=1 : length(sNu)
            if( i == 1 )% l_1
                Elength = expectedDistance(instance, tau, 0, instance.Q);
                a = 0;
            else% l_n | n > 1
                J0 = cost2goJ(instance, tau, l, sNu(j), x.q_l, 0);
                J1 = cost2goJ(instance, tau, l, sNu(j), x.q_l, 1);
                [Elength a] = min([J0 J1]);
                a = a - 1;
            end
            
            %Evaluate minimization
            if j==1                
                minElength = Elength;
                minTau = tau;
                minA = a;
            else
                if minElength > Elength
                    minElength = Elength;
                    minTau = tau;
                    minA = a;
                end
            end
            %tau = permuteTauElement( tau, sNu(j) , j);
            tau = shiftTauElement( tau );  
        end
        
        % --- General case including l = 0
        l = minTau(i+1);
        pi = [pi Control(l, minA)];
        %move to next state
        [x instance] = x.move2nextState(instance, pi(i));
        %remove minTau(2) of sNu
        if(x.r(l) == 0)
            sNu(sNu == minTau(i+1)) = [];
        end
        i = i + 1;
    end
end

function [ tau ] = permuteTauElement( iniTau, position, value )
    %Exchange            
    pos = find(iniTau == value);
    tau = iniTau;
    if any(pos) ~= 0
        tau(pos) = iniTau(position);
        tau(position) = value;
    else
        %Error: the customer 'value' not exist in iniTau
        err = MException('customerTau:OutOfRange', ...
            '(Rollout) the customer not exist in iniTau');
        throw(err)
    end            
end

%Cyclic heuristic
function [ tau ] = shiftTauElement(iniTau, pl)
%pl Position to shift remain vector pl \in [1,size(iniTau)]
    if(pl==1)
        tau = [0 iniTau(3:length(iniTau)-1) iniTau(2) 0];
    else
        if(pl>1 && pl <= length(iniTau)-2)
            tau = [0 tau2(2:pl) tau2(2+pl:length(tau2)-1) tau2(1+pl) 0];
        else
            error('error shifting in cyclic heuristic')            
        end 
    end
end
