function [ pi ] = rollout( tau0, instance, state )
%ROLLOUT Algorithm
%   Neuro dynamic programming
%   Tour tau0 = ( 0, l, l+1, ... , n, 0 )
%   instance: data VRPSD instance
%   state: x_l = (l, q_l, j_1, ... , j_n)

    pi = [];
    x = state;
    sNu = [1:instance.n]; %Set of nodes that still need to be visited (demand > 0)
    
    %First customer
    i=1;
    minTau = tau0;
    for j=1: instance.n
        tau = permuteTauElement( tau0, i+1, j );
        Elength = expectedDistance(instance, tau, 0, instance.Q);
        if j==1
            minElength = Elength;
            minTau = tau;
        else
            if minElength > Elength
                minElength = Elength;
                minTau = tau;
            end
        end
    end   
    pi = [pi Control(minTau(2), 0)]; % l_1
    %move to next state
    [ x instance sNu ] = moveState( x, instance, pi(length(pi)), sNu);
    
    %General case (customers after first customer visited)
    %I'm here
    tau = minTau;
    while ~isempty(sNu)
        i = i+1;
        l = x(1);
        ql = x(2);        
        for j=1 : length(sNu)
            m = sNu(j);
            
            %Compute J0            
            auxJ0 = instance.d(l+1, m+1);
            for k=0 : ql
                p = probabilityDemand(instance, x, m, k);
                auxJ0 = auxJ0 + p * expectedDistance(instance, tau, l, ql-k);
            end
            for k=ql+1 : instance.Cust(m).PD(2)
                p = probabilityDemand(instance, x, m, k);
                auxJ0 = auxJ0 + p * 2 * instance.d(1, m+1) + ...
                    expectedDistance(instance, tau, l, instance.Q-k);
            end            
            if(j == 1)
                J0 = auxJ0;
                argJ0 = m;
            end
            if(auxJ0 < J0)
                J0 = auxJ0;
                argJ0 = m;
            end
            %end Compute J0
            
            %Compute J1
            auxJ1 = instance.d(1, l+1) + instance.d(1, m+1);
            for k=0 : instance.Cust(m).PD(2)
                p = probabilityDemand(instance, x, m, k);
                auxJ1 = auxJ1 + p * ...
                    expectedDistance(instance, tau, l, instance.Q-k);
            end
            if(j == 1)
                J1 = auxJ1;
                argJ1 = m;
            end
            if(auxJ1 < J1)
                J1 = auxJ1;
                argJ1 = m;
            end
            %end Compute J1
        end
        %update control and move to next state
        if(J0 < J1)
            pi = [pi Control(argJ0, 0)];
            tau = permuteTauElement( tau, i+1, argJ0 );
        else            
            pi = [pi Control(argJ1, 1)];
            tau = permuteTauElement( tau, i+1, argJ1 );
        end
        %move to next state
        [ x instance sNu ] = moveState( x, instance, pi(length(pi)), sNu);
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
    %Important: tau can be shifted instead of exchanged
    %tau can be permuted randomly
end

function [ state instance rmNodes ] =  moveState ( iniState, auxInstance, ctrl, rmNodes)
%moveState: 
%   iniState: state before updating
%   auxInstance: VRPSD instance
%   ctrl: object of the class Control, decision (m,a)
%   Output>
%   state: state after updating movement
%   instance: demand of customer visited is updated
%   rmNodes (remaining nodes): set of nodes that still need to be visited
    state = iniState;
    instance = auxInstance;
    if(auxInstance.Cust(ctrl.m).D == -99)%Review
        instance.Cust(ctrl.m) = auxInstance.Cust(ctrl.m).getDemand;%Review
        state(ctrl.m + 2) = instance.Cust(ctrl.m).D;
    end    
    q = iniState(2);
    
    if ctrl.a == 0%directly
        if(state(ctrl.m + 2) <= q)
            state(2) = q - state(ctrl.m + 2);
            state(ctrl.m + 2) = 0;
        else
            state(2) = max( [0, q - state(ctrl.m + 2)]);
            state(ctrl.m + 2) = state(ctrl.m + 2) - q;
            %I'm not including replanishment here
        end
    else%replanishment
        state(2) = instance.Q - state(ctrl.m + 2);
        state(ctrl.m + 2) = 0;
    end    
    state(1) = ctrl.m;
    
    %List of nodes that still need to be visited
    %Remove node with demand 0
    if(state(ctrl.m + 2) == 0)
        rmNodes(rmNodes == ctrl.m) = [];
    end
end

function [ p ] = probabilityDemand ( instance, currentState, customer, demand )
%Assess probability of customer take a demand value
% Input arguments>>
% instance: VRPSD instance
% currentState: state
% customer: customer whose demand will be assessed
% demand: value of demand to assess probability

    if( currentState(customer + 2) == demand )
        p = 1;
    elseif( currentState(customer + 2) ~= -99 )
        p = 0;
    else
        p = 1/instance.Cust(customer).PD(2);    
    end
end