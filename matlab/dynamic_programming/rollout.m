function [ pi expPi ] = rollout( instance, state, tau )
%ROLLOUT Algorithm
%   Neuro dynamic programming
% Input:
%   (not exclusive) Cyclic Heuristic Tour tau0 = ( 0, l, l+1, ... , n, 1, ... , n-1, 0 )
%   instance: data VRPSD instance
%   state: x_l = (l, q_l, j_1, ... , j_n)
%   tau: e.g. (l,l+1,...,n,1,...n-1)
% Output:
% Policy pi improved
% Cache for J can improve exectuion time

    pi = [];
    x = state;
    sNu = 1:instance.n; % Set of nodes that still need to be visited (demand > 0)    
    l=0; % Current customer
    
    %computing the first customer to be visited by pi
    i = 0;
    minEd = inf;
    minl = 0;
    for l=1 : instance.n
        edl = backwardExpectedDistance([0 tau], instance);
        if(edl < minEd)
            minEd = edl;
            minl = l;
            min_tau = tau;
        end
        tau = circshift(tau, [1,1]); % cyclic heuristic
    end
    pi = [pi Control(minl, 0)];
    sNu(sNu == minl) = [];
    %move to next state
    x = x.move2nextState(instance, pi(i+1));
    i = i+1;
    tau = min_tau;
    
    %computing the remaining customers to be visited by pi
    while ( ~ x.isFinalState )
        for j = 1:length(sNu) % for each node unvisited ( sNu(j) ) asses expected distance for partial tour
            J0 = cost2goBackwardJ(instance, tau, i, x.q_l, 0);
            J1 = cost2goBackwardJ(instance, tau, i, x.q_l, 1);
            [edl a] = min([J0 J1]);
            if(edl < minEd)
                minEd = edl;
                minl = l;
                min_tau = tau;
            end
            tau = [tau(1:i-1) circshift(tau(i:instance.n), [1,1])]; % cyclic heuristic applied to tau subsequence
        end
        pi = [pi Control(minl, a-1)];
        sNu(sNu == minl) = [];
        %move to next state
        x = x.move2nextState(instance, pi(i+1));
        i = i+1;
        tau = min_tau;
    end    
    expPi = minEd;
end
