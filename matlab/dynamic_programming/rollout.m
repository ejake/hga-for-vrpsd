function [ pi ] = rollout( tau0, instance, state )
%ROLLOUT Algorithm
%   Neuro dynamic programming
% Input:
%   (not exclusive) Cyclic Heuristic Tour tau0 = ( 0, l, l+1, ... , n, 1, ... , n-1, 0 )
%   instance: data VRPSD instance
%   state: x_l = (l, q_l, j_1, ... , j_n)
% Output:
% Policy pi improved

    pi = [];
    x = state;
    sNu = [1:instance.n]; % Set of nodes that still need to be visited (demand > 0)    
    l=0; % Current customer
    
    %computing the first customer to be visited by pi
    minEd = inf;
    minl = 0;
    for l=1 : instance.n
        tau_l = [0 l:instance.n 1:l-1];
        edl = backwardExpectedDistance(tau_l, instance);
        if(edl < minEd)
            minEd = edl;
            minl = l;
        end
    end
    
end
