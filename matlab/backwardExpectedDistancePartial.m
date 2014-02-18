function Exp = backwardExpectedDistancePartial( tour, instance, l, q_l )
%BACKWARDEXPECTEDDISTANCEPARTIAL Dynamic programming algorithm to compute 
%   expected distance of a sub tour given
%   (2.17.2014)

    %Memorization
    J = zeros(instance.n, instance.Q+1);%size: n X Q
    %Base case
    J(instance.n,:) = backwardEd(instance.n, 0, instance, J, tour);

    for i = instance.n-1:-1:l
        for j = instance.Q+1:-1:1
            %fprintf('Computing gamma_%d(%d)\n',i,j-1);
            J(i,j) = backwardEd( i, j-1, instance, J, tour );
            if (i==l && j == q_l && l ~= 0)
                break
            end
        end
    end
    
    if l == 0
        Exp = backwardEd( 0, instance.Q, instance, J, tour );%if l = 1
    else
        Exp = J(l, q_l);
    end
    disp(J);

end

