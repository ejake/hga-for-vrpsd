function [ Exp J ] = backwardExpectedDistancePartial( tour, instance, l, q_l, J )
%BACKWARDEXPECTEDDISTANCEPARTIAL Dynamic programming algorithm to compute 
%   expected distance of a sub tour given
%   J_{ij}: Expected distance of a tour (subtour) since node at the
%   position i (in the tour) with current capacity j-1
%   (2.17.2014)

    %Memorization
    if( isempty(J) )
        J = NaN(instance.n, instance.Q+1);%size: n X Q
        %Base case
        J(instance.n,:) = backwardEd(instance.n, 0, instance, J, tour);
        
        for i = instance.n-1:-1:l
            for j = instance.Q+1:-1:1
                if( isnan(J(i,j)) )
                    %fprintf('Computing gamma_%d(%d)\n',i,j-1);
                    J(i,j) = backwardEd( i, j-1, instance, J, tour );
                end
                if (i==l && j == q_l + 1 && l ~= 0)
                    break
                end
            end
        end
    else
        J(l, q_l + 1) = backwardEd( l, q_l, instance, J, tour );
    end 
    
    
    if l == 0
        Exp = backwardEd( 0, instance.Q, instance, J, tour );%if l = 1
    else
        Exp = J(l, q_l+1);
    end
    %disp(J);

end

