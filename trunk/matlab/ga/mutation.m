function [ tmp_pop ] = mutation( type, tour )
%MUTATION Randomly alter a given tour
%   type 2: Flip a subsequence of tour
%   type 3: Swap two elements of tour
%   type 4: Slide/Shift a subsequence elements of tour
%   In other case, mutation return the same input tour
ins_pts = sort(ceil(length(tour)*rand(1,2)));
I = ins_pts(1);
J = ins_pts(2);

tmp_pop = tour;
switch type
    case 2 % Flip
        tmp_pop(I:J) = fliplr(tour(I:J));
    case 3 % Swap
        tmp_pop([I J]) = tour([J I]);
    case 4 % Slide - Shift
        tmp_pop(I:J) = tour([I+1:J I]);
    otherwise % Do Nothing
end

end

