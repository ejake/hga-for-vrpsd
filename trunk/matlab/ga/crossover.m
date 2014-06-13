function c = crossover(p1, p2, n)
	c = zeros(1,n);
	for i=1:ceil(n/2)
		c((i-1)*2+1) = p1(i);
		p2(p2 == p1(i)) = [];
        if i <= n/2
            c((i-1)*2+2) = p2(i);
            p1(p1 == p2(i)) = [];
        end		
	end
end
