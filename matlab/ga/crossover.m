function c = crossover(p1, p2, n)
	c = zeros(1,n);
	for i=1:n-1
		c(i) = p1(i);
		p2(p2 == p1(i)) = [];
		c(i+1) = p2(i);
		p1(p1 == p2(i)) = [];
	end
end
