function [V,R] = householder_qr(A)
	%Householder QR algorithm of A 
	%Returns Q as v1,..vk such that Q = (I-2v1v1')...(1-2vkvk')
	[m,n] = size(A);
	md    = min(m,n);
	V     = zeros(m,min(m,n));

	for j = 1:md
		[v] = house(A(j:end,j));
		A(j:end,j:end) = A(j:end,j:end) - 2*v*(v'*A(j:end,j:end));
		if j<=m
			V(j:m,j) = v;
		end
	end
	R = triu(A);
end
