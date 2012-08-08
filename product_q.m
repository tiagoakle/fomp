%Applies a householder Q in product form
function qx = product_q(V,mode,x)
	%applies a Q in product form
	k = size(V,2);
	qx = x;
	if mode == 2
		for j = 1:k
			mu = V(:,j)'*qx;
			qx = qx - 2*mu*V(:,j);
		end
	elseif mode ==1
		for j = k:-1:1
			mu = V(:,j)'*qx;
			qx = qx - 2*V(:,j)*(V(:,j)'*qx);
		end
	end
end


