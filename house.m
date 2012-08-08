%Calculates a householder vector to reflect x into e_1
function [v] = house(x)
	%Calculates the householder vectors
	if size(x) == 1
		v = 1;
		return
	end
	sigma = x(2:end)'*x(2:end);
	v  = [1;x(2:end)];
	if sigma == 0
		bet = 0;
	else
		mu = sqrt(x(1).^2+sigma);
		if x(1) <= 0
			v(1) = x(1)-mu;
		else
			v(1) = -sigma/(x(1)+mu);
		end
		bet = 2*v(1)^2/(sigma+v(1)^2);
		v   = v/v(1)*sqrt(bet)/sqrt(2);
	end
end

