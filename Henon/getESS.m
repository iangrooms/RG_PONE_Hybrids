function ESS = getESS(alpha,d)
% Note that the obs error is hard-coded
alpha = .5*(tanh(alpha)+1);
w = exp(-(alpha/2)*(d(1,:).^2 + 100*d(2,:).^2));
w = w/sum(w);
ESS = 1/sum(w.^2);
