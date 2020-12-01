function ESS = getESS_FK(alpha,X)
% Note that the obs error is hard-coded
alpha = (tanh(alpha)+1)/2;
Pp = cov(X');
R = diag([1;.01]);
K = alpha*Pp*inv(alpha*Pp+R);
y = [-4;.6];
nu = X + K*bsxfun(@plus,y,-X);
Q = (1/alpha)*K*R*(K');
w = mvnpdf(nu',y',Q+(R/(1-alpha)));
w = w/sum(w);
ESS = 1/sum(w.^2);
