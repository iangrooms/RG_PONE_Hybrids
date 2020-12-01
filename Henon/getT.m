function T = getT(w,X)
% A highly-inefficient way to solve the ETPF optimal transport problem
% Runs out of memory at very modest ensemble sizes
M = length(w);
f = zeros(M);
for jj=1:M
    for ii=jj+1:M
        f(ii,jj) = norm(X(:,ii)-X(:,jj),2).^2;
        f(jj,ii) = f(ii,jj);
    end
end
f = f(:);
A = -eye(M^2);
b = zeros(M^2,1);
Aeq = zeros(2*M,M^2);
beq = ones(2*M,1);
beq(M+1:end) = M*w(:);
for ii=1:M
    Aeq(ii,M*(ii-1) + (1:M)) = 1;
    Aeq(ii+M,ii-1+(1:M:M^2)) = 1;
end
T = linprog(f,A,b,Aeq,beq);
T = reshape(T,[M,M]);