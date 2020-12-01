Nt = 1000;
XT = [-4;.6];
obsErr = [1;.1];
Ne = 1E4;
y = XT + diag(obsErr)*randn(2,Nt);
XM = zeros(2,Nt);
CRPS = zeros(2,Nt);
alpha = 1;
for ii=1:Nt
    X = Henon(randn(2,Ne));
    d = y(:,ii)-X;
    w = exp(-(alpha/2)*(d(1,:).^2 + 100*d(2,:).^2));
    w = w/sum(w);
    % resample
    u = rand(1)/Ne + (0:Ne)/Ne;
    u(end) = 1;
    P = cumsum(w);
    jj=1;
    kk=1;
    rInd = 1:Ne;
    while (jj<=Ne)
        if (u(jj)<P(kk))
            rInd(jj)=kk;
            jj=jj+1;
        else
            kk=kk+1;
        end
    end
    X2S = X(:,rInd);
    % Get mean & spread
    XM(:,ii) = mean(X2S,2);

    % Compute CRPS at each point
    for jj=1:2
        CRPS(jj,ii) = getCRPS(X2S(jj,:),ones(Ne,1)/Ne,XT(jj));
    end
end