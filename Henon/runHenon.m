Nt = 1000;
XT = [-4;.6];
obsErr = [1;.1];
Ne = 100;
ESS0 = 10;

CRPST = zeros(2,Nt);
RTT = CRPST;
CRPSS = zeros(2,Nt);
RTS = CRPST;
CRPSF = zeros(2,Nt);
RTF = CRPST;
XMT = zeros(2,Nt);
XMS = XMT;
XMF = XMT;
ST = XMT;
SS = XMT;
SF = XMT;
alphas = zeros(Nt,1);
alphasFK = alphas;
CRPSK = zeros(2,Nt);
RTK = CRPSK;
CRPSP = zeros(2,Nt);
RTP = CRPSK;
XMK = zeros(2,Nt);
XMP = XMK;
SK = XMK;
SP = XMK;
ESS = zeros(Nt,1);

y = XT + diag(obsErr)*randn(2,Nt);
for ii=1:Nt
    X = Henon(randn(2,Ne));
    d = bsxfun(@plus,y(:,ii),-X);
    
%% ETPF Hybrid
    alpha = fsolve(@(alpha) abs(ESS0-getESS(alpha,d)),0.1);
    alpha = .5*(tanh(alpha)+1);
    alphas(ii) = alpha;
    w = exp(-(alpha/2)*(d(1,:).^2 + 100*d(2,:).^2));
    w = w/sum(w);
    T = getT(w,X);
    X1T = X*T;
    X2T = X1T;
    for jj=1:2
        indX = jj;
        % Ensemble perturbation matrix
        A = (X2T-mean(X2T,2))/sqrt(Ne-1);
        % Obs ensemble
        yn = X2T(indX,:);
        % obs perturbation matrix
        V = (yn - mean(yn))/sqrt(Ne-1);
        s2 = V*(V'); g2 = obsErr(jj)^2/(1-alpha);
        WHB = 1/(s2 + g2 + sqrt(s2*g2 + g2^2));
        X2T = X2T + (A*(V'))*(y(jj,ii)-mean(yn))*(1/(s2 + g2));
        X2T = X2T - WHB*sqrt(Ne-1)*A*(V')*V;
    end
    % Get mean & spread
    XMT(:,ii) = mean(X2T,2);
    ST(:,ii) = std(X2T,0,2);
    
    % Compute rank & CRPS at each point
    for jj=1:2
        RTT(jj,ii) = sum( X2T(jj,:)<XT(jj) ) +1;
        CRPST(jj,ii) = getCRPS(X2T(jj,:),ones(Ne,1)/Ne,XT(jj));
    end

%% MPRR Hybrid
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
    X1S = X(:,rInd);
    X2S = X1S;
    for jj=1:2
        indX = jj;
        % Ensemble perturbation matrix
        A = bsxfun(@plus,X2S,-mean(X2S,2))/sqrt(Ne-1);
        % Obs ensemble
        yn = X2S(indX,:);
        % obs perturbation matrix
        V = (yn - mean(yn))/sqrt(Ne-1);
        s2 = V*(V'); g2 = obsErr(jj)^2/(1-alpha);
        WHB = 1/(s2 + g2 + sqrt(s2*g2 + g2^2));
        X2S = bsxfun(@plus,X2S,(A*(V'))*(y(jj,ii)-mean(yn))*(1/(s2 + g2)));
        X2S = X2S - WHB*sqrt(Ne-1)*A*(V')*V;
    end
    TH = getMPRR(Ne);
    X2S = X2S*TH;
    % Get mean & spread
    XMS(:,ii) = mean(X2S,2);
    SS(:,ii) = std(X2S,0,2);

    % Compute rank and CRPS at each point
    for jj=1:2
        RTS(jj,ii) = sum( X2S(jj,:)<XT(jj) ) +1;
        CRPSS(jj,ii) = getCRPS(X2S(jj,:),ones(Ne,1)/Ne,XT(jj));
    end

%% FK13 Hybrid
    alpha = fsolve(@(alpha) abs(ESS0-getESS_FK(alpha,X)),0.1);
    alpha = (tanh(alpha)+1)/2;
    Pp = cov(X');
    R = diag([1;.01]);
    K = alpha*Pp*inv(alpha*Pp+R);
    nu = X + K*(y(:,ii)-X);
    Q = (1/alpha)*K*R*(K');
    w = mvnpdf(nu',y(:,ii)',Q+(R/(1-alpha)));
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
    X1F = nu(:,rInd) + (1/sqrt(alpha))*K*diag([1;.1])*randn(2,Ne);
    K = (1-alpha)*Q*inv((1-alpha)*Q+R);
    X2F = X1F + K*(y(:,ii) - X1F + (1/sqrt(1-alpha))*diag([1;.1])*randn(2,Ne));
    
    % Get mean & spread
    XMF(:,ii) = mean(X2F,2);
    SF(:,ii) = std(X2F,0,2);

    % Compute rank and CRPS at each point
    for jj=1:2
        RTF(jj,ii) = sum( X2F(jj,:)<XT(jj) ) +1;
        CRPSF(jj,ii) = getCRPS(X2F(jj,:),ones(Ne,1)/Ne,XT(jj));
    end
    
    %% ESRF
    X2K = X;
    for jj=1:2
        indX = jj;
        % Ensemble perturbation matrix
        A = (X2K-mean(X2K,2))/sqrt(Ne-1);
        % Obs ensemble
        yn = X2K(indX,:);
        % obs perturbation matrix
        V = (yn - mean(yn))/sqrt(Ne-1);
        s2 = V*(V'); g2 = obsErr(jj)^2;
        WHB = 1/(s2 + g2 + sqrt(s2*g2 + g2^2));
        X2K = bsxfun(@plus,X2K,(A*(V'))*(y(jj,ii)-mean(yn))*(1/(s2 + g2)));
        X2K = X2K - WHB*sqrt(Ne-1)*A*(V')*V;
    end
    % Get mean & spread
    XMK(:,ii) = mean(X2K,2);
    SK(:,ii) = std(X2K,0,2);
    
    % Compute rank & CRPS at each point
    for jj=1:2
        RTK(jj,ii) = sum( X2K(jj,:)<XT(jj) ) +1;
        CRPSK(jj,ii) = getCRPS(X2K(jj,:),ones(Ne,1)/Ne,XT(jj));
    end

    %% Pure ETPF
    w = exp(-(1/2)*(d(1,:).^2 + 100*d(2,:).^2));
    w = w/sum(w);
    ESS(ii) = 1/sum(w.^2);
    T = getT(w,X);
    X2P = X*T;
    % Get mean & spread
    XMP(:,ii) = mean(X2P,2);
    SP(:,ii) = std(X2P,0,2);
    % Compute rank & CRPS at each point
    for jj=1:2
        RTP(jj,ii) = sum( X2P(jj,:)<XT(jj) ) +1;
        CRPSP(jj,ii) = getCRPS(X2P(jj,:),ones(Ne,1)/Ne,XT(jj));
    end
    
    fprintf('Done with ii=%04d\n',ii)
end

save Henon.mat