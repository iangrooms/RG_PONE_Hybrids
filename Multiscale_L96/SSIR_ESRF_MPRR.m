Ne = 400;
locRad = 279;
ESS0 = 297;
ell = 1/20;
rInf = sqrt(1.06);

params = struct('F',8,'K',41,'J',128,'h',.38);
params.N = params.J*params.K;
Nt = 1500;
No = params.N/4;
dtObs = 1.2;
dxObs = 4;
obsErr = sqrt(0.5);

%% Get reference data and obs
rng('shuffle') % ensure different initial seeds for each run
[~,Y] = ode45(@(t,y) RHS(t,y,params),[0 linspace(9,9+dtObs*(Nt-1),Nt)],...
            randn(params.N,1));
XT = Y(2:end,:)';
clear Y
Y = XT(1:dxObs:end,:); % Observations
Y = Y + obsErr*randn(size(Y)); % Add obs error

%% Initialize ensemble
parfor jj=1:Ne
    [~,sol] = ode45(@(t,y) RHS(t,y,params),[0 1 9],...
        randn(params.N,1));
    X(:,jj) = sol(3,:)';
end

% Allocate space for results
FM =  zeros(size(XT));
FS = FM;
AM = FM;
AS = FM;
FCRPS = FM;
ACRPS = FM;

% localization
x = [0:params.N/2 -(params.N/2)+1:-1]';
loc1 = exp(-.5*(x/locRad).^2);
x = linspace(0,1,params.N+1);
x = x(1:end-1);

% Smoothing setup
Sk = [0:size(Y,1)/2 -(size(Y,1)/2)+1:-1]';

% Split setup
w = ones(Ne,1)/Ne;
aSplit = zeros(Nt,1);
Nr = aSplit; % Number of unique samples left after resampling
ESS = aSplit;
% SSIR-ESRF
for ii=1:Nt
    % Save forecast mean & spread
    FM(:,ii) = mean(X,2);
    FS(:,ii) = std(X,0,2);
    % Compute forecast CRPS at each point
    for jj=1:params.N
        FCRPS(jj,ii) = getCRPS(X(jj,:),ones(Ne,1)/Ne,XT(jj,ii));
    end
    % compute standardized innovations
    d = (1/obsErr)*bsxfun(@plus,Y(:,ii), -X(1:dxObs:end,:));
    % smooth innovations
    Sd = real(ifft(bsxfun(@times,1./(1+(ell*Sk).^2).^2,fft(d))));
    % To get the unsmoothed version just set Sd = d;
    % weights
    aSplit(ii) = bisect(@(alpha) (ESS0-getESS(Sd,alpha)),-6,0);
    w = -.5*aSplit(ii)*sum(Sd.^2);
    w = exp(w - max(w));
    w = w/sum(w);
    ESS(ii) = 1/sum(w.^2);
    if( aSplit(ii) > 0 )
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
        X = X(:,rInd);
        Nr(ii) = length(unique(rInd));
    else
        Nr(ii) = Ne;
    end
    % ESRF
    % Inflation
    A = rInf*(X-mean(X,2));
    X = mean(X,2) + A;
    for jj=1:length(Y(:,ii))
        indX = 1 + dxObs*(jj-1);
        % Ensemble perturbation matrix
        A = bsxfun(@plus,-mean(X,2),X)/sqrt(Ne-1);
        % Obs ensemble
        yn = X(indX,:);
        % obs perturbation matrix
        V = (yn - mean(yn))/sqrt(Ne-1);
        s2 = V*(V'); g2 = obsErr^2/(1-aSplit(ii));
        WHB = 1/(s2 + g2 + sqrt(s2*g2 + g2^2));
        X = bsxfun(@plus,X,circshift(loc1,[indX-1 0]).*(A*(V'))*(Y(jj,ii)-mean(yn))*(1/(s2 + g2)));
        X = X - WHB*sqrt(Ne-1)*bsxfun(@times,circshift(loc1,[indX-1 0]),A*(V')*V);
    end
    % Save analysis mean & spread
    AM(:,ii) = mean(X,2);
    AS(:,ii) = std(X,0,2);
    % MPRR
    TH = getMPRR(Ne);
    X = bsxfun(@plus,AM(:,ii),bsxfun(@plus,-AM(:,ii),X)*TH);
    % Compute analysis CRPS at each point
    for jj=1:params.N
        ACRPS(jj,ii) = getCRPS(X(jj,:),ones(Ne,1)/Ne,XT(jj,ii));
    end
    % forecast ensemble
    parfor jj=1:Ne
        [~,sol] = ode45(@(t,y) RHS(t,y,params),[0 0.6 1.2],X(:,jj));
        X(:,jj) = sol(3,:)';
    end
end
save('DATA.mat','FM','FS','rInf','locRad','Ne','ACRPS','FCRPS','AM','AS','XT','ESS','Nr','aSplit','ESS0')