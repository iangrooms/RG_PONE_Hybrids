Ne = 400;
locRad = 209;
rInf = sqrt(1.026);

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
    [~,sol] = ode45(@(t,y) RHS(t,y,params),[0 1 9],randn(params.N,1));
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

% ESRF
for ii=1:Nt
    % Save forecast mean & spread
    FM(:,ii) = mean(X,2);
    FS(:,ii) = std(X,0,2);
    % Compute forecast CRPS at each point
    for jj=1:params.N
        FCRPS(jj,ii) = getCRPS(X(jj,:),ones(Ne,1)/Ne,XT(jj,ii));
    end
    % Inflation
    A = rInf*(X-mean(X,2));
    X = mean(X,2) + A;
    for jj=1:length(Y(:,ii))
        indX = 1 + dxObs*(jj-1);
        % Ensemble perturbation matrix
        A = (X-mean(X,2))/sqrt(Ne-1);
        % Obs ensemble
        yn = X(indX,:);
        % obs perturbation matrix
        V = (yn - mean(yn))/sqrt(Ne-1);
        s2 = V*(V'); g2 = obsErr^2;
        WHB = 1/(s2 + g2 + sqrt(s2*g2 + g2^2));
        X = X + circshift(loc1,[indX-1 0]).*(A*(V'))*(Y(jj,ii)-mean(yn))*(1/(s2 + g2));
        X = X - WHB*sqrt(Ne-1)*bsxfun(@times,circshift(loc1,[indX-1 0]),A*(V')*V);
    end
    % Save analysis mean & spread
    AM(:,ii) = mean(X,2);
    AS(:,ii) = std(X,0,2);
    % MPRR
    TH = getMPRR(Ne);
    X = AM(:,ii) + (X-AM(:,ii))*TH;
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
save('DATA.mat','FM','FS','rInf','locRad','Ne','ACRPS','FCRPS','AM','AS','XT')
