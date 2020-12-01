%% Setup
y = [-4;.6];
obsErr = [1;.1];
R = diag(obsErr.^2);
Ne = 100;
ESS0 = 30;
X = Henon(randn(2,Ne));
d = bsxfun(@plus,y,-X);

%% ETPF Hybrid
% Find alpha/split to get desired ESS
alpha = fsolve(@(alpha) abs(ESS0-getESS(alpha,d)),0.1);
alpha = .5*(tanh(alpha)+1);
% Find weights, perform ETPF update
w = exp(-(alpha/2)*((1/obsErr(1)^2)*d(1,:).^2 + (1/obsErr(2)^2)*d(2,:).^2));
w = w/sum(w);
T = getT(w,X);
X1T = X*T;
% ESRF serial assimilation
X2T = X1T;
for jj=1:2
    indX = jj;
    % Ensemble perturbation matrix
    A = bsxfun(@plus,X2T,-mean(X2T,2))/sqrt(Ne-1);
    % Obs ensemble
    yn = X2T(indX,:);
    % obs perturbation matrix
    V = (yn - mean(yn))/sqrt(Ne-1);
    s2 = V*(V'); g2 = obsErr(jj)^2/(1-alpha);
    WHB = 1/(s2 + g2 + sqrt(s2*g2 + g2^2));
    X2T = bsxfun(@plus,X2T,(A*(V'))*(y(jj)-mean(yn))*(1/(s2 + g2)));
    X2T = X2T - WHB*sqrt(Ne-1)*A*(V')*V;
end

%% MPRR Hybrid - uses same weights and split as ETPF hybrid
% resample: `systematic' resampling
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
% ESRF serial assimilation
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
    X2S = bsxfun(@plus,X2S,(A*(V'))*(y(jj)-mean(yn))*(1/(s2 + g2)));
    X2S = X2S - WHB*sqrt(Ne-1)*A*(V')*V;
end
% Get mean preserving orthogonal matrix, then resample within posterior
TH = getMPRR(Ne);
X2S = X2S*TH;

%% FK13 Hybrid
% Find alpha/split to get desired ESS
alpha = fsolve(@(alpha) abs(ESS0-getESS_FK(alpha,X)),0.1);
alpha = .5*(tanh(alpha)+1);
% Find weights and perform Gaussian mixture update
Pp = cov(X');
K = alpha*Pp*inv(alpha*Pp+R);
nu = X + K*bsxfun(@plus,y,-X);
Q = (1/alpha)*K*R*(K');
w = mvnpdf(nu',y',Q+(R/(1-alpha)));
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
X1F = nu(:,rInd) + (1/sqrt(alpha))*K*R*randn(2,Ne);
% Perform final EnKF update
K = (1-alpha)*Q*inv((1-alpha)*Q+R);
X2F = X1F + K*(y - X1F + (1/sqrt(1-alpha))*R*randn(2,Ne));

%% ESRF
X2K = X;
for jj=1:2
    indX = jj;
    % Ensemble perturbation matrix
    A = bsxfun(@plus,X2K,-mean(X2K,2))/sqrt(Ne-1);
    % Obs ensemble
    yn = X2K(indX,:);
    % obs perturbation matrix
    V = (yn - mean(yn))/sqrt(Ne-1);
    s2 = V*(V'); g2 = obsErr(jj)^2;
    WHB = 1/(s2 + g2 + sqrt(s2*g2 + g2^2));
    X2K = bsxfun(@plus,X2K,(A*(V'))*(y(jj)-mean(yn))*(1/(s2 + g2)));
    X2K = X2K - WHB*sqrt(Ne-1)*A*(V')*V;
end

%% Pure ETPF
w0 = w;
w = exp(-(1/2)*((d(1,:)/obsErr(1)).^2 + (d(2,:)/obsErr(2)).^2));
w = w/sum(w);
T = getT(w,X);
X2P = X*T;

%% Plot
clf
subplot(3,2,5)
scatter(X(1,:),X(2,:),36,'filled')
hold on
scatter(X2K(1,:),X2K(2,:),12,'k','filled')
scatter(-4,.6,48,'g','filled')
xlabel('u'),ylabel('v'),title('ESRF')
axis([-12 4 -1 1])

subplot(3,2,3)
scatter(X(1,:),X(2,:),36,'filled')
hold on
scatter(unique(X2P(1,:),'stable'),unique(X2P(2,:),'stable'),12,'k','filled')
scatter(-4,.6,48,'g','filled')
xlabel('u'),ylabel('v'),title('ETPF')
axis([-12 4 -1 1])

subplot(3,2,4)
scatter(X(1,:),X(2,:),36,'filled')
hold on
scatter(X1T(1,:),X1T(2,:),16,'filled')
scatter(X2T(1,:),X2T(2,:),12,'k','filled')
scatter(-4,.6,48,'g','filled')
xlabel('u'),ylabel('v'),title('ETPF-ESRF Hybrid')
axis([-12 4 -1 1])

subplot(3,2,6)
scatter(X(1,:),X(2,:),36,'filled')
hold on
scatter(X1S(1,:),X1S(2,:),16,'filled')
scatter(X2S(1,:),X2S(2,:),12,'k','filled')
scatter(-4,.6,48,'g','filled')
xlabel('u'),ylabel('v'),title('SIR-ESRF Hybrid')
axis([-12 4 -1 1])

subplot(3,2,2)
scatter(X(1,:),X(2,:),36,'filled')
hold on
scatter(X1F(1,:),X1F(2,:),16,'filled')
scatter(X2F(1,:),X2F(2,:),12,'k','filled')
scatter(-4,.6,48,'g','filled')
xlabel('u'),ylabel('v'),title('GMM-EnKF Hybrid')
axis([-12 4 -1 1])

[X,Y] = ndgrid(linspace(-12,4),linspace(-1,1));
p = exp(-.5*((Y/.3)^2 + (X-1+(1.4/.09)*Y.^2).^2));
subplot(3,2,1)
h = pcolor(X,Y,p);set(h,'edgecolor','none','facecolor','interp')
colormap(cmocean('haline'))
l = exp(-.5*((X-y(1)).^2 + 100*(Y-y(2)).^2));
hold on
contour(X,Y,p.*l,.1:.2:.9,'k','linewidth',2)
xlabel('u'),ylabel('v'),title('Prior (color) & Posterior (contours)')