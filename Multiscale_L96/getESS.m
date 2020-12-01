function ESS = getESS(Sd,alpha)
alpha = 10^alpha;
%w = Sd(1,:);
% for ii=1:1200
%     w(ii) = -.5*alpha*Sd(:,ii)'*Sd(:,ii);
% end
w = -.5*alpha*sum(Sd.^2);
w = exp(w - max(w));
ESS = sum(w)^2/sum(w.^2);