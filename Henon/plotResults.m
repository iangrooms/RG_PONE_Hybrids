load Results.mat
ESS0 = [0:10:90 95 100];
subplot(2,2,1)
plot(ESS0,RMSEF(1,:),ESS0,RMSET(1,:),ESS0,RMSES(1,:))
hold on
plot([0 100],.85*[1 1],'--k')
xlabel('ESS'),ylabel('RMSE'),title('RMSE For u')
subplot(2,2,2)
plot(ESS0,RMSEF(2,:),ESS0,RMSET(2,:),ESS0,RMSES(2,:))
hold on
plot([0 100],.077*[1 1],'--k')
xlabel('ESS'),ylabel('RMSE'),title('RMSE For v')
subplot(2,2,3)
plot(ESS0,CMF(1,:),ESS0,CMT(1,:),ESS0,CMS(1,:))
hold on
plot([0 100],.36*[1 1],'--k')
xlabel('ESS'),ylabel('Median CRPS'),title('Median CRPS For u')
subplot(2,2,4)
plot(ESS0,CMF(2,:),ESS0,CMT(2,:),ESS0,CMS(2,:))
hold on
plot([0 100],.034*[1 1],'--k')
xlabel('ESS'),ylabel('Median CRPS'),title('Median CRPS For v')