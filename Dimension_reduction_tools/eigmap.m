function E = eigmap(X,energy,T)
E = [];
for j=1:19
    Y = X(:,51+100*(j-1):100*j);
    YE = energy(51+100*(j-1):100*j);
    [d,e] = diffmap3(Y,YE,T,0.5,0.1);
    E = [E,e];
end

E = E(2:6,:);
E = log(E);
subplot(2,3,[1,2,4,5])
for i=1:5
    plot(E(i,:),'-o','linewidth',2);
    hold on
end

legend('\lambda_1','\lambda_2','\lambda_3','\lambda_4','\lambda_5')
set(legend,'color','none');
hold off
subplot(2,3,3)
plot((E(1,:)-E(5,:))/(E(1,1)-E(5,1)),'-o','linewidth',2)
hold on
plot(sum(E,1)/sum(E(:,1)),'-o','linewidth',2)
legend('M(j)=\lambda_{max}^j-\lambda_{min}^j','M(j)=\Sigma_{i=1}^n\lambda_{i}^j')
set(legend,'color','none');
ylabel('M(j)/M(1)')
hold off
end