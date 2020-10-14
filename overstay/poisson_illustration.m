% overstay is modeled as poisson distribution
for lambda = 0.3:0.1:0.7
    x = 0:10;
    y = poisspdf(x,1/lambda);
    plot(x,y,'LineWidth',2)
    hold on
end

legend('y_1','y_2','y_3','y_4','y_5','FontSize',15)
xlabel('Overstay duration [hour(s)]','FontSize',15)
% 
% for lambda = 0.3:0.1:0.7
% %     plot(1/lambda, exp(-1/lambda)*(1/lambda)^(1/lambda)/factorial(round(1/lambda)),'o','MarkerSize',12)
%     plot(round(1/lambda), y(round(1/lambda)),'o','MarkerSize',12)
%     hold on
% end
title('Probability distribution of overstay w.r.t penalty $y$', 'FontSize',15, 'interpreter','latex')
hold off