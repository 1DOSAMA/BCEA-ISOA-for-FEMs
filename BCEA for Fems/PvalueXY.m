function [p,mu_x,mu_y,sigma_x,sigma_y,rho] = PvalueXY(x,y,xL,yL)

mu_x = mean(x);
mu_y = mean(y);
sigma_x = std(x);
sigma_y = std(y);
pxy = corrcoef(x,y);
rho = pxy(1,2);
% for i=1:numel(x)
%     for j=1:numel(y)
%         %t1 = (x(i)-ux)^2/sigx^2 - 2*pxy*(x(i)-ux)*(y(j)-uy)/sigx/sigy + (y(j)-uy)^2/sigy^2;
%         %t2 = exp(-t1/(1-pxy^2)/2);
%         %p(i,j)=t2/(2*pi*sigx*sigy*(1-pxy^2)^0.5);
%         %BVN=(exp(-((((x-mu_x).^2./sigma_x.^2)-(2.*rho.*(x-mu_x).*(y-mu_y)./(sigma_x.*sigma_y))+((y-mu_y).^2)./sigma_y.^2)./(2.*(1-rho.^2))))./(2*sigma_x*sigma_y.*pi.*sqrt(1-rho.^2)));
%     end
% end
p=(exp(-((((xL-mu_x).^2./sigma_x.^2)-(2.*rho.*(xL-mu_x).*(yL-mu_y)./(sigma_x.*sigma_y))+((yL-mu_y).^2)./sigma_y.^2)./(2.*(1-rho.^2))))./(2*sigma_x*sigma_y.*pi.*sqrt(1-rho.^2)));

end