function [p,mu_x,mu_y,sigma_x,sigma_y,rho] = PvalueXY(x,y,xL,yL)

mu_x = mean(x);
mu_y = mean(y);
sigma_x = std(x);
sigma_y = std(y);
pxy = corrcoef(x,y);
rho = pxy(1,2);
p=(exp(-((((xL-mu_x).^2./sigma_x.^2)-(2.*rho.*(xL-mu_x).*(yL-mu_y)./(sigma_x.*sigma_y))+((yL-mu_y).^2)./sigma_y.^2)./(2.*(1-rho.^2))))./(2*sigma_x*sigma_y.*pi.*sqrt(1-rho.^2)));

end