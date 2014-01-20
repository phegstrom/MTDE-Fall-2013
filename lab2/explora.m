function [t,x] = explora(NIA,theta)

rng(NIA)

t = linspace(0,1e-4,1001);
t = t(2:end);
r = 3e8 * t/2;

idx = round(400 + 200*rand(1));
idx = idx:idx+15;

theta_range1 = 6 * rand(1);
theta_range2 = theta_range1 + 20*pi/500;

rng('shuffle')

alpha = 5e5;
v = 1e-4;

x = normrnd(0,sqrt(v),1,length(r));

if (theta>=theta_range1 & theta<=theta_range2)
    
    x(idx) = x(idx) + alpha./(r(idx).^2);
    
end
