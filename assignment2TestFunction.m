% If h becomes greater than 10^-4 the function will not run
for i = 2:4
h=10.^(-i);
n=1/h;
% x on the nodes
xN = [0:h:1];
% x on the nodes excluding the end points
xi = xN(2:n);
% x on the cell centers
xC = [(h/2):h:(1-(h/2))];
mu = 0.1;
%mu = 10;
% creating sigma points
sig = (1+xC.^2)';
Sig = spdiags(sig, [0], n, n);
% 
diffusion = 4*pi*xi.*cos(2*pi*xi)-4*pi^2*(xi.^2 + 1).*sin(2*pi*xi);
advection = mu*2*pi*cos(2*pi*xi);
q =  diffusion + advection;
Pexp = convAdvec(h,Sig,mu,q);
Ptho = sin(2*pi*xC);
i+1;

figure
hold on
plot(Pexp)
plot(Ptho)
hold off
pause(10);
end