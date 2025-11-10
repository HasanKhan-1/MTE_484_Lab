

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This template is for the IOP control design with SPA, and is incomplete.
% You need to complete it by replacing every * with the correct code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set time step
T = 0.02764604299;
tfinal =3;
s = tf('s');
z = tf('z', T);
G_C = 1.998/(0.022*s^2 +s);
G_D = c2d(G_C, T);
G_D
%% Plant Poles and Coefficients in its Partial Fraction Decomposition
stableRealPlantPoles = [0.284609185783745];
stableComplexPlantPoles = [];
unstablePlantPoles = [1.000000000000000];
stablePlantPoles = [stableRealPlantPoles stableComplexPlantPoles];
qs = [stablePlantPoles unstablePlantPoles];
% coefficents go in order of the poles
cs = [-0.031445718629690 0.055236793894020];
n = length(qs);
nhat = length(stablePlantPoles);
nreal = length(stableRealPlantPoles);
ncomplex = length(stableComplexPlantPoles);
% verify that your plant is correct!
G = 0;
for k=1:n
   G = G + cs(k)/(z-qs(k));
end
G
%% Poles Chosen in the Simple Pole Approximation of W[z]
realWPoles = [];
complexWPoles = [
%  0.356930347157456 + 0.826196542765724i
%  0.356930347157456 - 0.826196542765724i
%  0.666437870114671 + 0.558444773703814i
%  0.666437870114671 - 0.558444773703814i
%  0.818568286012092 + 0.178734331160040i
%  0.818568286012092 - 0.178734331160040i
%  0.772183714907377 + 0.227447379480407i
%  0.772183714907377 - 0.227447379480407i
%  0.531829123858141 + 0.557815187151022i
%  0.531829123858141 - 0.557815187151022i
%  0.156352764579849 + 0.718020760847653i
%  0.156352764579849 - 0.718020760847653i
% -0.246138666287763 + 0.652239033605081i
% -0.246138666287763 - 0.652239033605081i
% -0.541821158388287 + 0.372061597484514i
% -0.541821158388287 - 0.372061597484514i
% -0.614171375570681 + 0.028169512412883i
% -0.614171375570681 - 0.028169512412883i
% -0.419743298364956 + 0.384467896550295i
% -0.419743298364956 - 0.384467896550295i
% -0.037728748040347 + 0.518243708665676i
% -0.037728748040347 - 0.518243708665676i
%  0.321582463138287 + 0.335536465085832i
%  0.321582463138287 - 0.335536465085832i
0.398371042488927 + 0.057450087083452i,
0.398371042488927 - 0.057450087083452i,
0.097390191562406 + 0.313871232493950i,
0.097390191562406 - 0.313871232493950i,
-0.213734070429800 + 0.091201683852363i,
-0.213734070429800 - 0.091201683852363i
];
ps = [realWPoles complexWPoles];
% ps = generate_poles(30,0.9,0);
mreal = sum(imag(ps) == 0);     % number of purely real poles
mcomplex = sum(imag(ps) > 0)*2; % count each conjugate pair twice
m = length(ps);
%% Calculation of alpha, beta, gamma, and gamma hat
alpha = zeros(m);
for i=1:m
   for k=1:n
       alpha(i,i) = alpha(i,i) + cs(k)/(ps(i)-qs(k));
   end
end
beta = zeros(n,m);
for i=1:m
   for k=1:n
       beta(k,i) = cs(k)/(qs(k)-ps(i));
   end
end
gamma = zeros(n-nhat,m);
for i=1:m
   for j=(nhat+1):n
       gamma(j-nhat,i) = cs(j)/(qs(j)-ps(i));
   end
end
gammaHat = zeros(n-nhat,nhat);
for k=1:nhat
   for j=(nhat+1):n
       gammaHat(j-nhat,k) = cs(j)/(qs(j)-qs(k));
   end
end
% % verify on a simple example that alpha, beta, gamma, and gammahat are correct!
% alpha
% beta
% gamma
% gammaHat
%% Determination of A and b matrices for IOP equations
A = [alpha eye(m) zeros(m,nhat);
    beta [zeros(nhat,m) eye(nhat);
          zeros(size(beta,1)-nhat,m+nhat)];
    zeros(size(gamma)) gamma gammaHat];
b = [zeros(m+size(beta,1),1);
    -cs((nhat+1):end)'];
%% Determination of step response matrices
% time horizon
K = 100;
stepry = zeros(K,m+nhat);
for k=1:K
   for i=1:m
       stepry(k,i) = -(1-ps(i)^k)/(1-ps(i));
   end
   for j=1:nhat
       stepry(k,m+j) = -(1-qs(j)^k)/(1-qs(j));
   end
end
stepru = zeros(K,m);
for k=1:K
   for i=1:m
       stepru(k,i) = (1-ps(i)^k)/(1-ps(i));
   end
end
% verify on a simple example that stepry and stepru are correct!
% stepry
% stepru
%% Determination of steady state vector
steadyState = zeros(1,m+nhat);
for i=1:m
   steadyState(i) = 1/(1-ps(i));
end
for k=1:nhat
   steadyState(m+k) = 1/(1-qs(k));
end
% verify on a simple example that steadyState is correct!
%steadyState
%% Defining the variables for the optimization
wreal = sdpvar(mreal,1,'full');
wcomplex = sdpvar(mcomplex/2,1,'full','complex');
w = wreal;
for i=1:(mcomplex/2)
   w = [w;
        wcomplex(i);
        conj(wcomplex(i))];
end
xreal = sdpvar(mreal,1,'full');
xcomplex = sdpvar(mcomplex/2,1,'full','complex');
x = xreal;
for i=1:(mcomplex/2)
   x = [x;
        xcomplex(i);
        conj(xcomplex(i))];
end
xhatreal = sdpvar(nreal,1,'full');
xhatcomplex = sdpvar(ncomplex/2,1,'full','complex');
xhat = xhatreal;
for i=1:(ncomplex/2)
   xhat = [xhat;
           xhatcomplex(i);
           conj(xhatcomplex(i))];
end
%% Defining the objective function and constraints for the optimization
Objective = 0;
% IOP constraint
Constraints = [A*[w; x; xhat] == b];
% input saturation constraint
Constraints = [Constraints, max(stepru * w) <= 6/1.4, min(stepru * w) >= -6/1.4];
% steady state constraint
Constraints = [Constraints,
               steadyState*[x; xhat]+1 == 0];
% overshoot constraint
Constraints = [Constraints, max(stepry*[x;xhat]) <= 1.05 * ...
   (-steadyState * [x;xhat])];
% settling time constraint
jhat = round(0.25/T);
Constraints = [Constraints,
              max(stepry(jhat:end,:)*[x; xhat]) <= ...
              1.02*(-steadyState*[x; xhat]),
              min(stepry(jhat:end,:)*[x; xhat]) >= ...
              0.98*(-steadyState*[x; xhat])];
%% Solving the optimization problem
% set some options for YALMIP and solver
options = sdpsettings('verbose',1,'solver','mosek');
% solve the problem
sol = optimize(Constraints,Objective,options);
% obtain the solution
wsol = value(w);
xsol = value(x);
xhatsol = value(xhat);
%% Plotting the solution
figure(1)
plot(T*(1:K),stepry*[xsol;xhatsol]);
% add green horizontal lines
yline(0.98,'g','LineWidth',1.5);
yline(1.02,'g','LineWidth',1.5);
xlabel('Time [s]');
ylabel('y[k]');
figure(2)
plot(T*(1:K),stepru*wsol);
xlabel('Time [s]');
ylabel('u[k]');
% use log scale for heat map?
log_scale_flag = 1;
% heat map
figure(3)
t = linspace(0,2*pi);
plot(cos(t),sin(t),'k--');
hold on;
if log_scale_flag
   scatter(real(ps),imag(ps),50,log(abs(wsol)),'filled');
   scatter(real(qs(1:nhat)),imag(qs(1:nhat)),50,log(abs(xhatsol)),'filled');
else
   scatter(real(ps),imag(ps),50,abs(wsol),'filled');
   scatter(real(qs(1:nhat)),imag(qs(1:nhat)),50,abs(xhatsol),'filled');
end
hold off;
colormap(jet);
colorbar;
%% Recover the transfer functions
z = tf('z',T);
% calculate W
W = 0;
for i=1:m
   W = W + wsol(i)/(z-ps(i));
end
% calculate X
X = 1;
for i=1:m
   X = X + xsol(i)/(z-ps(i));
end
for k=1:nhat
   X = X + xhatsol(k)/(z-qs(k));
end
% remove the imaginary coefficients in W
[num,den] = tfdata(W);
num{1} = real(num{1});
den{1} = real(den{1});
W = tf(num,den,T);
% remove the imaginary coefficients in X
[num,den] = tfdata(X);
num{1} = real(num{1});
den{1} = real(den{1});
X = tf(num,den,T);
% find the poles and zeros of W and X
zpk(W)
zero(W)
pole(W)
zpk(X)
zero(X)
pole(X)
D = minreal(W / X, 1e-6);  % Increase tolerance to remove near-zero poles
[numerator, denominator] = tfdata(D, 'v');
format long
%% Verify the transfer function D
[numD, denD] = tfdata(D, 'v');
fprintf('numD = [%s]\n', num2str(numD, '%.15f '));
fprintf('denD = [%s]\n', num2str(denD, '%.15f '));
%% Partial Fraction Decomposition of D
[r, p, k] = residuez(numD, denD);
disp('========================================');
disp('Partial Fraction Decomposition of D[z]:');
disp('========================================');
disp('Residues (r):');
disp(r);
disp('Poles (p):');
disp(p);
disp('Direct terms (k):');
disp(k);
%% Verify design in discrete time
j = sqrt(-1);
% compute Try and Tru by hand
Try = minreal(G_D * D / (1 + G_D * D), 1e-6);
Tru = minreal(D / (1 + G_D * D), 1e-6);
figure(4)
hold on;
plot(T*(1:K),stepry*[xsol;xhatsol],'Color',[0 0 1]);  % bluestep(Try,'r');
title('Closed-Loop Output y[k]')
hold off;
figure(5)
hold on;
plot(T*(1:K),stepru*wsol);
step(Tru,'r');
title('Control Output u[k]')
hold off;
