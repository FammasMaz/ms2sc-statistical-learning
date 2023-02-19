%% Dataset
clear all; close all; clc;
%rng('shuffle')
mu_tan = 0.5774;
mu_c=100;
q=250;
cv_tan=20;
cv_c=20;
sig_tan=mu_tan*0.2

N=5000;
u1=rand(1,N);
u2=rand(1,N);
% plot(u1,'o');
% %not correlation, so we can use the same
% hold on
% plot(u2,'o');
figure(1)
subplot(1,2,1);
histogram(u2, 'Normalization', 'cdf');
subplot(1,2,2);
plot(u1,u2,'o');


%% N[0,1]
x1test=norminv(u1,0,1); % Generation of N[0,1] 
x2test=norminv(u2,0,1);%N[0,1]

subplot(1,2,1);
histogram(x1test, 'Normalization', 'pdf');

subplot(1,2,2);
plot(x1test,x2test,'o')


%% Generation of variables
x1=norminv(u1,0,1); % Generation of N[0,1] 
x2=norminv(u2,0,1); % Generation of N[0,1];

tanphi = mu_tan + cv_tan*x1;
c = mu_c + cv_c*x2;


figure(2)
subplot(2,2,1);
histogram(tanphi, 'Normalization', 'pdf');
subplot(2,2,2);
plot(tanphi,c,'ro');
subplot(2,2,4);
histogram(c, 'Normalization', 'pdf');
subplot(2,2,3);
plot(c,tanphi,'ro');
%% Correlated variables

mu = [0, 0]; % mean vector
rho = 0.9; % correlation coefficient
sigma = [1, rho; rho, 1]; % covariance matrix
x = mvnrnd(mu, sigma, 500);
x1 = x(:, 1);
x2 = x(:, 2);% generate a 10-by-2 matrix of random numbers
rho = corrcoef(x); % compute the sample correlation coefficient
rho(1, 2) % display the correlation between the first and second column
subplot(1,2,1);
histogram(x1, 'Normalization', 'cdf');

subplot(1,2,2);
plot(x1,x2,'o')

%% Model
qu = 2.*c.*(tanphi+(1+tanphi.^2).^(.5));

%% Probability of failure pf
q = 250;
g = qu - q;
pt = find(g(1,:)<=0);
Ig = zeros(N,1);
Ig(pt, 1) = 1;
dat = Ig;
for j = 1:N
    pfN(1,j) = sum(dat(1:j,1))/j;
    qut(1,j) = sum(qu(1:j))/j;
    var_pf(1,j) = 1/(j*(j-1))*sum(dat(1:j, 1).^2) - 1/(j-1)*pfN(1, j)^2;
    cv_pf(1,j) = sqrt(var_pf(1,j))./pfN(1,j);
end
%% intervalle confiance
dat(:,1 ) = qu(1, :);
for i=2:N
     muN(1,i) = mean(dat(1:i,1));
     ICN(1,i) = tinv(0.975, i-1)*std(dat(1:i,1))/sqrt(i);
     stdN(1,i) = std(dat(1:i,1));
     Denom_p(1, i) = chi2inv(0.975, i-1);
     Denom_m(1, i) = chi2inv(1-0.975, i-1);
     Numera(1,i) = (i - 1)* var(dat(1:i,1));
     ICSTDN_p(1,i) = sqrt(Numera(1,i)./Denom_p(1,i));
     ICSTDN_m(1,i) = sqrt(Numera(1,i)./Denom_m(1,i));

end

muN(1,1) = dat(1,1);
muN(1,2) = mean(dat(1:2,1));
ICN(1,1) = ICN(1,2);
stdN(1,2) = std(dat(1:2,1));
stdN(1,1) = stdN(1,2);
ICSTDN_p(1,1) = ICSTDN_p(1,2);
ICSTDN_m(1,1) = ICSTDN_m(1,2);
%% N[mu_,sig_]
%phi=(ual-u)/standat deviation
c = mu_c + sig_c.*x1;
tan = mu_tan + sig_tan.*x2;


subplot(2,2,1);
histogram(tan, 'Normalization', 'pdf');

subplot(2,2,2);
plot(c,tan,'o');

subplot(2,2,3);
plot(tan,c,'o');

subplot(2,2,4);
histogram(c, 'Normalization', 'pdf');

%% Model
qu = 2.*c.*(tan + (1 + tan.^2).^(.5)); % maximum resistance for sol
g = qu - q; % model (performance function)
plot(g,tan,'or')

%% Propability density function
%pt trouve les indices de tous les éléments de la matrice g qui sont 
% inférieurs ou égaux à zéro dans la première ligne
pt = find(g(1,:) <= 0); % pf -> Igx Pr(g(x)<0)
%% probability of failure
Ig = zeros(N,1);
dat = zeros(N,1);
Ig(pt,1) = 1; %attribu 1 dans les celules ou la valeur de la fonction est négative ou nulle (pt est un vecteur)
dat(:,1) = Ig(:,1);
%%
for j=2:N
    pfN(1,j) = sum(dat(1:j,1))/j;%probability of failure
    var_pf(1,j) = 1/(j*(j - 1))*sum(dat(1:j,1).^2) - 1/(j - 1)*pfN(1,j)^2;% Var pr
    cv_pf(1,j) = sqrt(var_pf(1,j))./pfN(1,j);% coeficient de variance
end
pfN(1,1) = dat(1,1);
var_pf(1,1) = var_pf(1,2);
cv_pf(1,1) = cv_pf(1,2);
%%
%%=== intervalle confiance
dat(:,1) = qu(1,:);
for i=2:N
    muN(1,i) = mean(dat(1:i,1));
    ICN(1,i) = tinv(0.975,i-1)*std(dat(1:i,1))/sqrt(i);
    stdN(1,i) = std(dat(1:i,1));
    Denom_p(1,i) = chi2inv(0.975,i-1);
    Denom_m(1,i) = chi2inv(0.025,i-1);
    Numera(1,i) = (i - 1)* var(dat(1:i,1));
    ICSTDN_p(1,i) = sqrt(Numera(1,i)./Denom_p(1,i));
    ICSTDN_m(1,i) = sqrt(Numera(1,i)./Denom_m(1,i));
end
muN(1,1) = dat(1,1);
muN(1,2) = mean(dat(1:2,1));
ICN(1,1) = ICN(1,2);
stdN(1,2) = std(dat(1:2,1));
stdN(1,1) = stdN(1,2);
ICSTDN_p(1,1) = ICSTDN_p(1,2);
ICSTDN_m(1,1) = ICSTDN_m(1,2);

%% - figures matrix
%X = [u1;u2];
%X = [c; tan];
%X = [c; tan; qu];
X = [c; tan; g];
figure(9)
plotmatrix(X,'or');
%%

