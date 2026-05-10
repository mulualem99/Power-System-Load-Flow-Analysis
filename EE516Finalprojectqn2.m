% EE516 Power Systems Final Project 
% Name: Mulualem Ayena Student ID: 5237753
% Problem 2/Example 7 using MIT formulation
clc; clear; close all;
% Line impedances: z12, z13, z23
Z = [0.02 + j*0.04;
     0.01 + j*0.03;
     0.0125 + j*0.025];

line = [
    1 2;
    1 3;
    2 3
];

% MIT formulation
YL = diag(1./Z);
% Bus incidence matrix A
A = [ 1  1  0;
     -1  0  1;
      0 -1 -1];

Ybus = A*YL*A.';

% Bus data
nbus = 3;
swing = 1;

% Power injections on 100 MVA base
S = zeros(nbus,1);
S(2) = -2.566 - j*1.102;
S(3) = -1.386 - j*0.452;

% Initial voltages
V = ones(nbus,1);
V(1) = 1.05 + j*0;   % Slack bus
tol = 1e-10;
max_iter = 1000;

% Gauss-Seidel iteration
for iter = 1:max_iter
    Vold = V;

    for k = 1:nbus
        if k == swing
            continue;
        end

        sumYV = Ybus(k,:)*V - Ybus(k,k)*V(k);
        V(k) = (conj(S(k))/conj(V(k)) - sumYV)/Ybus(k,k);
    end

    if max(abs(V - Vold)) < tol
        break;
    end
end

fprintf('\nProblem 2 - Example 7 Results\n');
fprintf('Converged in %d iterations\n', iter);
fprintf('\nBus Voltages(pu):\n');
fprintf('Bus      |V|        Angle(deg)\n');

for k = 1:nbus
    fprintf('%d      %.5f      %.4f\n', k, abs(V(k)), angle(V(k))*180/pi);
end

% Bus power injections
Ibus = Ybus*V;
Sbus = V.*conj(Ibus);

fprintf('\nBus Power Injections(pu):\n');
fprintf('Bus        P          Q\n');

for k = 1:nbus
    fprintf('%d      %.3f      %.3f\n', k, real(Sbus(k)), imag(Sbus(k)));
end

% Line flows and losses
fprintf('\nLine Flows and Losses (MW/MVAR)\n\n');
Sbase = 100;
sgn = @(x) char('+'*(x>=0) + '-'*(x<0));

for k = 1:length(Z)
    i = line(k,1); 
    m = line(k,2);

    Sij = V(i)*conj((V(i)-V(m))/Z(k))*Sbase;
    Sji = V(m)*conj((V(m)-V(i))/Z(k))*Sbase;
    SL  = Sij + Sji;
    fprintf('S%d%d = %.2f %c j%.2f\n', i, m, real(Sij), sgn(imag(Sij)), abs(imag(Sij)));
    fprintf('S%d%d = %.2f %c j%.2f\n', m, i, real(Sji), sgn(imag(Sji)), abs(imag(Sji)));
    fprintf('SL%d%d = %.2f %c j%.2f\n', i, m, real(SL), sgn(imag(SL)), abs(imag(SL)));
end