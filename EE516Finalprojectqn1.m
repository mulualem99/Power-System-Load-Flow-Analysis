% EE516 Power Systems Final Project 
% Name: Mulualem Ayena Student ID: 5237753
% Problem 1/Gauss-Seidel Program
clc; clear; close all;

for case_num = 1:3

    % Line data: [from to R X]
    line = [
        1 2 0.005 0.10;   % Line 1
        5 2 0.010 0.10;   % Line 2
        1 5 0.005 0.15;   % Line 3
        5 6 0.001 0.05;   % Line 4
        5 3 0.005 0.10;   % Line 5
        3 6 0.005 0.20;   % Line 6
        6 4 0.010 0.30;   % Line 7
        6 2 0.005 0.05    % Line 8
    ];

    if case_num == 3
        line(1,:) = [];     % remove Line 1
    end
    nbus = 6;
    Ybus = zeros(nbus);

    % Build Ybus
    for k = 1:size(line,1)
        i = line(k,1);
        m = line(k,2);
        z = line(k,3) + 1j*line(k,4);
        y = 1/z;
        Ybus(i,i) = Ybus(i,i) + y;
        Ybus(m,m) = Ybus(m,m) + y;
        Ybus(i,m) = Ybus(i,m) - y;
        Ybus(m,i) = Ybus(m,i) - y;
    end

    % Scheduled bus powers
    S = zeros(nbus,1);
    S(1) =  2.0 + j*0.0;
    S(2) = -2.0 - j*0.5;
    S(3) =  1.0 + j*0.0;   % PV bus, Q unknown
    S(4) =  0.0 + j*0.0;   % swing bus
    S(5) =  0.0 + j*0.0;
    S(6) = -1.0 + j*0.0;

    if case_num == 2
        S(5) = 0.0 + j*0.5;   % inject Q = 0.5 pu at bus 5
    end

    swing = 4;
    PV = 3;
    V = ones(nbus,1);
    V(swing) = 1 + j*0;
    V(PV) = 1 + j*0;
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

            if k == PV
                Scalc = V(k)*conj(Ybus(k,:)*V);
                Qcalc = imag(Scalc);
                Suse = real(S(k)) + j*Qcalc;
                Vnew = (conj(Suse)/conj(V(k)) - sumYV)/Ybus(k,k);
                V(k) = Vnew/abs(Vnew);    
            else
                V(k) = (conj(S(k))/conj(V(k)) - sumYV)/Ybus(k,k);
            end
        end

        if max(abs(V - Vold)) < tol
            break;
        end
    end

    if case_num == 1
        fprintf('Part (a):\n');
    elseif case_num == 2
        fprintf('\n'); 
        fprintf('Part (b): Q = 0.5 pu injected at Bus 5\n');
    else
        fprintf('\n');  
        fprintf('Part (c): Line 1 Removed\n');
    end
    fprintf('\nBus Voltages (pu):\n');
    fprintf('Bus      |V|        Angle(deg)\n');

    for k = 1:nbus
        fprintf('%d      %.3f      %.3f\n', ...
            k, abs(V(k)), angle(V(k))*180/pi);
    end

    Ibus = Ybus*V;
    Sbus = V.*conj(Ibus);
    fprintf('\nComplex Power at the buses (pu):\n');
    fprintf('Bus        P          Q\n');

    for k = 1:nbus
        fprintf('%d      %.3f      %.3f\n', ...
            k, real(Sbus(k)), imag(Sbus(k)));
    end
    fprintf('\nLine Current Magnitudes (pu):\n');
    for k = 1:size(line,1)
        i = line(k,1);
        m = line(k,2);
        z = line(k,3) + j*line(k,4);
        Iline = abs((V(i)-V(m))/z);
        fprintf('Line(%d) = %.3f\n', k, Iline);
    end
end