%% Part (1) --- Solving using KKT conditions

s = 10; %number of subsystems
n = 20; % dimensions of x
H =zeros(n,n*s);
global x_star
[H,h] = generate_H(n,s);


sigma_H = zeros(n,n);
for i = 1:n:(n*s)
    sigma_H = sigma_H + H(:,i:(i+19));
end

sigma_h = sum(h,2);

%x = quadprog(sigma_H,sigma_h);

x_star = (sigma_H)\(-sigma_h);
OBJ_VAL = cost_F(sigma_H,sigma_h,x_star);
%% Part (2) --- Solving using D-ADMM 

% using previous variables  
global z % This will be broadcasted to all subsystems = {1,2,..,s} 
rho = 0.1:0.2:2.3;

f = figure('Name','Objective value vs Number of Iterations','NumberTitle','off');
f_1 = figure('Name','||xk - zk||2','NumberTitle','off');
f_2 = figure('Name','||xk - x*||2','NumberTitle','off');

for i=1:1:length(rho)
    z = zeros(n,1);
    
    % [x, history] = solve_admm(P, q, rho)
    [x_sol(:,i), history(i)] = solve_dadmm(H, h, rho(i));
    K(i) = length(history(i).objval);
    
    figure(f);
    subplot(4,3,i);
    plot(1:K(i), history(i).objval, 'k');
    ylabel('f(x^k)'); xlabel('iter (k)'); title(['\rho :', num2str(rho(i))]);
    
    figure(f_1);
    subplot(4,3,i);
    plot(1:K(i), history(i).xk_zk_norm, 'k');
    ylabel('norm(xk - zk)'); xlabel('iter (k)'); title(['\rho :',num2str(rho(i))]);
    
    figure(f_2);
    subplot(4,3,i);
    plot(1:K(i), history(i).xstar_xk_norm, 'k');
    ylabel('norm(xstar - xk)'); xlabel('iter (k)'); title(['\rho :',num2str(rho(i))]);
    
    
end

f_3 = figure('Name','Impact of rho on Number of Iterations required','NumberTitle','off');
plot(rho,K, 'k', 'MarkerSize', 10, 'LineWidth', 2);
    ylabel('Iterations required'); xlabel('rho'); title('Effect of changing Value of \rho');
    




