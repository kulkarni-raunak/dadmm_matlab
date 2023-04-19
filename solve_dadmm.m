function [x_sol, history] = solve_dadmm(H, h, rho)


t_start = tic;
global z x_star

Output_to_screen= 0;
MAX_ITER = 1000;

n = size(H,1);
s = size(h,2);
x = zeros(n,s);

y = zeros(n,s);
sigma_H = zeros(n,n);
for i = 1:n:(n*s)
    sigma_H = sigma_H + H(:,i:(i+19));
end

if ~Output_to_screen
    fprintf('%3s\t%10s\t%10s\t%10s\n', 'iter','||xk-zk||2','||xstar-xk||2', 'objective');
end

for k = 1:MAX_ITER
    col =1;
    for i = 1:1:s
        
        % x -update step
        x(:,i) = (H(:,col:(col+19)) + rho*eye(n))\(rho*(z - y (:,i))-h(:,i));
        
        % global update for z
        z = sum(x+((1/rho)*y),2)/s;
        
        % y - update step
        y(:,i) = y(:,i) + rho*(x(:,i) - z);
        
        % value of objective function
        history.objvali(k,i)  = cost_F(H(:,col:(col+19)), h(:,i), x(:,i));
        
        col = col + n; % to access element inside H
    end
    
    % x - average
    x_sol = sum(x,2)/s;
    
    
    history.objval(k) = cost_F(sigma_H,sum(h,2),x_sol);
    history.xk_zk_norm(k) = norm(x_sol -z);
    history.xstar_xk_norm(k) = norm(x_sol - x_star);



    if ~Output_to_screen
        fprintf('%3d\t%10.4f\t%10.4f\t%10.2f\n', k,history.xk_zk_norm(k),history.xstar_xk_norm(k), history.objval(k));
    end


    if history.xk_zk_norm(k) < 1e-6
        break;
    end
end

if ~Output_to_screen
    toc(t_start);
end

end


