function obj = cost_F(H, h, x)

%---Just to compute value of our objective function------%

   obj = 0.5*x'*H*x + h'*x;
   
   
end