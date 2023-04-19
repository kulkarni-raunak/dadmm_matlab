function [H,h] = generate_H(n,s)
%---This function will generate H(n,n*s) and h(n,s) which will be-----%
%---required and used as need. Note: This is only used to ------------%
%---generate and no calculations are done here------------------------%

H =[];
h =[];


for i=1:1:s     
    for_Hi = rand(n);
    Hi = symposdef(for_Hi);
    H = [H,Hi];
    hi = rand(n,1);
    h = [h,hi];

end


end