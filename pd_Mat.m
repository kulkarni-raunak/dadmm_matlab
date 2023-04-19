% function to convert a matrix to positive definite matrix
function Aspd= pd_Mat(A)
if nargin==0 
    disp('No argument provided');
    disp('Please provide a symmetric matrix... (2x2 or 3x3 ...) as an argument');
    return;
end
% Making sure matrix A is symmetric 
size_A = size(A);
isSym=@(x) isequal(x,x.'); % Check weather the matrix is symmetric
sn=isSym(A);
if size_A(1) ~= size_A(2)  || max(size_A) < 2 || sn == 0
    
    disp('Please provide a symmetric matrix... (2x2 or 3x3 ...)');
    Aspd = [];
    return;
    
end
[R,p]=chol(A);
if p > 0
    [Vec,Val] = eig(A);
    Val(Val==0) = 1; % making zero eigenvalues non-zero
    size_Val = size(Val);
    for  i = 1:1:size_Val(1)
        for j = 1:1:size_Val(2)
            if i ~= j
                Val(i,j) = 0;
            end
        end
    end
            
    nidx = find(Val<0);
    Val(nidx) = -Val(nidx); % making negative eigenvalues positive
    
    Aspd = int(Vec * Val * Vec');
    
else
    Aspd = A;
    disp('Provided matrix is aleready positive Definite');
end
return