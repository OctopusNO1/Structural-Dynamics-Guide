% 
%   Stochastic subspace identification of system matices 
%
%           [A,C,G,L0] = sto_ac(y,i);
% 
%   Inputs:
%           y: matrix of measured outputs
%           i: number of block rows in Hankel matrices 
%              (i * #outputs) is the max. order that can be estimated 
%              Typically: i = 2 * (max order)/(#outputs)
%           
%   Outputs:
%                A,C: stochastic state space system
%               G,L0: covariance model
%           
%                 x_{k+1} = A x_{k} + w_{k}
%                   y_{k} = C x_{k} + v_{k}
%         E[y_k * (y_k)'] = L0
%       E[x_k+1 * (y_k)'] = G
%                
%   Optional:
%
%           [A,C,G,L0,ss] = sto_ac(y,i,n,W);
%   
%           W:    optional weighting flag
%                    CVA: canonical variate analysis (default)
%                    PC:  principal components
%                    UPC: unweighted principal components
%           n:    the system order
%           ss:   column vector with singular values
%           
%   Copyright:
%   
%           FENG Zhouquan, Feb. 2009
%
%

function [A,C,G,L0,ss] = sto_ac(y,i,n,W)

% Check the arguments
if (nargin < 3);error('sto_ac needs at least two arguments');end
if (nargin < 4);W = 'CVA';end

% Turn the data into row vectors and check the weighting matrices
[l,ny] = size(y);if (ny < l);y = y';[l,ny] = size(y);end
if (i < 0);error('Number of block rows should be positive');end
if (l < 0);error('Need a non-empty output vector');end
if ((ny-2*i+1) < (2*l*i));error('Not enough data points');end
Wn = 0;
if (length(W) == 3) 
  if (strcmp(W, 'CVA') || strcmp(W, 'cva') || strcmp(W, 'Cva'));Wn = 1;end 
  if (strcmp(W, 'UPC')|| strcmp(W, 'upc') || strcmp(W, 'Upc'));Wn = 3;end
end    
if (length(W) == 2) 
  if (strcmp(W, 'PC') || strcmp(W, 'pc') ||strcmp(W, 'Pc') );Wn = 2;end 
end
if (Wn == 0);error('W should be CVA, PC or UPC');end
W = Wn;

% Determine the number of columns in Hankel matrices
j = ny-2*i+1;

% Compute the R factor
  Y = blkhank(y/sqrt(j),2*i,j); 		% Output block Hankel
  disp('      Computing ... R factor');
  R = triu(qr(Y'))'; 			% R factor
  R = R(1:2*i*l,1:2*i*l); 		% Truncate
  clear Y
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                  BEGIN ALGORITHM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% **************************************
%               STEP 1 
% **************************************

% First compute the orthogonal projections Ob and Obm
  
Ob  = R(l*i+1:2*l*i,1:l*i);
  
% **************************************
%               STEP 2 
% **************************************

% Compute the SVD
  disp('      Computing ... SVD');
  % Compute the matrix WOW we want to take an svd off
  % W == 1 (CVA), W == 2 (PC), W == 3 (UPC)
  if (W == 1)
    W1i = triu(qr(R(l*i+1:2*l*i,1:2*l*i)'));
    W1i = W1i(1:l*i,1:l*i)';
    WOW = W1i\Ob;
  end
  if (W == 2)
    WOW = R(l*i+1:2*l*i,1:l*i)*R(1:l*i,1:l*i)';
  end
  if (W == 3)
    WOW = Ob;
  end
  [U,S,V] = svd(WOW);
  if (W == 1);U = W1i*U;end 		% CVA
  ss = diag(S);
  clear V S WOW
  
% **************************************
%               STEP 3 
% **************************************

% Determine the order from the singular values

U1 = U(:,1:n); 				% Determine U1

% **************************************
%               STEP 4 
% **************************************

% Determine gam and gamm
gam  = U1*diag(sqrt(ss(1:n)));
gamm = U1(1:l*(i-1),:)*diag(sqrt(ss(1:n)));
% And their pseudo inverses
gam_inv  = pinv(gam);
gamm_inv = pinv(gamm);
clear gam gamm

% **************************************
%               STEP 5 
% **************************************

% Determine the states Xi and Xip
Xi  = gam_inv  * Ob;
Xip = gamm_inv * R(l*(i+1)+1:2*l*i,1:l*(i+1));
clear gamm_inv

% **************************************
%               STEP 6 
% **************************************

% Determine the state matrices A and C
disp(['      Computing ... System matrices A,C (Order ',num2str(n),')']); 
Rhs = [       Xi , zeros(n,l)  ]; 	% Right hand side
Lhs = [      Xip   ;  R(l*i+1:l*(i+1),1:l*(i+1))]; % Left hand side

% Solve least squares
sol = Lhs/Rhs;

% Extract the system matrices
A = sol(1:n,1:n);
C = sol(n+1:n+l,1:n);

% **************************************
%               STEP 7 
% **************************************

% Determine delta
disp(['      Computing ... System matrices G,L0 (Order ',num2str(n),')']); 
delta = gam_inv*(R(l*i+1:2*l*i,1:l*i)*R(1:l*i,1:l*i)');
G = delta(:,l*(i-1)+1:l*i); 		% G = last l columns


% **************************************
%               STEP 8 
% **************************************

% Determine L0
L0 = R(l*i+1:l*(i+1),:)*R(l*i+1:l*(i+1),:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                  END ALGORITHM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










