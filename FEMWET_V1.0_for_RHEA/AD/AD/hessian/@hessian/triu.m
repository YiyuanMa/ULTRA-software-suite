function a = triu(a,k)
%TRIU         Implements  triu(a,k)  for hessians
%
%   c = triu(a,k)
%
% functionality as Matlab function triu for matrices
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  if nargin==1
    k = 0;
  end

  a.x = triu(a.x,k);
  index = ( triu( ones(size(a.x)) , k ) == 0 );
  a.dx(:,index) = 0;
  a.hx(:,index) = 0;
