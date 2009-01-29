function [d, i] = getdiag(x, k)
%GETDIAG Get diagonal of a matrix.
%
%   D = GETDIAG(X) returns the main diagonal of X.
%
%   D = GETDIAG(X, K) returns the K-th diagonal of X, where K = 0 is the main
%   diagonal, K > 0 is above the main diagonal, and K < 0 is below the main
%   diagonal.
%
%   [D, I] = GETDIAG(...) also returns the indices of the elements.
%
%   See also DIAG.

%   Author:      Peter J. Acklam
%   Time-stamp:  2004-09-21 08:34:35 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   error(nargchk(1, 2, nargin));

   if ndims(x) ~= 2
      error('X must be a 2D matrix.');
   end
   if nargin < 2
      k = 0;
   else
      if ~isnumeric(k) | any(size(k) ~= 1)
         error('K must be a scalar integer.');
      end
      k = double(k);
      if k ~= fix(k)
         error('K must be an integer.');
      end
   end

   [m, n] = size(x);

   if k == 0
      % main diagonal
      i = 1 : m+1 : m*min(m,n);
   elseif k > 0
      % diagonal above main diagonal
      i = m*k+1 : m+1 : m*min(m+k, n);
   elseif k < 0
      % diagonal below main diagonal
      i = 1-k : m+1 : m*min(m+k, n);
   else
      % we should never get here
      error('Internal error. Please tell the author what you did.');
   end

   i = i.';
   d = x(i);

