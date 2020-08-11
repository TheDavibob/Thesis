function [Z] = Zeros_Box(f,x,y,varargin)
% for a function f(z) which is analytic in the box [x1 x2] x [y1 y2], finds
% the number of zeros.

% N_opt is optional, allows for varying the integration precision.
% prec_opt allows definition of how close to the real line a "real" zero must
% be, again optional, and how precise the numerical derivative is.

% Z has components:
%   Z.z number of zeros
%   Z.r number of real zeros
%   Z.U number of zeros with postiive imaginary part
%   Z.L number of zeros with negative imaginary part

% x and y are two-vectors

% created 2018-02-07

N = 500;
prec = 1e-2;

if nargin == 4
	N=varargin(1) ;
elseif nargin == 5
    N=varargin{1} ; prec = varargin{2};
end

prec2 = prec*1e-5;

x1=min(x) ; x2 = max(x);
y1=min(y) ; y2 = max(y);


F.f = @(z) f(z);
F.df = @(z) (f(z+prec2) - f(z-prec2)) ./ (2*prec2) ;

F.dfdivf = @(z) F.df(z) ./ (2*pi*1i *F.f(z)) ;

[t,wt] = lgwt(N,0,1);
t=t(end:-1:1).'; wt=wt(end:-1:1).';

z.ll = x1 + 1i*y1;
z.lu = x1 + 1i*y2;
z.ul = x2 + 1i*y1;
z.uu = x2 + 1i*y2;

% Total zeros
I.xl = sum(F.dfdivf( z.ll + (z.lu - z.ll)*t ).*(z.lu-z.ll).*wt,2);
I.xu = sum(F.dfdivf( z.ul + (z.uu - z.ul)*t ).*(z.uu-z.ul).*wt,2);
I.yl = sum(F.dfdivf( z.ll + (z.ul - z.ll)*t ).*(z.ul-z.ll).*wt,2);
I.yu = sum(F.dfdivf( z.lu + (z.uu - z.lu)*t ).*(z.uu-z.lu).*wt,2);

Z.z = round(I.yl + I.xu - I.yu - I.xl,1) ;
disp([num2str(Z.z),' zeros in box']);

% Real zeros

if (0<y2) && (0>y1)
    z.ll = x1 - 1i*prec;
    z.lu = x1 + 1i*prec;
    z.ul = x2 - 1i*prec;
    z.uu = x2 + 1i*prec;
    
    I.xl = sum(F.dfdivf( z.ll + (z.lu - z.ll)*t ).*(z.lu-z.ll).*wt,2);
    I.xu = sum(F.dfdivf( z.ul + (z.uu - z.ul)*t ).*(z.uu-z.ul).*wt,2);
    I.yl = sum(F.dfdivf( z.ll + (z.ul - z.ll)*t ).*(z.ul-z.ll).*wt,2);
    I.yu = sum(F.dfdivf( z.lu + (z.uu - z.lu)*t ).*(z.uu-z.lu).*wt,2);

    Z.r = round(I.yl + I.xu - I.yu - I.xl,1) ;
    disp([num2str(Z.r),' real zeros in box']);
else
    Z.r = [] ;
    disp('Real line not included in the box');
end

if 0<y2
    z.ll = x1 + 1i*prec;
    z.lu = x1 + 1i*y2;
    z.ul = x2 + 1i*prec;
    z.uu = x2 + 1i*y2;
      
    I.xl = sum(F.dfdivf( z.ll + (z.lu - z.ll)*t ).*(z.lu-z.ll).*wt,2);
    I.xu = sum(F.dfdivf( z.ul + (z.uu - z.ul)*t ).*(z.uu-z.ul).*wt,2);
    I.yl = sum(F.dfdivf( z.ll + (z.ul - z.ll)*t ).*(z.ul-z.ll).*wt,2);
    I.yu = sum(F.dfdivf( z.lu + (z.uu - z.lu)*t ).*(z.uu-z.lu).*wt,2);

    Z.U = round(I.yl + I.xu - I.yu - I.xl,1) ;
    disp([num2str(Z.U),' zeros in box with positive imaginary part']);
else
    Z.U = [] ;
    disp('No zeros in box with positive imaginary part - box lies in lower half plane');
end

if 0>y1
    z.ll = x1 + 1i*y1;
    z.lu = x1 - 1i*prec;
    z.ul = x2 + 1i*y1;
    z.uu = x2 - 1i*prec;
    
    I.xl = sum(F.dfdivf( z.ll + (z.lu - z.ll)*t ).*(z.lu-z.ll).*wt,2);
    I.xu = sum(F.dfdivf( z.ul + (z.uu - z.ul)*t ).*(z.uu-z.ul).*wt,2);
    I.yl = sum(F.dfdivf( z.ll + (z.ul - z.ll)*t ).*(z.ul-z.ll).*wt,2);
    I.yu = sum(F.dfdivf( z.lu + (z.uu - z.lu)*t ).*(z.uu-z.lu).*wt,2);

    Z.L = round(I.yl + I.xu - I.yu - I.xl,1) ;
    disp([num2str(Z.L),' zeros in box with negative imaginary part']);
else
    Z.L = [] ;
    disp('No zeros in box with negative imaginary part - box lies in upper half plane');
end

if Z.r + Z.U + Z.L ~= Z.z
    disp('Warning: Zeros dont add up? Possible precision faliure')
end
end