function [z] = ZeroFinder_Recursive(f,z0,varargin)
% Finds zeros using a NR method, with initial guesses given by z0

% Agressively makes sure it finds all the zeros, firstly locating them and
% their complex sign using Zeros_Box, encorporated as structure Z, if it is
% chosen to be used (as the final input). Otherwise, searches for zeros as
% in z0.

init.r = z0(imag(z0)==0) ; % real z0
init.U = z0(imag(z0)>0) ;
init.L = z0(imag(z0)<0) ;

if nargin == 3
    Z = varargin{1};
else
    Z.z = numel(z0);
    Z.r = numel(init.r);
    Z.U = numel(init.U);
    Z.L = numel(init.L);
end

if numel(init.r) ~= Z.r
    disp('Warning: guess has different number of real zeros to the function') ;
elseif numel(init.U) ~= Z.U
    disp('Warning: guess has different number of positive imaginary zeros to the function') ;
elseif numel(init.L) ~= Z.L
    disp('Warning: guess has different number of negative imaginary zeros to the function') ;
end

z.z=zeros(1,min(numel(z0),Z.z));
if numel(z0) < Z.z
    disp('Warning: initial guess lower than the total number of zeros, some will be missed')
elseif numel(z0) > Z.z
    disp('More elements in initial guess than zeros of the function: ignoring latter guesses')
end

prec_diff = 1e-10; % differentiation precision
prec_real = 1e-5; % how real a real number is

F.f=@(k) f(k);
F.df = @(k) ( f(k+prec_diff) - f(k-prec_diff) ) ./(2*prec_diff) ;

% Number of each sort of zero to be found
N.r = min(Z.r,numel(z0)) ;
N.U = min(Z.U,numel(z0)-N.r);
N.L = min(Z.L,numel(z0)-N.r-N.U);
N.z=N.r+N.U+N.L ;
if N.z ~= min(numel(z0),Z.z)
    error('Code error 1: something is wrong with the number counting')
end

% Real zeros first of all
z.r = zeros(1,N.r) ; % Real zeros
j=1;
if numel(z.r) == numel(init.r)
    guess.r = sort(init.r);
elseif numel(z.r) > numel(init.r)
    guess.r = sort([init.r,init.r(end)*ones(1,numel(z.r)-numel(init.r))]);
else
    guess.r = sort(init.r(1:numel(z.r)));
end

count = 0 ; % an infinite loop check

while j<N.r+1
    z.r(j) = NewtonRaphson(F.f,F.df,guess.r(j)) ;
    if isnan(z.r(j)) == 1
        guess.r(j) = 1.1*guess.r(j) ; % give it a nudge
        count = count + 1;
        if count > 10
            error('A zero is being too NaNNy')
        end
    elseif abs(imag(z.r(j)))>prec_real
        guess.r(j) = 1.1*guess.r(j) ; % give it a nudge, don't want complex answers
        count = count + 1;
        if count > 10
            error('A real zero is trying to be complex')
        end
    else
        M=Multiplicity(F,z.r(j)); % multiplicity check, needs to be an integer
        if M<1
            error('Something fishy about the multiplicity')
        elseif M ==1
            F.f = @(k) F.f(k)./(k-z.r(j)); % removes the zero
            F.df = @(k) F.df(k)./(k-z.r(j))-F.f(k)./((k-z.r(j)).^2) ;
            disp(['Found ',num2str(j),' real zeros'])
            j=j+1;
            count = 0;
        elseif M>1
            z.r(j+1:1:j+M-1) = z.r(j) ; % assigns multple zeros
            F.f = @(k) F.f(k)./((k-z.r(j)).^M); % removes all the zeros
            F.df = @(k) F.df(k)./((k-z.r(j)).^M)-M*F.f(k)./((k-z.r(j)).^(M+1)) ;
            disp(['Found ',num2str(j+M-1),' real zeros'])
            j=j+M;
            count = 0;
        end
    end
end

if numel(z.r) == 0
    disp('Found 0 real zeros') ;
end

% Zeros in the UHP
% Note: F is still reduced by the removal of real zeros
z.U = zeros(1,N.U) ; % Imaginary positive zeros
j=1; % reset j
if numel(z.U) == numel(init.U)
    guess.U = sort(init.U);
elseif numel(z.U) > numel(init.U)
    guess.U = sort([init.U,init.U(end)*ones(1,numel(z.U)-numel(init.U))]);
else
    guess.U = sort(init.U(1:numel(z.U)));
end

count = 0 ; % an infinite loop check

while j<N.U+1
    z.U(j) = NewtonRaphson(F.f,F.df,guess.U(j)) ;
    if isnan(z.U(j)) == 1
        guess.U(j) = 1.1*guess.U(j) ; % give it a nudge
        count = count + 1;
        if count > 10
            error('A UHP zero is being too NaNNy')
        end
    elseif imag(z.U(j))<prec_real
        guess.U(j) = 1.1*guess.U(j) ; % give it a nudge, don't want complex answers
        count = count + 1;
        if count > 10
            error('A UHP zero is giving LHP answers')
        end
    else
        M=Multiplicity(F,z.U(j)); % multiplicity check, needs to be an integer
        if M<1
            error('Something fishy about the multiplicity')
        elseif M ==1
            F.f = @(k) F.f(k)./(k-z.U(j)); % removes the zero
            F.df = @(k) F.df(k)./(k-z.U(j))-F.f(k)./((k-z.U(j)).^2) ;
            disp(['Found ',num2str(j),' zeros in the UHP'])
            j=j+1;
            count = 0;
        elseif M>1
            z.U(j+1:1:j+M-1) = z.U(j) ; % assigns multple zeros
            F.f = @(k) F.f(k)./((k-z.U(j)).^M); % removes all the zeros
            F.df = @(k) F.df(k)./((k-z.U(j)).^M)-M*F.f(k)./((k-z.U(j)).^(M+1)) ;
            disp(['Found ',num2str(j+M-1),' zeros in the UHP'])
            j=j+M;
            count = 0;
        end
    end
end

if numel(z.U) == 0
    disp('Found 0 zeros in the UHP') ;
end

% Zeros in the LHP
z.L = zeros(1,N.L) ; % Imaginary negative zeros
j=1; % reset j
if numel(z.L) == numel(init.L)
    guess.L = sort(init.L);
elseif numel(z.L) > numel(init.L)
    guess.L = sort([init.L,init.L(end)*ones(1,numel(z.L)-numel(init.L))]);
else
    guess.L = sort(init.L(1:numel(z.L)));
end

count = 0 ; % an infinite loop check

while j<N.L+1
    z.L(j) = NewtonRaphson(F.f,F.df,guess.L(j)) ;
    if isnan(z.L(j)) == 1
        guess.L(j) = 1.1*guess.L(j) ; % give it a nudge
        count = count + 1;
        if count > 10
            error('A LHP zero is being too NaNNy')
        end
    elseif imag(z.L(j))<prec_real
        guess.L(j) = 1.1*guess.L(j) ; % give it a nudge, don't want complex answers
        count = count + 1;
        if count > 10
            error('A LHP zero is giving UHP answers')
        end
    else
        M=Multiplicity(F,z.L(j)); % multiplicity check, needs to be an integer
        if M<1
            error('Something fishy about the multiplicity')
        elseif M ==1
            F.f = @(k) F.f(k)./(k-z.L(j)); % removes the zero
            F.df = @(k) F.df(k)./(k-z.L(j))-F.f(k)./((k-z.L(j)).^2) ;
            disp(['Found ',num2str(j),' zeros in the UHP'])
            j=j+1;
            count = 0;
        elseif M>1
            z.L(j+1:1:j+M-1) = z.L(j) ; % assigns multple zeros
            F.f = @(k) F.f(k)./((k-z.L(j)).^M); % removes all the zeros
            F.df = @(k) F.df(k)./((k-z.L(j)).^M)-M*F.f(k)./((k-z.L(j)).^(M+1)) ;
            disp(['Found ',num2str(j+M-1),' zeros in the UHP'])
            j=j+M;
            count = 0;
        end
    end
end
if numel(z.L) == 0
    disp('Found 0 zeros in the LHP') ;
end
        
z.z = [z.r,z.U,z.L];


end

function M = Multiplicity(F,z0)
    N_int = 100 ;
    prec_mult = 1e-3; % multiplicity check precision
    [t,wt] = lgwt(N_int,0,2*pi) ;
    t = t(end:-1:1).' ; wt = wt(end:-1:1).';
    s = prec_mult.*exp(1i*t);
    M = round(sum(wt.*s.*F.df(z0+s)./(2*pi*F.f(z0+s)),2),1);
end
