function [ZP]=ZerosPolesBranchCuts(fun,dfun,x,y)

% computes Z-P for function fun in the box x1<Re(z)<x2, y1<Im(z)<y2;

% That is, zeros-poles (with multiplicity)
% Analytic input of fun and its derivative (currently)

% If there is a branch cut, makes every attempt to integrate around it.

Coarsity=0.001;

Zh=ContourSweep(fun,dfun,x(1)+1i*(y(1):Coarsity:y(2)),1,x(2)-x(1),Coarsity);
%Zv=ContourSweep(fun,dfun,x(1):Coarsity:x(2)+1i*y(1),1i,y(2)-y(1),Coarsity);

%Z=[Zh(:); Zv(:)];
Z=Zh(:); %for now ignoring any horizontal branch cuts.

f=@(z) dfun(z)./(2*pi*1i*fun(z));

if numel(Z)==0
    ZP= integral(f,x(1)+1i*y(1),x(1)+1i*y(1),'Waypoints',[x(2)+1i*y(1),x(2)+1i*y(2),x(1)+1i*y(2)]);
else
    B = BranchCuts(Z,100*Coarsity);
    if numel(B,1)>1
        ZP=3; % This is a fudge to ensure that the boxes are halved until
        % only a single branch point is contained within them
    else
        ZP1= integral(f,x(1)+1i*y(1),x(1)+1i*y(1),'Waypoints',[x(2)+1i*y(1),x(2)+1i*y(2),x(1)+1i*y(2)]);
        % probably need to redefine ZP1 as integrating over the branch cut
        % here
        [ZP2,ZP3]=BranchCutIntegral(f,B,Coarsity*0.1);
        ZP=ZP1+ZP2+ZP3;
    end
end
    

ZP=real(ZP); % ZP is pretty much always real unless something's gone wrong.
%Also mostly integral, though branch points etc can make that skewy.

end