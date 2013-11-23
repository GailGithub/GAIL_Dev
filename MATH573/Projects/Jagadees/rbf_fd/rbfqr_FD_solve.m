function coeffs = rbfqr_FD_solve(x,y,ep,alpha,M)
% function rbfqrOBJ = rbfqr_FD_solve(x,y,ep,alpha,M)
% This function accepts required inputs x, y and ep
% x is a Nxd vector of input data points of the form
%     x = [x1';x2';...;xN'], x1 is a column vector
% y is a Nx1 vector of data values at the x locations
% ep is the traditional RBF shape parameter
%
% There is also a recommended input alpha which is the
% global scale parameter which determines the
% orthogonality enjoyed by the eigenfunctions
% The default is set in rbfsetup, which can be used by passing []
% For the time being alpha is only a single value
% In the future we may allow for different alpha in each dimension
%
% Note that only positive values of ep and alpha are used
%
% Optional input is M>N which is the length of the
% expansion.  If you don't pass M, this will choose an
% M based on the eigenvalues
%
% What is returned from this function is rbfqrOBJ which is
% a large object that encapsulates everything you need to
% do RBF-QR.  It is packaged like this to make it easier for
% you to call other rbfqr functions.

	% Import global parameters from rbfsetup
	global GAUSSQR_PARAMETERS
	if ~isstruct(GAUSSQR_PARAMETERS)
		error('GAUSSQR_PARAMETERS does not exist ... did you forget to call rbfsetup?')
	end
	alertuser = GAUSSQR_PARAMETERS.WARNINGS_ON;

	if nargin<3
		error('Insufficient inputs')
	end
	rbfqrOBJ.warnid = '';
	rbfqrOBJ.warnmsg = '';

	[N,d] = size(x);
	if sum(N~=size(y,1))
		error('Different number of input and output points (x and y)')
	elseif size(y,2)~=1
		error('You can only pass a 1D output vector y')
	end

	if nargin==3
		[ep,alpha,Marr,lam] = rbfsolveprep(0,x,ep);
	elseif nargin==4
		[ep,alpha,Marr,lam] = rbfsolveprep(0,x,ep,alpha);
	else
		[ep,alpha,Marr,lam] = rbfsolveprep(0,x,ep,alpha,M);
	end
	phiMat = rbfphi(Marr,x,ep,alpha);
	if d==2
		phiMatL = rbfphi(Marr,x,ep,alpha,[2,0]) + rbfphi(Marr,x,ep,alpha,[0,2]); %Laplacian in 2D
	elseif d==3
		phiMatL = rbfphi(Marr,x,ep,alpha,[2,0,0]) + rbfphi(Marr,x,ep,alpha,[0,2,0]) + rbfphi(Marr,x,ep,alpha,[0,0,2]); %Laplacian in 3D
	end

	[Q,R] = qr(phiMat);
	R1 = R(:,1:N);
	R2 = R(:,N+1:end);

	lastwarn('')
	warning off MATLAB:divideByZero
	iRdiag = diag(1./diag(R1));
	[warnmsg,msgid] = lastwarn;
	if strcmp(msgid,'MATLAB:divideByZero')
		rbfqrOBJ.warnid = 'GAUSSQR:zeroQRDiagonal';
		rbfqrOBJ.warnmsg = 'At least one value on the R diagonal was exactly 0';
	end
	warning on MATLAB:divideByZero

	R1s = iRdiag*R1;
	opts.UT = true;

	lastwarn('')
	warning off MATLAB:singularMatrix
	if strcmp(rbfqrOBJ.warnid,'GAUSSQR:zeroQRDiagonal')
		Rhat = linsolve(R1s,iRdiag*R2,opts);
	else
		Rhat = linsolve(R1s,iRdiag*R2,opts);
		[warnmsg,msgid] = lastwarn;
		if strcmp(msgid,'MATLAB:singularMatrix')
			rbfqrOBJ.warnid = 'GAUSSQR:singularR1invR2';
			rbfqrOBJ.warnmsg = 'Computing inv(R1)R2 ... R1 singular to working precision';
		end
	end

	% Here we apply the eigenvalue matrices
	% Note that the -d term in the power goes away because it
	% appears in both Lambda_2 and Lambda_1^{-1}
	Ml = size(Marr,2);
	D = lam.^(repmat(sum(Marr(:,N+1:end),1)',1,N)-repmat(sum(Marr(:,1:N),1),Ml-N,1));
	Rbar = D.*Rhat';

	psi = phiMat*[eye(N);Rbar];
	psiL = phiMatL*[eye(N);Rbar];
	D1 = lam.^(sum(Marr(:,1:N),1));
	D1 = diag(D1);
	h = (pinv(D1*(R1')*(Q'))')*ones(N,1);
	A = [psi' h; ones(1,N) 0];
	B = [psiL(1,:)';0];

	warning('off', 'MATLAB:singularMatrix');
	warning off all;
	coeffs = linsolve(A,B);
	warning on all;
	if 0
		[coef,recipcond] = ranksolve(Rhat,Rbar,linsolve(R1s,iRdiag*(Q'*y),opts));
		if (recipcond<eps || isnan(recipcond)) && strcmp(rbfqrOBJ.warnid,'')
			rbfqrOBJ.warnid = 'GAUSSQR:illConditionedRanksolve';
			rbfqrOBJ.warnmsg = sprintf('ranksolve encountered an ill-conditioned system, rcond=%g',recipcond);
		end

		warning('on', 'MATLAB:singularMatrix');

		rbfqrOBJ.reg   = false;
		rbfqrOBJ.ep    = ep;
		rbfqrOBJ.alpha = alpha;
		rbfqrOBJ.N     = N;
		rbfqrOBJ.coef  = coef;
		rbfqrOBJ.Rbar  = Rbar;
		rbfqrOBJ.Marr  = Marr;

		if alertuser && ~strcmp(rbfqrOBJ.warnid,'')
			warning(rbfqrOBJ.warnid,rbfqrOBJ.warnmsg)
		end
	end

% Developer's note: I should throw something in about handling the
% ill-conditioning.  Users should be alerted when things are bad.
% Specifically it is possible to get some of the diagonal values of R
% exactly equal to zero, which is obviously a problem to which the user
% should be alerted.
% Relevant warning codes include MATLAB:illConditionedMatrix
%                                MATLAB:divideByZero
%                                MATLAB:nearlySingularMatrix
%                                MATLAB:singularMatrix
% I'll need to use warning query to check along the way
%
% Eventually I'll need to find a way to allow users to simultaneously pick
% M and alpha.  Right now I have them guess at M in the rbfsetup file and
% then use that guess to come up with alpha, which then comes up with the
% real M.  This works because there is an acceptable alpha range, but
% there's a better way to do this, I just don't know it yet.