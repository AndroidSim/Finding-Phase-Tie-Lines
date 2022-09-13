function [kp,chisq,bases,p_theory] = kp_bd_fit(S,c)

% S = spectra
% c = compositions of the spectra

warning off
options = optimset('display','off');
nc = length(c);
[np,ns] = size(S);
kp = repmat(NaN,[nc nc]);
chisq = repmat(NaN,[nc nc]);

% B = basis set 
% p = coefficients (fraction of probe in each phase)
% S(:,i) = spectrum
        
% choose basis
% B = [S(:,1) S(:,20)];

if ns ~= nc
    error('the number of spectra must = the number of compositions');
end

for i = 2:nc-3
    alpha = c(i);
    
    for k = i+2:nc-1
        beta = c(k);
        
        if alpha == beta
            continue
        end
        
        % choose basis
        B = [S(:,i) S(:,k)];
        
        % linear least squares with constraints on p
        
        % constraints
        Aeq = [1 1];
        beq = 1;
        A = [-1 0;0 -1];
        b = [0;0];
        
        for q = 1:nc
            p(:,q) = lsqlin(B,S(:,q),A,b,Aeq,beq,[],[],[],options);% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options)
        end
        
        p = p';
        % a = find(p > 1);
        % p(a) = 1;
        % a = find(p < 0);
        % p(a) = 0;
        bases{i,k} = p;
        
        [x,resnorm] = lsqcurvefit(@kp_bd_fit_fun,1,c,p(:,1),0.1,10,options,alpha,beta); % c(i:k),p(i:k,1)
        
        p_theory{i,k} = feval(@myfun,x,c,alpha,beta);
        kp(i,k) = x; 
        chisq(i,k) = resnorm;
        clear p;
    end
end

warning on
return


% linear no constraints
%for i = 1:20
%x(:,i) = B\S(:,i);
%end

% linear least squares with constraints on x
%for i = 1:20
%y(:,i) = lsqlin(B,S(:,i),[],[],[1 1],1);
%end

%x = lsqcurvefit(@myfun,[5 0.5 0.5],[1:17]',y_ratio(4:end),[1 0 0],[10 1 1]);