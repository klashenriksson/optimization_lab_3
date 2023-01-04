    % Get points on the circle.
pts=circdata(1,11,0);

% Weight matrix
W=eye(numel(pts));

% Common constants
epsR=1e-6;
maxIter=50;
mu=0.1;
alphaMin=1e-3;
epsC=1e-8;
nu0=0.1;

for explicit=[true,false]
    if explicit
        fig=1;

        % Create a naive estimate of the initial values.
        [c,r,theta]=circle_x0naive(pts);
        % Initial value for the explicit formulation
        x0e=[c(:);r;theta(:)];

        [x,code,n,X,alphas]=gaussn_niclas_damped(@circle_r_explicit,x0e,epsR, ...
                                                 maxIter,alphaMin,{pts});
        code,n

        % Re-compute the residual and Jacobian at the solution.
        [r,J]=circle_r_explicit(X(:,end),pts);

        % Extract c, r for explicit formulation.
        ce=X(1:2,end);
        re=X(3,end);
        % Compute the standard deviation of unit weight.
        sigma0e=sqrt(norm(r)^2/(size(J,1)-size(J,2)));
        % Compute the posterior covariance matrix.
        Cxxe=sigma0e^2*inv(J'*W*J);
        % Compute the posterior standard deviations for the
        % global parameters
        sxxe=sqrt(diag(Cxxe(1:3,1:3)))';
        % Compute the correlation matrix for the circle parameters
        rxxe=corrmat(Cxxe(1:3,1:3),true);

        % Compute redundancy numbers
        Pe=inv(J'*W*J)*J'*W;
        Ue=J*Pe;
        Re=eye(size(Ue))-Ue;
        riie=diag(Re);
        rPtsE=sum(reshape(riie,2,[]));
        rPts=rPtsE;
    else
        fig=2;

        % Create a naive estimate of the initial values.
        [c,r,theta]=circle_x0naive(pts);
        % Initial value for the implicit formulation
        x0i=[c;r;pts(:)];
        
        [x,n,code,l,X,alphas,C,L,nus]=sqpsq(@circle_r_implicit,...
                                            @circle_c,x0i,epsR, ...
                                            epsC,maxIter,{pts},nu0,mu);
        
        code,n

        % Extract c, r for implicit formulation.
        ci=X(1:2,end);
        ri=X(3,end);

        % Re-compute residual, constraint + jacobians at solution.
        [r,J]=circle_r_implicit(x,pts);
        [c,A]=circle_c(x,pts);

        % Compute the standard deviation of unit weight.
        sigma0i=sqrt(norm(r)^2/(size(J,1)-size(J,2)+length(c)));

        % Form the system matrix
        S=[J'*W*J,A';A,zeros(size(A,1))];
        % Compute the posterior covariance matrix.
        % We need part of the inverse for Cxx.
        N=inv(S);
        % Matrix to extract first n elements
        E=eye(size(J,2),size(J,2)+length(c));
        % Projection matrix.
        Pi=E*N*E'*J'*W;
        % Posterior covariance matrix
        Cxxi=sigma0i^2*Pi*Pi';
        % Compute the posterior standard deviations for the
        % global parameters
        sxxi=sqrt(diag(Cxxi(1:3,1:3)))';
        % Compute the correlation matrix for the circle parameters
        rxxi=corrmat(Cxxi(1:3,1:3),true);

        % Compute redundancy numbers
        Ui=J*Pi;
        Ri=eye(size(Ui))-Ui;
        riii=diag(Ri);
        rPtsI=sum(reshape(riii,2,[]));
        rPts=rPtsI;
    end
    
    % Plot measured points with redundancy numbers
    figure(fig)
    h=subplot(2,2,3,'parent',fig);
    plot(h,pts(1,:),pts(2,:),'k*');
    axis(h,'equal')
    
    for i=1:size(pts,2)
        text(h,pts(1,i),pts(2,i),sprintf('%.2f',rPts(i)),'vertical', ...
             'top','horizontal', 'center');
    end

    title(h,'Redundancy numbers')

    % Compute the objective function value, constraint value, and the
    % estimated circles for each iteration. Each layer estPts(:,:,i)
    % contain the estimated points for iteration i.
    estPts=nan(2,size(pts,2),size(X,2));
    resNorm=nan(1,size(X,2));
    conVal=nan(1,size(X,2));
    for i=1:size(X,2)
        x=X(:,i);
        c=x(1:2);
        r=x(3);
        if explicit
            theta=x(4:end);
            resNorm(i)=norm(circle_r_explicit(x,pts));
            estPts(:,:,i)=circle_g(c,r,theta);
        else
            resNorm(i)=norm(circle_r_implicit(x,pts));
            estPts(:,:,i)=reshape(x(4:end),2,[]);
            conVal(i)=max(abs(circle_c(x,pts)));
        end
    end
    
    % Plot residual evolution
    h=subplot(3,2,2,'parent',fig);
    semilogy(h,0:size(X,2)-1,resNorm,'x-')
    title(h,'norm(r)')
    xlim(h,[0,size(X,2)-1])
    
    % Plot constraint values
    h=subplot(3,2,4,'parent',fig);
    semilogy(h,0:size(X,2)-1,conVal,'x-')
    title(h,'max(abs(c))')
    xlim(h,[0,size(X,2)-1])
    
    % Plot step lengths
    h=subplot(3,2,6,'parent',fig);
    plot(h,1:length(alphas),alphas,'x-')
    ylim(h,[0,1.1])
    xlim(h,[0,length(alphas)])
    title(h,'alphas')
    
    
    % Plot observed points again for iterations
    h=subplot(2,2,1,'parent',fig);
    plot(h,pts(1,:),pts(2,:),'k*');
    axis(h,'equal')

    % Plot successive circles
    for i=1:size(X,2)
        if i==1
            % Dark green for initial estimate
            color=[0,0.5,0];
        elseif i==size(X,2)
            % Red color for final estimates
            color='r';
        else
            % Blue for intermediate values.
            color='b';
        end
        x=X(:,i);
        c=x(1:2);
        r=x(3);
        t=linspace(0,2*pi,101);
        % Points on the full circle
        circlePtsFull=circle_g(c,r,t);
    
        % Plot estimated points
        line(h,estPts(1,:,i),estPts(2,:,i),'marker','x', ...
             'linestyle','none', 'color',color)
        % Plot full circle
        line(h,circlePtsFull(1,:),circlePtsFull(2,:),'marker', ...
             'none','color',color)
        
        if explicit
            title(h,'Circle estimates, explicit version');
        else
            title(h,'Circle estimates, implicit version');
        end
    end

end

% Define some useful functions...
absErr=@(A,B)abs(A-B);
relErr=@(A,B)absErr(A,B)./abs(A);
