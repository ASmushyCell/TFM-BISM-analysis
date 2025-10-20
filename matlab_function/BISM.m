%% Bayesian Inversion Stress Microscopy
%% MATLAB script
%% Vincent NIER, 15/11/16

%If you use BISM please cite the following reference in your work: 
%Nier, V. et al., Inference of internal stress in a cell monolayer, Biophysical Journal, 110(7), 1625-1635, 2016

% Input: Traction_field.mat, traction force data
% Output: Stress_BISM.mat, internal stress tensor
% Parameters: C, number of data points along the x direction, between xmin and xmax
%             R, number of data points along the y direction, between ymin and ymax		      
%		     kM, number of time steps				
% The regularization parameter Lambda can be estimated by either (meth_Lambda)
%             Maximal A Posteriori (MAP) optimization
%             L-curve method 
% Free stress boundary conditions can be imposed in the prior (BC).
%
%%%%%%% MODIFs by JdA
% To reduce the loading time, now 'ForceName' is a Matlab file containing
% kM variables, with each eval(['frame_' num999(k0)]) variable
% corresponding to the k0th frame.
%
% Also the building of the discretised divergence operator has been changed
% to prevent memory overload.

% clear
close all;

% coeff = 0.64; % conversion coefficient from pixels to micrometers

%% Traction force data
% Path = path;
% ForceName = [Path 'Stress\TFM_to_BISM_' nom_image '_dx=' num2str(dx) '_overlap=' num2str(overlap*100) '%.mat'];
ForceName = [Path 'tractionforce.mat'];  % changes made by LB for ease of reading file
load(ForceName,'frame_001');
R=size(frame_001.Tx,1); % number of rows (y direction)
C=size(frame_001.Tx,2); % number of columns (x direction)
% R=size(frame_001.Tx,1); % number of rows (y direction)
% C=size(frame_001.Tx,2); % number of columns (x direction)

k1=KMIN; % starting frame which to start, usually set to 1 but in some cases set to a later frame
kM=KMAX; % number of time steps
N=R*C; % system size for traction force field
Ninf=4*N+2*(C+R); % system size for stress tensor field
 
%% Define spatial grid (rectangular)

% load ROI parameters
% load([Path 'Stress\TFM_to_BISM_' nom_image '_dx=' num2str(dx) '_overlap=' num2str(overlap*100) '%.mat'],'xmin_pix','xmax_pix','ymin_pix','ymax_pix');
load(ForceName,'xmin_pix','ymin_pix','xmax_pix','ymax_pix');
xmin = coeff*xmin_pix;
xmax = coeff*xmax_pix;
ymin = coeff*ymin_pix;
ymax = coeff*ymax_pix;
% uniform and isotropic spatial resolution (lx=ly=l)
l=(xmax-xmin)/(C-1); 

x=xmin:l:xmax; % x coordinate (column from left to right)
y=ymin:l:ymax; % y coordinate (row from top to bottom)
 
%% Parameters
BC=JEDI; 
% BC = 1 for a confined tissue
% BC = 0 otherwise
meth_Lambda=CYRILLE; 
% method to determine the regularization parameter Lambda
% 1: MAP 
% 2: L-curve 
% 3: fixed value
noise_value=0.0; 
% noise amplitude of the traction force data 
% (same unit as traction force data)
% set to zero when unknown
fplot=0; 
% Graphics
% 0: no figures
% 1: plot figures
mult=10^(0); 
% multiplicative coefficient for a better graphical representation of the
% stress tensor

%% Matrices A and B
 
%% Building 2N*Ninf matrix A such that A*\sigma=T (discretized divergence operator)
% Building matrix Ax (derivative wrt x)
Ax1=diag(ones(C-1, 1), 1)-diag(ones(C, 1), 0); %bulk (order 2)
Ax1=[Ax1,[zeros(C-1, 1); 1]];
Ax=sparse(Ax1);
for i=1:R-1
    Ax=blkdiag(Ax,Ax1); %iteration R times for all the rows
end
clear Ax1;

%%%%% Memory issue => trying to build too large Ay1 matrix (not sparse)
% % Building matrix Ay (derivative wrt y)
% Ay1=diag(ones((R-1)*C,1), C)-diag(ones(R*C,1 ), 0); %bulk (order 2)
% Ay=sparse([Ay1, [zeros((R-1)*C,C); eye(C)]]);
% clear Ay1;
%%%%%
%%%%%%% Replaced by 
Ay1 = speye(R*C);
Ay = [-Ay1 zeros(R*C,C)] + [zeros(R*C,C) Ay1];
%%%%%%%%

% Building A from Ax and Ay
A=[Ax, sparse(N,N+C) Ay, sparse(N,N+R); sparse(N,N+R) Ay, sparse(N,N+C), Ax]/l;
clear Ax Ay; 

%% Building Ninf*Ninf B covariance matrix of the prior, S_0^{-1}=B/s_0^2
%% Stress norm regularization
B=speye(Ninf);
 
%% Add to the prior a term enforcing the equality of shear components
%%%%% Same issue as before, solved in a similar way. JdA, 2017/02/27
%d1=diag(ones((R-1)*C, 1), C)+diag(ones(R*C,1), 0); %bulk (order 2)
%%%%%
%%%%%% Replaced by
d0 = speye(R*C);
d1 = [d0 zeros(R*C,C)] + [zeros(R*C,C) d0];
d1 = d1(:,1:R*C);
%%%%%%

bd1=[d1,[zeros((R-1)*C, C); eye(C)]];
clear d1;
d2=-diag(ones(C-1, 1), 1)-diag(ones(C, 1),0); %bulk (order 2)
d2=sparse([d2,[zeros(C-1, 1); -1]]);
bd2=d2;
% loop on rows
for i=1:R-1
   bd2=blkdiag(bd2, d2); 
end
clear d2;
Bdiff=sparse([sparse(N, 2*N+C+R), bd1, bd2]);
clear bd1 bd2;
alpha_xy=10^3; 
% hyperparameter for the equality of shear components (has to be large)
B=B+alpha_xy^2*sparse(Bdiff'*Bdiff);
clear Bdiff;
 
%% Add to the prior a term enforcing the free stress boundary conditions
if BC==1
    Bbcx=zeros(2*C, C*(R+1));
    Bbcy=zeros(2*R, (C+1)*R);   
    Bbcx(1:C,1:C)=speye(C);
    Bbcx(C+1:2*C, 1+R*C:(R+1)*C)=speye(C);
    for i=1:R
         Bbcy(i, 1+(i-1)*(C+1))=1;
         Bbcy(i+R, C+1+(i-1)*(C+1))=1;
    end
    Bbc=sparse([Bbcy, sparse(2*R,(C+1)*R+2*C*(R+1)); sparse(2*C,(C+1)*R), Bbcx, sparse(2*C,(C+1)*R+C*(R+1)); sparse(2*C,(C+1)*R+C*(R+1)), Bbcx, sparse(2*C,(C+1)*R); sparse(2*R,(C+1)*R+2*C*(R+1)), Bbcy]); %xx, yy, yx and xy components
    clear Bbcx Bbcy; 
    % hyperparameter for the free stress boundary condition (has to be large)
    alpha_BC=10^3; 
    B=B+alpha_BC^2*sparse(Bbc'*Bbc);
    clear Bbc;
end
 
%% Traction force
% Recording the L-curve dynamics
DDL = cell(kM,1);
SIG_NO = cell(kM,1);
% loop on the time step
for k0=k1:kM 

    %load(sprintf('%s',ForceName)); % load the traction force data
    %f=sprintf('frame%d', k0); % load the k0 time step
    %T=traction.(f);
    
    % Joseph's data handling: load one single frame instead of loading
    % everything and clearing
    f=['frame_' num2str(k0,'%03d')]; % load the k0 time step
    load(ForceName,f);
    
    eval(['mTx= ' f '.Tx;']) % Tx in a R*C matrix form
    eval(['mTy= ' f '.Ty;']) % Ty in a R*C matrix form 
    clear (f); 
    clear Traction_field;
 
    vTx=reshape(mTx', N, 1); % Tx in a N vector form
    vTy=reshape(mTy', N, 1); % Ty in a N vector form 
    T=[vTx; vTy]; % vector T (2N components)
    [xg, yg]=meshgrid(x,y);
    
%% Figure of the traction force field
    if fplot==1
        figure(1000+k0) %field vector T
        quiver(xg,yg,mTx,mTy,'b','LineWidth',2);
        set(gca, 'FontSize', 18, 'fontName','Times');
        set(gcf,'Color','w')
        xlabel('x (\mum)', 'Fontsize', 18)
        ylabel('y (\mum)', 'Fontsize', 18)
        title('t', 'Fontsize', 18)
        axis([xmin-0.01*(xmax-xmin) xmax+0.01*(xmax-xmin) ymin-0.01*(ymax-ymin) ymax+0.01*(ymax-ymin)])
    end
    clear xg yg;

%% Hyperparameters

% Lambda
    if meth_Lambda==1
    % MAP estimation with Jeffreys' hyperprior
    % iterative methods 
        % number of steps in the iteration
        stepM=10; 
 
        % initial condition
        step=1;
        Lambda_MAP(step)=10^(-3);
 
        % iteration loop
        while step<stepM
            sigmaMAP=(Lambda_MAP(step)*B+l^2*A'*A)\(l^2*A'*T);
            sigmaMAP_norm(step)=sigmaMAP'*B*sigmaMAP;
            res_norm(step)=(T-A*sigmaMAP)'*(T-A*sigmaMAP);
            step=step+1;
            s0_2=sigmaMAP_norm(step-1)/(4*N+2*(C+R)+2);
            s_2=res_norm(step-1)/(2*N+2);
            Lambda_MAP(step)=l^2*s_2/s0_2;
        end
 
        % graphics
        figure(1)
        semilogy(0:stepM-1,Lambda_MAP,'LineWidth',2)
        set(gca, 'FontSize', 18, 'fontName','Times');
        set(gcf,'Color','w')
        xlabel('Step','Fontsize',18)
        ylabel('\Lambda', 'Fontsize',18)
        title('MAP estimate', 'Fontsize',18)
 
        % hyperparameter value   
        Lambda=Lambda_MAP(stepM);
        
        % MAP estimate of the noise amplitude
        noise_value_MAP=sqrt(s_2); 

    elseif meth_Lambda==2
    % (parametric) L-curve
        k=1;
        for lc=-10:0.1:0
             Lambda_Lcurve(k)=10^lc;
             Lsigma_post=(Lambda_Lcurve(k)*B+l^2*A'*A)\(l^2*A'*T);
             sigma_norm(k)=Lsigma_post'*B*Lsigma_post;
             res_norm(k)=(T-A*Lsigma_post)'*(T-A*Lsigma_post);
             k=k+1;
        end
 
        % local slope of the L-curve in logarithmic scale
        dL=(log(sigma_norm(3:end))-log(sigma_norm(1:end-2)))./(log(res_norm(3:end))-log(res_norm(1:end-2)));
        % local curvature of the L-curve
        ddL=2*(dL(3:end)-dL(1:end-2))./(log(res_norm(4:end-1))-log(res_norm(2:end-3)));

        % graphics
        % plot of the L-curve
        figure(1) 
        loglog(res_norm.^0.5,sigma_norm.^0.5,'LineWidth',2)
        set(gca, 'FontSize', 18, 'fontName','Times');
        set(gcf,'Color','w')
        xlabel('$\mathrm{Residual\,norm\,\Vert\vec{T}-\textbf{A}\vec{\sigma}\Vert\,(kPa)}$','interpreter','latex','Fontsize',18)
        ylabel('Prior norm(kPa.\mum)', 'Fontsize',18)
        title('L curve', 'Fontsize',18)
 
        % plot of the local curvature as a function of Lambda
        figure(2)
        semilogx(Lambda_Lcurve(3:end-2),ddL,'LineWidth',2)
        set(gca, 'FontSize', 18, 'fontName','Times');
        set(gcf,'Color','w')
        xlabel('$\Lambda$','interpreter','latex','Fontsize',18)
        ylabel('$\mathrm{\frac{d^2\,log(\vec{\sigma}^T\textbf{B}\vec{\sigma})^{1/2}}{d\,log^2\Vert\vec{T}-\textbf{A}\vec{\sigma}\Vert} (\mu m.kPa^{-1})}$', 'interpreter','latex','Fontsize',18) 
        title('Curvature of the L curve', 'Fontsize',18)
        
        % Recording
        SIG_NO{k0} = sigma_norm;
        DDL{k0} = ddL;
     
        % hyperparameter
        % Lambda at the point of maximal curvature
        [max_curv,index_max]=max(ddL);
        Lambda=Lambda_Lcurve(index_max+2); 

        % L-curve estimate of the noise amplitude
        noise_value_Lcurve=sqrt(res_norm(index_max+2)/(2*N+2)); 
        
    elseif meth_Lambda==3
    % set the value of Lambda
        Lambda = LAMBDA;
        % Default Lambda=10^(-6);
    end

%% Stress and covariance matrix estimations
    sigma_post_obs{k0}=(Lambda*B+l^2*A'*A)\(l^2*A'*T);

    if noise_value ~= 0
        Spost{k0}=noise_value^2*l^2*((Lambda*B+l^2*A'*A)\speye(Ninf));
    end
    
end
clear B;
 
%% Data interpolation and stress representation
for k0=1:kM

    if noise_value ~= 0
        % Inferred variances of xx, yy, xy and yx stress components
        Spost_xx=reshape(diag(Spost{k0}(1:N+R, 1:N+R)), C+1, R)'; 
        Spost_yy=reshape(diag(Spost{k0}(N+R+1:2*N+C+R, N+R+1:2*N+C+R)), C, R+1)'; 
        Spost_xy=reshape(diag(Spost{k0}(2*N+C+R+1:3*N+2*C+R, 2*N+C+R+1:3*N+2*C+R)), C, R+1)'; 
        Spost_yx=reshape(diag(Spost{k0}(3*N+2*C+R+1:4*N+2*(C+R), 3*N+2*C+R+1:4*N+2*(C+R))), C+1, R)'; 
    end        
    vsigma_post_xx=sigma_post_obs{k0}(1:N+R); % inferred stress xx in N+R-vector form
    vsigma_post_yy=sigma_post_obs{k0}(N+R+1:2*N+C+R); % inferred stress yy in N+C-vector form
    vsigma_post_xy=sigma_post_obs{k0}(2*N+C+R+1:3*N+2*C+R); % inferred stress xy in N+C-vector form
    vsigma_post_yx=sigma_post_obs{k0}(3*N+2*C+R+1:4*N+2*(C+R)); % inferred stress yx in N+R-vector form
 
% Interpolation on the grid of the force data
    for i=1:R
        for j=1:C       
            vsigma_xx((i-1)*C+j, 1)=(vsigma_post_xx((i-1)*(C+1)+j)+vsigma_post_xx((i-1)*(C+1)+j+1))/2;
            vsigma_yx((i-1)*C+j, 1)=(vsigma_post_yx((i-1)*(C+1)+j)+vsigma_post_yx((i-1)*(C+1)+j+1))/2;
            vsigma_yy((i-1)*C+j, 1)=(vsigma_post_yy((i-1)*C+j)+vsigma_post_yy((i-1)*C+j+C))/2;
            vsigma_xy((i-1)*C+j, 1)=(vsigma_post_xy((i-1)*C+j)+vsigma_post_xy((i-1)*C+j+C))/2; 
            if noise_value ~= 0
                Ssigma_xx(i,j)=(Spost_xx(i, j)+Spost_xx(i, j+1))/4;
                Ssigma_yx(i,j)=(Spost_yx(i, j)+Spost_yx(i, j+1))/4;
                Ssigma_yy(i,j)=(Spost_yy(i, j)+Spost_yy(i+1, j))/4;
                Ssigma_xy(i,j)=(Spost_xy(i, j)+Spost_xy(i+1, j))/4;
            end
        end
    end
    vsigma_post_xx=vsigma_xx; % inferred stress xx in N-vector form
    vsigma_post_yy=vsigma_yy; % inferred stress yy in N-vector form
    vsigma_post_xy=(vsigma_yx+vsigma_xy)/2; % inferred stress xy in N-vector form
    % Reshape the stress vector in a RxC-matrix form
    sigma_post_xx=reshape(vsigma_post_xx, C, R)'; 
    sigma_post_yy=reshape(vsigma_post_yy, C, R)'; 
    sigma_post_xy=reshape(vsigma_post_xy, C, R)'; 

    if noise_value ~= 0 % estimates of error bars
        Spost_xx=Ssigma_xx;
        Spost_yy=Ssigma_yy;
        Spost_xy=(Ssigma_xy+Ssigma_yx)/4;
        clear Ssigma_xx Ssigma_yy Ssigma_xy Ssigma_yx Spost;
        spostxx=Spost_xx.^0.5; 
        spostyy=Spost_yy.^0.5;
        spostxy=Spost_xy.^0.5;
    end
 
 
%% Figure of the inferred stress tensor field
    if fplot==1
        for i=1:R
            for j=1:C
                % calculation of stress eigenvalues and eigenvectors 
                [V,D]=eig([sigma_post_xx(i,j) sigma_post_xy(i,j); sigma_post_xy(i,j) sigma_post_yy(i,j)]); 
                ab1=[x(j)-mult*D(1,1)*V(1,1),x(j)+mult*D(1,1)*V(1,1)];
                or1=[y(i)-mult*D(1,1)*V(2,1),y(i)+mult*D(1,1)*V(2,1)];
                if D(1,1)>=0
                    % tensile stress in blue
                    figure(2000+k0)
                    hold on
                    plot(ab1,or1,'b','LineWidth',2);
                else
                    % compressive stress in red
                    figure(2000+k0)
                    hold on
                    plot(ab1,or1,'r','LineWidth',2);
                end
                ab2=[x(j)-mult*D(2,2)*V(2,1),x(j)+mult*D(2,2)*V(2,1)];
                or2=[y(i)+mult*D(2,2)*V(1,1),y(i)-mult*D(2,2)*V(1,1)];
                if D(2,2)>=0
                    % tensile stress in blue
                    figure(2000+k0)
                    hold on
                    plot(ab2,or2,'b','LineWidth',2);
                else
                    % compressive stress in red
                    figure(2000+k0)
                    hold on
                    plot(ab2,or2,'r','LineWidth',2);
                end 
            end
        end
        set(gca, 'FontSize', 18, 'fontName','Times');
        set(gcf,'Color','w')
        xlabel('x (\mum)', 'Fontsize', 18)
        ylabel('y (\mum)', 'Fontsize', 18)
        title('\sigma^{inf}', 'Fontsize', 18)
        axis([xmin-0.01*(xmax-xmin) xmax+0.01*(xmax-xmin) ymin-0.01*(ymax-ymin) ymax+0.01*(ymax-ymin)])
    end
 
%% Accuracy of the inference
    Tinf=A*sigma_post_obs{k0}; 
    % divergence of the inferred stress Tinf =A*sigma_inf
    clear sigma_post_obs{k0};
    Tinfx=Tinf(1:N); 
    Tinfy=Tinf(N+1:2*N);
    clear Tinf;
    R2x=1-((vTx'-Tinfx')*(vTx-Tinfx))/((vTx'-mean(vTx))*(vTx-mean(vTx)));
    R2y=1-((vTy'-Tinfy')*(vTy-Tinfy))/((vTy'-mean(vTy))*(vTy-mean(vTy)));
    % R2: coefficient of determination between the data and the traction forces 
    % obtained by calculating the divergence of the inferred stress 
    % (should be close to 1)
    R2=(R2x+R2y)/2; 

% Inferred stress mean values
    sxxm=sum(vsigma_post_xx)/N;
    syym=sum(vsigma_post_yy)/N;
    sxym=sum(vsigma_post_xy)/N;
 
% Inferred hyperparameter
    stress.Lambda{k0}=Lambda;
    
% Inferred stress values
    stress.sxx{k0}=sigma_post_xx;
    stress.syy{k0}=sigma_post_yy;
    stress.sxy{k0}=sigma_post_xy;
    stress.x{k0}=x;
    stress.y{k0}=y;
    stress.R2_T{k0}=R2;
    stress.mean_from_sigma{k0}=[sxxm,syym,sxym];   

% Stress mean values from traction force
    if BC==1
        sxx_mean=-sum(sum(reshape(vTx, C, R)'.*(ones(R, 1)*(x-(xmax-xmin)/2))))/N;
        syy_mean=-sum(sum(reshape(vTy, C, R)'.*((y-(ymax-ymin)/2)'*ones(1, C))))/N;
        sxy_mean=-sum(sum(reshape(vTy, C, R)'.*(ones(R, 1)*(x-(xmax-xmin)/2))+reshape(vTx, C, R)'.*((y-(ymax-ymin)/2)'*ones(1,C))))/(2*N);
        stress.mean_from_t{k0}=[sxx_mean,syy_mean,sxy_mean];
    end
    
% Values of error bars on the inferred stress
    if noise_value~=0
        stress.error_sxx{k0}=spostxx;
        stress.error_syy{k0}=spostyy;
        stress.error_sxy{k0}=spostxy;
    end
end

% save output
save([Path 'Stress_BISM.mat'],'stress');


%% Additional figures 
if fplot==1
    clear stress;
    clear V D ab1 or1 ab2 or2;
    StressName='Stress_field.mat';

    for k0=1:kM 
        
        % Figure of the true stress tensor field
        load(sprintf('%s',StressName)); % load the computed stress 
        f=sprintf('frame%d', k0); % load the k0 time step
        s=stress.(f);
        true_stress_xx=s.sxx; % sxx in a R*C matrix form
        true_stress_yy=s.syy; % syy in a R*C matrix form 
        true_stress_xy=s.sxy; % sxy in a R*C matrix form
        clear s; 
        clear Stress_field;
    
        for i=1:R
            for j=1:C
                % calculation of stress eigenvalues and eigenvectors 
                [V,D]=eig([true_stress_xx(i,j) true_stress_xy(i,j); true_stress_xy(i,j) true_stress_yy(i,j)]); 
                ab1=[x(j)-mult*D(1,1)*V(1,1),x(j)+mult*D(1,1)*V(1,1)];
                or1=[y(i)-mult*D(1,1)*V(2,1),y(i)+mult*D(1,1)*V(2,1)];
                if D(1,1)>=0
                    % tensile stress in blue
                    figure(3000+k0)
                    hold on
                    plot(ab1,or1,'b','LineWidth',2);
                else
                    % compressive stress in red
                    figure(3000+k0)
                    hold on
                    plot(ab1,or1,'r','LineWidth',2);
                end
                ab2=[x(j)-mult*D(2,2)*V(2,1),x(j)+mult*D(2,2)*V(2,1)];
                or2=[y(i)+mult*D(2,2)*V(1,1),y(i)-mult*D(2,2)*V(1,1)];
                if D(2,2)>=0
                    % tensile stress in blue
                    figure(3000+k0)
                    hold on
                    plot(ab2,or2,'b','LineWidth',2);
                else
                    % compressive stress in red
                    figure(3000+k0)
                    hold on
                    plot(ab2,or2,'r','LineWidth',2);
                end 
            end
        end
        set(gca, 'FontSize', 18, 'fontName','Times');
        set(gcf,'Color','w')
        xlabel('x (\mum)', 'Fontsize', 18)
        ylabel('y (\mum)', 'Fontsize', 18)
        title('\sigma^{num}', 'Fontsize', 18)
        axis([xmin-0.01*(xmax-xmin) xmax+0.01*(xmax-xmin) ymin-0.01*(ymax-ymin) ymax+0.01*(ymax-ymin)])
        
     
        vtrue_stress_xx = reshape(true_stress_xx',N,1);
        vtrue_stress_yy = reshape(true_stress_yy',N,1);
        vtrue_stress_xy = reshape(true_stress_xy',N,1);

        % Figure comparing the true and inferred stress tensor fields
        figure(4000+k0) 
        subplot(3,1,1) 
        hold on
        plot(vtrue_stress_xx,vsigma_post_xx,'b+','MarkerSize',2)
        plot([min(vtrue_stress_xx) max(vtrue_stress_xx)],[min(vtrue_stress_xx) max(vtrue_stress_xx)],'r') %linear case
        if noise_value~=0
            errorbar(vtrue_stress_xx,vsigma_post_xx,reshape(Spost_xx',N,1).^0.5,'.b','MarkerSize',2)
        end
        xlabel('\sigma_{xx}^{num} (Pa \mum)', 'Fontsize', 18)
        ylabel('\sigma_{xx}^{inf}', 'Fontsize', 18)
        
        
        subplot(3,1,2) 
        hold on
        plot(vtrue_stress_yy,vsigma_post_yy,'b+','MarkerSize',2)
        plot([min(vtrue_stress_yy) max(vtrue_stress_yy)],[min(vtrue_stress_yy) max(vtrue_stress_yy)],'r') %linear case
        if noise_value~=0
            errorbar(vtrue_stress_yy,vsigma_post_yy,reshape(Spost_yy',N,1).^0.5,'.b','MarkerSize',2)
        end
        xlabel('\sigma_{yy}^{num} (Pa \mum)', 'Fontsize', 18)
        ylabel('\sigma_{yy}^{inf}', 'Fontsize', 18)
        
        subplot(3,1,3) 
        hold on
        plot(vtrue_stress_xy,vsigma_post_xy,'b+','MarkerSize',2)
        plot([min(vtrue_stress_xy) max(vtrue_stress_xy)],[min(vtrue_stress_xy) max(vtrue_stress_xy)],'r') %linear case
        if noise_value~=0
            errorbar(vtrue_stress_xy,vsigma_post_xy,reshape(Spost_xy',N,1).^0.5,'.b','MarkerSize',2)
        end
        xlabel('\sigma_{xy}^{num} (Pa \mum)', 'Fontsize', 18)
        ylabel('\sigma_{xy}^{inf}', 'Fontsize', 18)
        
    end
end