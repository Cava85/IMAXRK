% function [k] = consraints_delta_dogs(Var, Method, fun, con,lob,upb)
% The Delta_Dogs algorithm with nonlinear costraints.
% Simple bound costraints for the control parameter is imposed.
% bnd1<x<bnd2

% fun: The objective function
% con{i}= The costraint function

% n:dimension of the problem, ms: number of state constraints, m: number of
% control constraints

% Search.Method=1 constant_K
% Search.Method=2 Adaptive_K
clear all
close all
clc
% run('../addTopath.m')
global n m ms bnd1 bnd2 Ain bin acon Search r xi tri
r=2;
% parameters
n=2; ms=1; m=2*n; 
% iter_max=10;
iter_max=18;
RES = 1e-2;
% RES = 1e-3;
x_star=[0.7;0.1];
% if nargin <2 
% Search.method=2;
% Search.constant=0;
%%
Method = 2;
 Var = 0.25;
% Method =1;
%  Var = 1;

Ex=1;
%%
if Ex==1
b=1;
% specify the objective function
% % % % % % % % % % % % % % % fun =@(x)x(2,:);
fun=@(x) (x(1,:)).^2+(x(2,:)).^2;
%%
% specify the contraints functions
 con{1}=@(x) (-(x(2,:)-b*x(1,:).^2 - 0.5));
%  con{2}=@(x) ((x(2,:)-b*x(1,:).^2 - 0.5));
end
    lob=zeros(n,1); upb=ones(n,1);

Search.method = Method;
Search.constant = Var;
bnd1 = lob;
bnd2 = upb;
% Search.method=1;
% Search.constant=1;
% interpolaion strategy
inter_method=1;
% upper and lower bouds
% Input the equality constraints
Ain=[eye(n);-eye(n)];
bin=[bnd2 ;-bnd1];
% Calculate the initial points
  xi=bounds(bnd1,bnd2, n);
% Calculate the function evaluation at initial points
for ii=1:m
    acon{ii}=[];
end
for ii=1:size(xi,2)
    ind=1:2*n;
    ind=ind(Ain*xi(:,ii)-bin>-0.01);
    for jj=1:length(ind)
    acon{ind(jj)}=[acon{ind(jj)} ii];
    end
    yi(ii)=fun(xi(:,ii));
    for jj=1:ms
        C{jj}(ii)=con{jj}(xi(:,ii));
    end
end
x_prev = xi(:,1);
y_prev = yi(:,1);
delta_tol=0.2;

xi(1,1)=xi(1,1)+0.00001;
        for k=1:iter_max
            tri=delaunayn(xi.');
            
      inter_par_p= interpolateparametarization(xi(:,C{1}-RES<0),yi(C{1}-RES<0),inter_method);
%         inter_par_p= interpolateparametarization(xi,yi,inter_method);
  for jj=1:ms
      inter_par_g{jj}= interpolateparametarization(xi,C{jj},inter_method);
  end 
%   keyboard
            [xm ym(k) Cs(k)]= tringulation_search_constraints(inter_par_p,inter_par_g,xi,tri);
           % [y_star(k),indm]=search_fun_val(x_star, inter_par_p,inter_par_g,xi, tri);
            % feasible constraint projection
% % % % %            [xm]=feasible_constraint_box(xm,xi,length(xi)+1);
            xi=[xi xm];
            yi=[yi fun(xm)];
            for jj=1:ms
               C{jj}=[C{jj} con{jj}(xm)];
            end
%%
    tt = 0:0.01:1;
      impro=mindis(xm,xi(:,1:end-1));
if impro<RES|| k == iter_max +1
    disp('it converged ... at step    ')
    disp(k)
    iter_1 = k;
    rootpath = './figures/';
     fpath = strcat(rootpath,'/Method_',num2str(Method), '_const_', num2str(Var*100),'_example');
     fName = strcat( fpath ,'/Method_',num2str(Method), '_const_', num2str(Var*100), 'ex_iterations', num2str(33));
    % keyboard
%      save(strcat(fName,'_workspace','.mat')) 
%   csvwrite(strcat(fName,'.csv'), iter_1);
   break
end
figure(1)
scatter(xi(1,:), xi(2,:))
%% plot the feasible domain at each iteration
for Iv = 1
clear U_f p_f G_c C_l Err TT SS G1 G2
xv=0:0.01:1; 

for ii=1:length(xv)
    for jj=1:length(xv)
%         U_f(ii,jj)=fun([xv(ii) ;xv(jj)]);
         P_f(ii,jj)=interpolate_val([xv(ii) ;xv(jj)],inter_par_p);
         C_l(ii,jj)=con{1}([xv(ii) ;xv(jj)]);
%         G_c(ii,jj)=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{1});
% keyboard
        G1=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{1});
%         G2=interpolate_val([xv(ii) ;xv(jj)],inter_par_g{2});
       [Err(ii,jj),TT(ii,jj),SS(ii,jj)] = direct_uncer([xv(ii) ;xv(jj)],xi,inter_par_p,inter_par_g,tri);
%        CC(ii,jj) = max([G1, G2])- Search.constant * Err(ii,jj);
            CC(ii,jj) = G1- Search.constant * Err(ii,jj);
       if CC(ii,jj) >0
           CC(ii,jj) = nan;
       end
    end
end


h=figure(13); clf;
set(h,'units','normalized','outerposition',[0 0 0.9 0.9])
%default values
op.fontsize = 12;
op.quiet = 0;
% op.dirpath = conf.outputDir;
op.dirpath = './';
op.width = 1500;
op.height = 900;
op.visible = 1;
op.title = 1;
op.cache = 0;
op.output = 0;
op.embed = 0;
op.save = 0;
op.label = 1; %Settign to zero removes labels on figures
ah(1) = axes('Position', [0.015 0.48 0.3 0.47]);
contourf(Err.')
colorbar
title('e(x)')
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		ah(2) = axes('Position', [0.34 0.5 0.33 0.47]);
        ah(2) = axes('Position', [0.34 0.48 0.30 0.47]);
        hold on
  contourf(xv,xv,(CC).')
   plot(xi(1,:),xi(2,:),'sr') 
   plot(xi(1,C{1}<0),xi(2,C{1}<0),'kx')
caxis([-0.5,0])
colorbar
title('g(x) - K e(x)')
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
% ah(3) = axes('Position', [0.675 0.5 0.33 0.47]);
ah(3) = axes('Position', [0.675 0.48 0.30 0.47]);
contourf(xv,xv,(P_f-Var*Err).')
caxis([0,1])
colorbar
title('p(x) - K e(x)')
axis off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		ah(4) = axes('Position', [0.045 0.015 0.36 0.48]);
        	ah(4) = axes('Position', [0.045 0.05 0.45 0.38]);
    plot(1:length(yi(1:end)),yi,'-')
    title('objective function')
    grid on
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		ah(5) = axes('Position', [0.45 0.05 0.53 0.41]);
        	ah(5) = axes('Position', [0.55 0.05 0.38 0.38]);
        hold on
    triplot(tri,xi(1,:),xi(2,:)) 
    if Ex ==1
%        plot(tt, b*sin(pi*tt), 'r', 'linewidth', 3)
    elseif Ex==2
       plot(tt, rastriginn(tt), 'r', 'linewidth', 3)  
    elseif Ex ==3
        
       tri=delaunayn(xi.');
       hold on
       tt=0:0.01:1;
for ii=1:length(tt)
    for jj=1:length(tt)
        U(ii,jj)=rastriginn2([tt(jj) ;tt(ii)]);
    end
end
contourf(tt,tt,U,0:0.5:0.5)
colormap gray
% colormap bone
brighten(0.5)
    end
        triplot(tri,xi(1,:),xi(2,:))
    plot(xi(1,:),xi(2,:),'sr') 
    hold off
      if Method ==1
    title(strcat('Position of the points',' for K = ', num2str(Search.constant)))
    elseif Method ==2
   title(strcat('Position of the points',' for y_0 = ', num2str(Search.constant) ))
    end
    xlim([0 1])
    ylim([0 1])
    
		set(ah, 'FontSize', op.fontsize);		        
%        F(k) = getframe(h); 
%         saveas(h,strcat('./images/', num2str(k,'%02d')), 'png')
%        imagerDir = './images/';
%        videoDir = '.'
% 	   mp4filename = sprintf('%s/Videos/MP4/%d.mp4', imagerDir,k);
% 		encodeFramesToVideo(videoDir, 'png', mp4filename);


end

drawnow
disp( [' finished iteration ' ,num2str(k), ' ... '] )
% keyboard
        end
        
       
        