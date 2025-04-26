%                               Amir Asgharnia--------Mahdi Aliverdinia
%        Numerical study of hydrodynamic stability of Blasius boundary layer flow using Spectral method
%%_________________________________________________________________________________________________________________________
clc; 
clear;
close all;
%%_________________________explanation for Blasius equation________________

%The equation we wish to solve is f''' + (1/2)*f*f'' with f(0) = 0, f'(0) = 0,
%f'(inf) = 1. We recast this problem as a system of first-order ODEs: y =
%[f; f'; f''] = [y(1); y(2); y(3)] so that dy/dEta = y' = [f'; f''; f'''] =
%[y(2); y(3); -(1/2)*y(1)*y(3)] with y(1)(0) = 0, y(2)(0) = 0, y(2)(inf) = 1.

%%_______(1)__________________caling Shooting method and ploting result____

[x,y] = shooting;
figure; hold on;
plot(y(:,1),x,'k-','Linewidth',2)
plot(y(:,2),x,'r-','Linewidth',2)
plot(y(:,3),x,'b-','Linewidth',2)
figformat

%%_______(2)__________________BLASIUS numerical solotion output____________

f=y(:,1);
fp=y(:,2);
fpp=y(:,3);
eta=x;
fppp=diff(fpp)./diff(eta);

%%_______(3)__________________ parameteres input for  stability of Blasius study_______

N=20;
eta_inf=9.55;
mm=(0:N+1);
alfaa=0.01:0.001:0.5;
S=2/eta_inf;
Ree=100:5:3000;
A=zeros(N+2,N+2);
B=zeros(N+2,N+2);

%%_______(4)__________________inline functions for T_______________________  

T=@(m,zita)cos(m*acos(zita));
Tp=@(m,zita)S*((m*sin(m*acos(zita)))/(1 - zita^2)^(1/2));
Tpp=@(m,zita)(S^2)*(((m^2*cos(m*acos(zita)))/(zita^2 - 1) + (m*zita*sin(m*acos(zita)))/(1 - zita^2)^(3/2)));
Tppp=@(m,zita)(S^3)*((m*sin(m*acos(zita)))/(1 - zita^2)^(3/2) - (m^3*sin(m*acos(zita)))/(1 - zita^2)^(3/2) - (3*m^2*zita*cos(m*acos(zita)))/(zita^2 - 1)^2 + (3*m*zita^2*sin(m*acos(zita)))/(1 - zita^2)^(5/2));
Tpppp=@(m,zita)(S^4)*((m^4*cos(m*acos(zita)))/(zita^2 - 1)^2 - (4*m^2*cos(m*acos(zita)))/(zita^2 - 1)^2 + (9*m*zita*sin(m*acos(zita)))/(1 - zita^2)^(5/2) + (15*m^2*zita^2*cos(m*acos(zita)))/(zita^2 - 1)^3 - (6*m^3*zita*sin(m*acos(zita)))/(1 - zita^2)^(5/2) + (15*m*zita^3*sin(m*acos(zita)))/(1 - zita^2)^(7/2));

%%_______(5)__________________making zita vector___________________________ 

zitaa=zeros(1,N);
for J=1:N
zitaa(J)=cos(((J-1)*pi)/(N-1));
end

%%_______(6)__________________making dfzita and df3zita____________________

eta_mo=(eta_inf*(zitaa+1))/2;
dfzita=spline(eta,fp,eta_mo);  
df3zita=spline(eta(1:end-1),fppp,eta_mo);

%%_______(7)__________________boundry for Matrix A_________________________

j=1;
for m=mm
     A(1,j)=1;
     A(2,j)=S*m^2;
     A(N+1,j)=cos(pi*m);
      A(N+2,j)=S*(-m^2*cos(pi*m));
      j=j+1;
    
end

%%_______(8)__________________making A and B Matrix________________________

Re_span=[];
alfa_span=[];
cimax2200=[];
cimax520= [];
for Re=Ree
 
    
 for alfa=alfaa
        
        
%-------(8.1)---------------------making A Matrix(internal nodes)---------
i=3;
j=1;
for e=2:N-1
    
   for m=mm
       coff=(-alfa^2*dfzita(e)-df3zita(e)-((alfa^4)/(1i*alfa*Re)))*T(m,zitaa(e))+((dfzita(e)+((2*alfa^2)/(1i*alfa*Re))))*Tpp(m,zitaa(e))-(1/(1i*alfa*Re))*Tpppp(m,zitaa(e));
       A(i,j)=coff;
       j=j+1;
   end
   j=1;
   i=i+1;
end
%-------(8.2)---------------------making B Matrix(internal nodes)---------
i=3;
j=1;
 for zita=zitaa(2:end-1)
   for m=mm
       coff=(Tpp(m,zita))-((alfa^2)*T(m,zita));
       B(i,j)=coff;
       j=j+1;
   end
   j=1;
   i=i+1;
 end

%---------------------------finding Eigen Values--------------------------
C=eig(A,B);

%-------(8.3)-----------------filtering-----------------------------------
outlierC=isoutlier(imag(C));
CO=[];
for i=1:N
    if outlierC(i)==0
        CO(end+1)=imag(C(i));
    end
    
end
%-------(8.4)-----------------finding max ci-------------------------------
 cimax=max(CO);
%-------------------------------for c and e part in projrct    ------------
     if Re==2200
        cimax2200(end+1)=cimax;
     end
     %%%
      if Re==520
        cimax520(end+1)=cimax;
     end
     
%----------------------finding roots---------------------------------------
 
 if (cimax>0 && cimax<0.0001)
     Re_span(end+1)=Re;
     alfa_span(end+1)=alfa;

 end
%--------------------------------------------------------------------------
             
 end %end of loop alfaa
      
end %end of loop Ree
%___________________ploting result_________________________________________
figure('name','neutral stability Curves','NumberTitle','off');
scatter(Re_span,alfa_span,'LineWidth',2);
xlabel('Re');
ylabel('\alpha');
grid on;
%--------------------------------------------------------------------------
figure('name','Re = 2200','NumberTitle','off');
scatter(alfaa, cimax2200,'LineWidth',2);
xlabel('\alpha');
ylabel('ci-max at reynolds = 2200');
grid on;
%--------------------------------------------------------------------------
figure('name','Re critical = 520','NumberTitle','off');
scatter(alfaa, cimax520,'LineWidth',2);
xlabel('\alpha');
ylabel('ci-max at critical reynolds = 520');
grid on;
%--------------------------------------------------------------------------


%__________________________________________________________________________
%% Supplementary Functions

function [x,y] = shooting
% Use fsolve to ensure the boundary function is zero. The result is the
% unknown initial condition.
opt = optimset('Display','off','TolFun',1E-20);
F = fsolve(@(F) eval_boundary(F),[0,0,0.33],opt);

% Solve the ODE-IVP with the converged initial condition
[x,y] = solve_ode(F);
end

function [x,y] = solve_ode(F)
% Solve the ODE-IVP with initial condition F on [0 100] (arbitrary upper
% bound)
[x,y] = ode45(@(x,y) [y(2); y(3); -0.5*y(1)*y(3)],[0 9.55],F); %solve BVP                
end

function [g] = eval_boundary(F)
% Get the solution to the ODE with inital condition F
[x,y] = solve_ode(F);

% Get the function values (for BCs) at the starting/end points
f_start = y(1,1); %f(0) = 0
df_start = y(1,2); %f'(0) = 0
df_end = y(end,2); %f'(inf) - 1 = 0

% Evaluate the boundary function
g = [f_start
     df_start
     df_end - 1];
end

function figformat
% This function simply formats the figure
xlabel('\it{f, f^{(1)}, f^{(2)}}');
ylabel('\eta');
xlim([0 2]);
xticks(0:0.5:2);
ylim([0 10]);

h = legend('\it{f}','\it{f}^{(1)}','\it{f^{(2)}}','Location','NorthEast');
legend boxoff;
set(h,'FontSize',18);

axis square
set(gca, 'box', 'on'); %creates a box border
set(gca,'FontWeight','normal','linewidth',1.05); fontname = 'Arial'; %the next few lines format the plot
set(0,'defaultaxesfontname',fontname); set(0,'defaulttextfontname',fontname);
fontsize = 20;
set(0,'defaultaxesfontsize',fontsize); set(0,'defaulttextfontsize',fontsize);
end