clear all,close all,clc;
addpath(genpath('C:\Users\Artun\Documents\MATLAB\YALMIP-master'));
%%



A_11=zeros(3,3);
A_12=eye(3);
A_21=[-2,1,0;
      1,-2,1;
      0,1,-1];
A_22=[-1.6,0.8,0;
      0.8,-1.6,0.8;
      0,0.8,-0.8];
A=[A_11,A_12;A_21,A_22];

B=zeros(6,1); B(4)=1;
%% check the stability using YALMIP
yalmip('clear');
P=sdpvar(6);
F=[P>=0.001*eye(6)];
F=[F;A'*P+P*A<=0];
optimize(F);
Pfeasible = value(P);

% constraint-1: P>0 --> P>=eps*I
% constraint-2: A'P+PA<0

% eig(Pfeasible)
% eig(A'*Pfeasible+Pfeasible*A)
%%
Q=eye(6);
vQ=-Q(:);
tempP=inv([kron(eye(6),A')+kron(A',eye(6))])*vQ;

P=reshape(tempP,6,6);

% constraint-1: P>0
% constraint-2: A'P+PA<0

% eig(P)
% eig(A'*P+P*A)






%%
yalmip('clear');
P=sdpvar(6);
Z=sdpvar(1,6,'full');
F=[P>=0.001*eye(6);
   (A*P)+(A*P)'+(B*Z)+(B*Z)'<=-0.001*eye(6)];
optimize(F);
Pfeasible = value(P);

P=double(P);
Z=double(Z);
K=Z*inv(P);

eig(A+B*K)




%% FMINCON stability test lyapunov
% % % % P>=0
% % % % A'P+PA<0
% % % 
% % % obj1 = @(x)[1];
% % % % % options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% % % % % options = optimoptions('fmincon','Display','iter','Algorithm','trust-region-reflective');
% % % % % options = optimoptions('fmincon','Display','iter','Algorithm','sqp-legacy');
% % % % % options = optimoptions('fmincon','Display','iter','Algorithm','active-set');
% % % % % options = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
% % % options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunctionEvaluations',[1e4]);
% % % x1 = fmincon(obj1,zeros(6,6),[],[],[],[],[],[],@(x)nonlcon(x,A),options)
% % % 
% % % P=x1;
% % % 
% % % [~,flag] = chol(P);
% % % if flag==0
% % %     disp('SYM-POS-DEF-->GOOD!');
% % % else
% % %     disp('NOT SYM-POS-DEF-->BAD!');
% % % end
% % % 
% % % [~,flag] = chol(-(A'*P+P*A));
% % % if flag==0
% % %     disp('SYM-POS-DEF-->GOOD!');
% % % else
% % %     disp('NOT SYM-POS-DEF-->BAD!');
% % % end
% % % 
% % % % constraint function
% % % function [c,ceq] = nonlcon(x,A)
% % % small_number=1e-3;
% % % P=x;
% % % c=[];
% % % c=[c;
% % %    0.0001-eig(P);          % <=0         % P>=0
% % %    0.0001+eig(A'*P+P*A);
% % %    0.0001-trace(P)];   % <=0         % A'P+PA<0
% % % ceq=[ norm(P-P','fro');
% % %       vec(P)-vec(P')];   % ==0
% % % end
%%



%% CVX version
% % % clear all,close all,clc;
% % % 
% % % 
% % % 
% % % 
% % % A_11=zeros(3,3);
% % % A_12=eye(3);
% % % A_21=[-2,1,0;
% % %       1,-2,1;
% % %       0,1,-1];
% % % A_22=[-1.6,0.8,0;
% % %       0.8,-1.6,0.8;
% % %       0,0.8,-0.8];
% % % A=[A_11,A_12;A_21,A_22];
% % % B=zeros(6,1); B(4)=1;
% % % 
% % % 
% % % 
% % % cvx_begin
% % %     variable P(6,6)  
% % %     minimize(1)
% % %     subject to
% % %         -eye(6)-(A'*P+P*A) == semidefinite(6)
% % %         %trace(P) >= 1
% % %         P == semidefinite(6)
% % % cvx_end
% % % 
% % % eig(P)
% % % eig(A'*P+P*A)
%%

















%