clear all;
close all;
clc;
fov = 1;
k0  = 9703*3600;
H_rxn = 5960;
E_act = 11843;
rhocp = 500;
Tf    = 298;
caf   = 10;
UAoV  = 150;
R     = 1.987;
Tj    = 298;
% al rededor de este punto de operacion:
Temp=311.2;  
ca=8.5636;
%
ks = k0*exp(-E_act/(R*Temp));
ksprime = ks*(E_act/(R*Temp*Temp));
%
A(1,1) = -fov - ks;
A(1,2) = -ca*ksprime;
A(2,1) = ks*H_rxn/rhocp;
A(2,2) = -fov - UAoV/rhocp + (H_rxn/rhocp)*ca*ksprime;

B(1,1)=0;
B(1,2)=0;
B(2,1)=UAoV/rhocp;
B(2,2)=fov;

C=[1 0 ; 0 1];
D=zeros(2,2);
t = 0:0.001:9;
u = zeros(size(t,2),2);
x0 = [3 330];
%%
%%OPENLOOP
sys = ss(A,B,C,D);
[y_ol,t_ol,x_ol] = lsim(sys,u,t,x0);
%%
%OBSERVER1
op1 = -100;
op2 = -102;
L = place(A',C',[op1 op2])';
At = [ A     zeros(size(A))
     L*C  A-L*C ];
Bt = [    B
      B ];
Ct=eye(4);
sys = ss(At,Bt,Ct,0);
[y_ob,t_ob,x_ob]=lsim(sys,zeros(size(t,2),2),t,[x0 x0]); % orden 4
%%
%OBSERVER2
op1 = -100;
op2 = -102;
L = place(A',C',[op1 op2])';
At = [ A     zeros(size(A))
     L*C  A-L*C ];
sys = ss(At,Bt,Ct,0);
[y_ob1,t_ob1,x_ob1]=lsim(sys,zeros(size(t,2),2),t,[x0 -14 0]); % orden 4
%%
%OBSERVER3
op1 = -1;
op2 = -2;
L = place(A',C',[op1 op2])';
At = [ A     zeros(size(A))
     L*C  A-L*C ];
sys = ss(At,Bt,Ct,0);
[y_ob2,t_ob2,x_ob2]=lsim(sys,zeros(size(t,2),2),t,[x0 -14 0]); % orden 4

%%
%OBSERVER
op1 = -1+30i;
op2 = -2-30i;
L = place(A',C',[op1 op2])';
At = [ A     zeros(size(A))
     L*C  A-L*C ];
sys = ss(At,Bt,Ct,0);
[y_ob3,t_ob3,x_ob3]=lsim(sys,zeros(size(t,2),2),t,[x0 -14 0]); % orden 4
%%
%CLOSED LOOP
p1 = -5;
p3 = -6;
K = place(A,B,[p1 p3]);
sys = ss(A-B*K,B,C,0);
[y_cl,t_cl,x_cl]=lsim(sys,zeros(size(t,2),2),t,x0); 

%%
%CLOSED LOOP WITH OBSERVER1
op1 = -100;
op2 = -102;
L = place(A',C',[op1 op2])';
Atclo = [ A-B*K     zeros(size(A))
     L*C  A-L*C ];
 sys = ss(Atclo,Bt,Ct,0);
[y_clo,t_clo,x_clo]=lsim(sys,zeros(size(t,2),2),t,[x0 x0]);
%%
%CLOSED LOOP WITH OBSERVER
op1 = -100;
op2 = -102;
L = place(A',C',[op1 op2])';
Atclo = [ A-B*K     zeros(size(A))
     L*C  A-L*C ];
 sys = ss(Atclo,Bt,Ct,0);
[y_clo1,t_clo1,x_clo1]=lsim(sys,zeros(size(t,2),2),t,[x0 -14 0]);
%%
%%
%CLOSED LOOP WITH OBSERVER
op1 = -1;
op2 = -2;
L = place(A',C',[op1 op2])';
Atclo = [ A-B*K     zeros(size(A))
     L*C  A-L*C ];
 sys = ss(Atclo,Bt,Ct,0);
[y_clo2,t_clo2,x_clo2]=lsim(sys,zeros(size(t,2),2),t,[x0 -14 0]);
%%
%%
%CLOSED LOOP WITH OBSERVER
op1 = -1+30i;
op2 = -2-30i;
L = place(A',C',[op1 op2])';
Atclo = [ A-B*K     zeros(size(A))
     L*C  A-L*C ];
 sys = ss(Atclo,Bt,Ct,0);
[y_clo3,t_clo3,x_clo3]=lsim(sys,zeros(size(t,2),2),t,[x0 -14 0]);
%%
%%
%PLOTS

figure,
subplot(2,1,1),plot(t_ol,[y_ob(:,1),y_ob(:,3),y_ob1(:,3),y_ob2(:,3),y_ob3(:,3)]),legend({'OpenLoop','Best','InCond','BadPole','ImgPole'}),grid,ylabel('Conc[A/B]')
subplot(2,1,2),plot(t_ol,[y_ob(:,2),y_ob(:,4),y_ob1(:,4),y_ob2(:,4),y_ob3(:,4)]),legend({'OpenLoop','Best','InCond','BadPole','ImgPole'}),grid,ylabel('Temp[K]')
xlabel('Tiempo[s]');

figure,
subplot(2,1,1),plot(t_ol,[y_cl(:,1),y_clo(:,3),y_clo1(:,3),y_clo2(:,3),y_clo3(:,3)]),xlim([0 2]),legend({'ClosedLoop','Best','InCond','BadPole','ImgPole'}),grid,ylabel('Conc[A/B]')
subplot(2,1,2),plot(t_ol,[y_cl(:,2),y_clo(:,4),y_clo1(:,4),y_clo2(:,4),y_clo3(:,4)]),xlim([0 2]),legend({'ClosedLoop','Best','InCond','BadPole','ImgPole'}),grid,ylabel('Temp[K]')
xlabel('Tiempo[s]');
