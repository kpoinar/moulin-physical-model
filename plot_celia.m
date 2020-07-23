t = res{1,2}.t;
dt = res{1,2}.dt;
Vadd_E = res{1,2}.Vadd_E;
Vadd_C = res{1,2}.Vadd_C;
Vadd_T = res{1,2}.Vadd_turb;
Vadd_OC = res{1,2}.Vadd_oc;
Q_E = Vadd_E;
Q_C = Vadd_C;
Q_T = Vadd_T;
Q_OC = Vadd_OC;
Qin = res{1,2}.Qin;
Qin_comp = res{1,2}.Qin_compensated;

percent_E = Q_E./Qin;
percent_C = Q_C./Qin;
percent_T = Q_T./Qin;
percent_OC = Q_OC./Qin;


clf
figure(20)
%plot(res{1,2}.t,res{1,2}.Qin)
plot(t,Vadd_E)
hold on
plot(t,Vadd_C)
plot(t,Vadd_T)
plot(t,Vadd_OC)
%plot(res{1,2}.t,res{1,2}.Vadd_p)
%plot(res{1,2}.t,res{1,2}.Qin_compensated)
%legend('Qin','Qvadd')
legend('Vadd E','Vadd C','Vadd turb','Vadd oc')%,'Vadd p')



figure(21)
plot(t,percent_E,t,percent_C,t,percent_T,t,percent_OC)
legend('E','C','turb','oc')%,'Vadd p')
title('Percent of Qin. (Vadd/dt /Qin)')

figure(22)
plot(t,Qin,t,Qin_comp)
legend('Qin','Qin_comp')