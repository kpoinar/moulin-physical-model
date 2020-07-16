clf
figure(20)
%plot(res{1,2}.t,res{1,2}.Qin)
plot(res{1,2}.t,res{1,2}.Vadd_E)
hold on
plot(res{1,2}.t,res{1,2}.Vadd_C)
plot(res{1,2}.t,res{1,2}.Vadd_turb)
plot(res{1,2}.t,res{1,2}.Vadd_oc)
%plot(res{1,2}.t,res{1,2}.Vadd_p)
%plot(res{1,2}.t,res{1,2}.Qin_compensated)
%legend('Qin','Qvadd')
legend('Vadd E','Vadd C','Vadd turb','Vadd oc')%,'Vadd p')



