ref_data = csvread('figure/data_fixedfixed/paidoussis.csv',2,0);

figure();
subplot(1,2,1);
hold on;
title('(a)');
xlabel('U_R');
ylabel('\Bar{\omega}_1');
for i=1:3
    vectX = ref_data(:,i*2-1);
    vectY = ref_data(:,i*2);
    vectX(vectX==0) = [];
    vectY(vectY==0) = [];
    plot(vectX,vectY);
end
plot(U_bar,real(EigFreq_bar(1,:)),'k','LineWidth',1);
plot(U_bar,imag(EigFreq_bar(1,:)),'k:','LineWidth',1.2);
plot(U_bar,imag(EigFreq_bar(2,:)),'k:','LineWidth',1.2);

subplot(1,2,2);
hold on;
title('(b)');
ylabel('\Bar{\omega}_1');
xlabel('U_R');
ylabel('\Bar{\omega}_1');
for i=4:5
    vectX = ref_data(:,i*2-1);
    vectY = ref_data(:,i*2);
    vectX(vectX==0) = [];
    vectY(vectY==0) = [];
    plot(vectX,vectY);
end
plot(U_bar,real(EigFreq_bar(3,:)),'k','LineWidth',1);
plot(U_bar,imag(EigFreq_bar(3,:)),'k:','LineWidth',1.2);
plot(U_bar,imag(EigFreq_bar(4,:)),'k:','LineWidth',1.2);