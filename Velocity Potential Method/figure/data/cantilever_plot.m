ref_data = csvread('paidoussis.csv',2,0);
pot_data = csvread('potential_method_results.csv',2,0);

figure();
subplot(1,2,1);
hold on;
for i=1:3
    vectX = ref_data(:,i*2-1);
    vectY = ref_data(:,i*2);
    vectX(vectX==0) = [];
    vectY(vectY==0) = [];
    plot(vectX,vectY);
end
for i=1:3
    vectX = pot_data(:,i*2-1);
    vectY = pot_data(:,i*2);
    vectX(vectX==0) = [];
    vectY(vectY==0) = [];
    plot(vectX,vectY);
end

subplot(1,2,2);
hold on;
for i=4:5
    vectX = ref_data(:,i*2-1);
    vectY = ref_data(:,i*2);
    vectX(vectX==0) = [];
    vectY(vectY==0) = [];
    plot(vectX,vectY);
end
for i=4:5
    vectX = pot_data(:,i*2-1);
    vectY = pot_data(:,i*2);
    vectX(vectX==0) = [];
    vectY(vectY==0) = [];
    plot(vectX,vectY);
end