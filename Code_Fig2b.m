%code for Fig2.b
clear all;
clc;
S=load('viotf1tf1.mat');
z=S.VIO;
x=linspace(0,1,1001);
y=linspace(0,1,1001);
imagesc(x,y,z);
axis xy;

caxis([-0.15,0.06]);




colorbar;
xlabel('x');
ylabel('y');

n=256;
z_min=-0.15;
z_max=0.06;
z_zero=0;
zero_pos=(z_zero-z_min)/(z_max-z_min);
n_neg=round(zero_pos*n);
n_pos=n-n_neg+1;



blue=[linspace(0,1,n_neg)',linspace(0,1,n_neg)',ones(n_neg,1)];
red=[ones(n_pos,1),linspace(1,0,n_pos)',linspace(1,0,n_pos)'];
cmap=[blue;red(2:end,:)];
colormap(cmap);
cb=colorbar;
cb.Label.String='-\langle\sigma\rangle';
cb.Label.Rotation=0;

hold on;
contour(x,y,z,[0,0],'k-','LineWidth',1.5);
xlabel('\epsilon');
ylabel('d');
zlabel('-\sigma');
hold on;
plot([0.2 0.2],[0 1],'k--','LineWidth',1.5);
hold on;
plot([0 1],[0.2 0.2],'k--','LineWidth',1.5);
hold off;