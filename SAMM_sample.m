%% === load sample data ===
load('vel_sample.mat'); % velocity
load('kulite1_sample.mat'); % kulite 1
load('kulite2_sample.mat'); % kulite 2
Kulite(1).pp = kulite1; % structure
Kulite(2).pp = kulite2;

%% === Run SAMM ===
f_range = [33 76 121 170]; % NFFT = 2048
nfft = 2048;
fs = 25600;
method = 'RR';
blocks = 100; % blocks for SMAMM-SPOD (reduce this number to speed up)
[freq,eigval,modes] = SAMM_methods(vel,Kulite,nfft,fs,f_range,method);
%%
[nx,ny] = size(X);
dx = X(1,2)-X(1,1);
figure('position',[0 800 1100 135])
set(gcf, 'DefaultAxesFontName', 'Times New Roman');
set(gcf, 'DefaultTextFontName', 'Times New Roman');
count = 1;

for fi = 1:length(f_range)

    if strcmp(method,'RR')
    phi_temp = modes(1:nx*ny,1,fi);
       elseif strcmp(method,'SPOD')
    phi_temp = modes(1:nx*ny,fi);
    end
    phi_temp = phi_temp./norm(phi_temp);

    ang_ref = 0;
    vel_amp = abs(reshape(phi_temp,[nx,ny]))/dx;
    vel_ang = angle(reshape(phi_temp,[nx,ny]));
    vel_real = vel_amp.*(cos(vel_ang+(ang_ref-1)*pi/45));

    subplot(1,4,count)
    contourf(X,Y,vel_real,50,'linestyle','none');
    colormap(jet)
    caxis([-0.5 0.5])
    axis image
    if fi == length(f_range)
    colorbar 
    end
    hold on
    xlim([0 6])
    set(gca,'fontsize',16);
    xlabel('$x/D$','fontsize',16,'interpreter','latex');
    if fi == 1
    ylabel('$y/D$','fontsize',16,'interpreter','latex');
    else
        set(gca,'YTick',[]);
    end
    set(gcf,'renderer','Painters')
    if fi == 1
        set(gca,'position',[0.049 0.255 0.21 0.88]);
    elseif fi == 2
        set(gca,'position',[0.275 0.255 0.21 0.88]);
    elseif fi == 3
        set(gca,'position',[0.501 0.255 0.21 0.88]);
    else
        set(gca,'position',[0.727 0.255 0.21 0.88]);
    end
    count = count+1;
end
