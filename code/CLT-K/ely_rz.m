% read the particle data
partdat = './p5383670.dat';
fid = fopen(partdat,'rb');
pdata = fread(fid,'double');
fclose('all');

%%
figure;
ng = 256;
jp = 1;
dr = 3 / (ng - 1);
dz = 3 / (ng - 1);
z0 = -1.5 : dz : 1.5;
r0 = 2.5 : dr : 5.5;
[r, z] = meshgrid(r0, z0);
elytot = zeros(ng);
    
for it = 2 : 1 : length(pdata) / 8
%     figure; 
    ely = zeros(ng);
    n0 = ones(ng);
%     set(gca,'nextplot','replacechildren');
    set(gcf,'DefaultAxesFontSize',15);
%     set(gcf,'Position',get(0,'ScreenSize'));
    set(gcf,'Position',[200 200 800 700]);
    for i = 1 : 1 : 1
        eval(['R = pdata','(1 : 8 : end);']);
        eval(['Z = pdata','(3 : 8 : end);']);
        eval(['timep = pdata','(8 : 8 : end);']);

        if it > length(R)
            continue;
        else
            gbr = floor(abs(R(it - 1) - min(r0))/ dr) + 1; 
            gfr = gbr + 1;
            wbr = 1-abs(abs(R(it - 1) - min(r0)) / dr - gbr);
            wfr = 1 - wbr;
            gbz = floor((Z(it - 1) - min(z0)) / dz) + 1; 
            gfz = gbz + 1;
            wbz = 1-abs((Z(it - 1) - min(z0)) / dz - gbz);
            wfz = 1 - wbz;
            rr = R(it) - R(it - 1);
            rz = Z(it) - Z(it - 1);
            rly = sqrt(rr .* rr + rz .* rz);
            ely0 = max(rly);

            ely(gbr,gbz) = ely(gbr,gbz) + ely0 * wbr * wbz;
            ely(gbr,gfz) = ely(gbr,gfz) + ely0 * wbr * wfz;
            ely(gfr,gbz) = ely(gfr,gbz) + ely0 * wfr * wbz;
            ely(gfr,gfz) = ely(gfr,gfz) + ely0 * wfr * wfz;
            n0(gbr,gbz) = n0(gbr,gbz) + 1;
            n0(gbr,gfz) = n0(gbr,gfz) + 1;
            n0(gfr,gbz) = n0(gbr,gbz) + 1;
            n0(gfr,gfz) = n0(gbr,gbz) + 1;
        end
        elytot = elytot + ely ./ n0 /length(pdata) * 8;
        pcolor(r,z,elytot');shading interp;
        xlim([3, 5]);ylim([-1, 1]);
        xlabel('R','fontsize',18);
        ylabel('Z','fontsize',18);
        title(['Lyapunov exponent t = ',num2str(timep(it - 1))]);
%         caxis([-2,2]);
        colorbar();
        drawnow;
        axis equal;
    end
    Frz(jp) = getframe(gcf);
    jp = jp + 1;
end
