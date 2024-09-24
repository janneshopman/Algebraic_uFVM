function plot2DIncBo(field, Nx, Ny)

plotField = zeros(Nx+2, Ny+2);

Nc = Nx*Ny;

%TGV
% bNorth = field(Nc+1:Nc+Ny);
% bSouth = field(Nc+Nx+1:Nc+2*Nx);
% bWest = field(Nc+2*Nx+1:Nc+2*Nx+Ny);
% bEast = field(Nc+2*Nx+Ny+1:Nc+2*Nx+2*Ny);

%LDC
bNorth = field(Nc+1:Nc+Ny);
bWest = field(Nc+Nx+1:Nc+2*Nx);
bEast = field(Nc+2*Nx+1:Nc+2*Nx+Ny);
bSouth = field(Nc+2*Nx+Ny+1:Nc+2*Nx+2*Ny);

for iy = 1:Ny+2
    for ix = 1:Nx+2
        ic = (iy-1)*(Nx+2)+ix;
        if iy == Ny+2 && ix ~= 1 && ix ~= Nx+2
            plotField(ix, iy) = bNorth(ix-1);
        elseif iy == 1 && ix ~= 1 && ix ~= Nx+2
            plotField(ix, iy) = bSouth(ix-1);
        elseif ix == 1 && iy ~= 1 && iy ~= Ny+2
            plotField(ix, iy) = bWest(iy-1);            
        elseif ix == Nx+2 && iy ~= 1 && iy ~= Ny+2
            plotField(ix, iy) = bEast(iy-1);
        elseif ix ~= 1 && ix ~= Nx+2 && iy ~= 1 && iy ~= Ny+2
            plotField(ix, iy) = field((iy-2)*Nx+(ix-1));
        else
        end
    end
end

imagesc(rot90(plotField));
colorbar;
