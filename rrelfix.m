%RRELFIX Read output file from relfix, an exe-file
%	 in the DUT GPS s/w package

% Written by Kai Borre
% October 28, 2001

fid = fopen('dut.tex','rt');
line = fgetl(fid); % reads ---
line = fgetl(fid); % reads headings
i = 0;

while feof(fid) == 0
    i = i+1;
    line = fgetl(fid);
    [n,r] = strtok(line,',');
    for j = 1:5
	[n,r] = strtok(r,',');
    end
    lent(i) = str2num(n);
    [n,r] = strtok(r,',');
    east(i) = str2num(n);
    [n,r] = strtok(r,',');
    north(i) = str2num(n);
    [n,r] = strtok(r,',');
    up(i) = str2num(n);
    line = fgetl(fid);	% reads ----
end
plot(1:i,[east-east(1); north-north(1); up-up(1)]*1000,'linewidth',2)
title('Differential Position Estimates from RELFIX')
xlabel('Epoch [1 s interval]')
ylabel('Corrections to Initial Position [mm]')
legend('Easting','Northing','Upping')
fclose(fid);
print -depsc relfix
%%%%%%%%%%%%%%%%%%%%%%% end rrelfix.m  %%%%%%%%%%%%%%%%%%%%%%%


