function [FEA_RF_T] = get_reactfrcFEA(FEA_RF_ADD)
%----------------------------------------------------------------------%
%---------      Extract reaction force from rigid wall-----------------%
%----------------------------------------------------------------------%

fid=fopen(FEA_RF_ADD); % FEA_RF_ADD is the reaction force output file address
Rwforce_Matlab = zeros(101,6); % 101 is the time step number from LS-Dyna ouput
Num = 0;
while 1
    Num = Num + 1;    
    tline = fgetl(fid);
    if ~ischar(tline )
        break;
    end    
    if Num > 11
       Rwforce_Matlab(Num-11,:) = str2num(tline);
    end    
end
fclose(fid);

FEA_RF = zeros(101,1);
for i = 1:101
    FEA_RF(i,1) = 4* Rwforce_Matlab(i,3);% Since it is 1/8 model so we times 4
end

FEA_RF_T = FEA_RF';
%----------------------------------------------------------------------%
end