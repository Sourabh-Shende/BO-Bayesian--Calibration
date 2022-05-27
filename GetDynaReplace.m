function [Line_ID] = GetDynaReplace1(K_Address, cm1,cm4,cm6)

%----------------------------------------------------------------------%
%---------   Update C10 and C20 values for FEA input file  ------------%
%----------------------------------------------------------------------%  


fid = fopen(K_Address,'r');        % Open File to read
cm2 = 0;
cm3 = 0;
cm5 = 0;

C = sprintf('%10.1f %9.3e %9.1f %9.1f %9.1f %9.1f',cm1,cm2,cm3,cm4,cm5,cm6)
replaceline = C;

i = 1;
tline = 's';
A = {[]};
while ischar(tline)
    tline = fgetl(fid);
    if ~isempty(strfind(tline,'         2    1000.0     0.495         0 ')) % Locate to the key position which is above the C10 and C20
        A{i}=tline;
        A{i+1} = replaceline;               % Replacce the old value with new value
        tline = fgetl(fid);
        i = i+1;
        i;
    else
        A{i}=tline;
    end
    i = i+1;
end
Line_ID = i;
fclose(fid);
fid2=fopen(K_Address,'w');            % Open file to write
for i=1:length(A)-1
    fprintf(fid2,['%s',char([13,10])],A{i});
end
fclose(fid2); % close the file
%----------------------------------------------------------------------%  

end

