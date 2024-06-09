d.Value = .5;
d.Message = 'Determining connected cameras';

adaptors=imaqhwinfo;
installed_adaptors = adaptors.InstalledAdaptors;

if length(installed_adaptors)~=0
for i=1:1:length(installed_adaptors)
info = imaqhwinfo(installed_adaptors{1,i});
installed_adaptors_vector(i)=length(info.DeviceInfo);    
end

indx = find(installed_adaptors_vector==1);

disp(['Detected cameras: '])

for i=1:1:length(indx)
disp(installed_adaptors{1,indx(i)})
end

pause(1)
else

disp(['Detected cameras: none']);    

end
%%
d.Value = .75;
d.Message = 'Determining SLM resolution';

Screen_positions = get(0,'MonitorPositions');
No_of_screens = size(Screen_positions,1);

if No_of_screens==1
    disp(['Detected SLMs: 0'])    
else
    No_of_screens=No_of_screens-1;    
    disp(['Detected SLMs: ' num2str(No_of_screens) ' Resolution: ' num2str(Screen_positions(2,3)) 'x' num2str(Screen_positions(2,4))])    
end

pause(1)