function ADCS_save(SIM,code_source)


if code_source == 1
    str = 'm';
else
    str = 'c';
end

c = 1;

for i = 1:22
    saveas(figure(i),['plots_',str,'/fig',num2str(i),'_',str],'jpg')
end

disp('Plots saved')