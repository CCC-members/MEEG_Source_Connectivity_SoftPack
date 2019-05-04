function [] = figure_occ_harmonics(F_task,PSD_task,freq)
figure; 
set(gcf,'Position',[50 50 1400 700]);
plot(F_task,PSD_task([9 10],:)); xlabel('frequency'); ylabel('PSD'); legend('O1','O2'); title('power-spectral-density task-segments');
[tmp,subh1]   = min(abs(F_task-(1/2)*freq));
[tmp,mainh]   = min(abs(F_task-freq));
[tmp,superh1] = min(abs(F_task-2*freq));
txt           = ['\leftarrow','1/2','X',num2str(freq),'Hz'];
text(F_task(subh1),PSD_task(9,subh1),txt);
txt           = ['\leftarrow',num2str(freq),'Hz','(main-harmonic)'];
text(F_task(mainh),PSD_task(9,mainh),txt);
txt           = ['\leftarrow','2','X',num2str(freq),'Hz'];
text(F_task(superh1),PSD_task(9,superh1),txt);
end
