% Load files
addpath('utils/hyper_pro/');
dir_data = 'lab7/data/';

% Load black files : HyperPro in buoy mode with the black float and the cone
% Group 1
[wv_Es_1, Es_1] = hp_import([dir_data 'HP_float_black_1-HSE187k.dat'], [dir_data 'HP_float_black_1-HED187k.dat'], 'ES');
[wv_Lu_1, Lu_1] = hp_import([dir_data 'HP_float_black_1-HPL174g.dat'], [dir_data 'HP_float_black_1-PLD174g.dat'], 'LU');

% Group 2
[wv_Es_2, Es_2] = hp_import([dir_data 'HP_float_black_2-HSE187k.dat'], [dir_data 'HP_float_black_2-HED187k.dat'], 'ES');
[wv_Lu_2, Lu_2] = hp_import([dir_data 'HP_float_black_2-HPL174g.dat'], [dir_data 'HP_float_black_2-PLD174g.dat'], 'LU');

% Group 3
[wv_Es_3, Es_3] = hp_import([dir_data 'HP_buoy_black_3-HSE187k.dat'], [dir_data 'HP_buoy_black_3-HED187k.dat'], 'ES');
[wv_Lu_3, Lu_3] = hp_import([dir_data 'HP_buoy_black_3-HPL174g.dat'], [dir_data 'HP_buoy_black_3-PLD174g.dat'], 'LU');

figure(1); clf(1, 'reset');
plot(wv_Es_1, Es_1, wv_Es_2, Es_2, wv_Es_3, Es_3);
legend('Group 1', 'Group 2', 'Group 3');
title('HyperPro in buoy mode with the cone');
xlabel('Wavelength (nm)'); ylabel('Es (uW/cm^2/nm)');
xlim([350 800]);
figure(2); clf(2, 'reset'); 
plot(wv_Lu_1, Lu_1, wv_Lu_2, Lu_2, wv_Lu_3, Lu_3);
legend('Group 1', 'Group 2', 'Group 3');
title('HyperPro in buoy mode with the cone');
xlabel('Wavelength (nm)'); ylabel('Lu (uW/cm^2/nm)');
xlim([350 800]);

% Interpolate wavelength
Es_1 = interp1(wv_Es_1, Es_1, wv_Lu_1, 'linear', 'extrap');
Es_2 = interp1(wv_Es_2, Es_2, wv_Lu_2, 'linear', 'extrap');
Es_3 = interp1(wv_Es_3, Es_3, wv_Lu_3, 'linear', 'extrap');

%% Compute Rrs
Rrs_1 = Lu_1 ./ Es_1;
Rrs_2 = Lu_2 ./ Es_2;
Rrs_3 = Lu_3 ./ Es_3;

figure(3); clf(3, 'reset');
plot(wv_Lu_1, Rrs_1,'r', wv_Lu_2, Rrs_2,'c', wv_Lu_3, Rrs_3,'b');
legend('Group 1', 'Group 2', 'Group 3');
title('Lee method for Rrs');
xlim([350 800]);
xlabel('Wavelength (nm)'); ylabel('Rrs (sr^{-1})');

%% Compute SeaWIFS
wide = [20, 20, 20, 20, 20, 20, 40]; % nm
centered = [412, 443, 490, 510, 555, 670, 765]; % nm
i = {};
SeaWIFS_Rrs_1 = []; SeaWIFS_Rrs_2 = []; SeaWIFS_Rrs_3 = [];
for i=1:size(centered,2);
  j = find(centered(i) - wide(i) / 2 <= wv_Lu_1 & wv_Lu_1 <= centered(i) + wide(i) / 2);
  SeaWIFS_Rrs_1(i) = mean(Rrs_1(j));
  j = find(centered(i) - wide(i) / 2 <= wv_Lu_2 & wv_Lu_2 <= centered(i) + wide(i) / 2);
  SeaWIFS_Rrs_2(i) = mean(Rrs_2(j));
  j = find(centered(i) - wide(i) / 2 <= wv_Lu_3 & wv_Lu_3 <= centered(i) + wide(i) / 2);
  SeaWIFS_Rrs_3(i) = mean(Rrs_3(j));
end;

hold('on');
hl=plot(centered, SeaWIFS_Rrs_1, 'ro', centered, SeaWIFS_Rrs_2, 'co', centered, SeaWIFS_Rrs_3, 'ob');
title('Lee method for Rrs', 'FontSize', 18);
xlabel('Wavelength (nm)'); ylabel('Rrs (sr^{-1})');
legend('Group 1', 'Group 2', 'Group 3','SeaWIFS Bands');
set(hl, 'LineStyle', 'none');
set(gca, 'FontSize', 18);

%% Save data
return
info = 'Corrected for dark, Rrs computed according to Lee et al. (2015)';
save('lab7/lab7_black', 'info', 'Es_1', 'Es_2', 'Es_3', 'Lu_1', 'Lu_2', 'Lu_3', 'Rrs_1', 'Rrs_2', 'Rrs_3'); 

%% Save figures
% with eps