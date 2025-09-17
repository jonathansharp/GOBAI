idxsurf = predictor_matrix(:,3) == 10;
figure;
histogram2(predictor_matrix(idxsurf,1),predictor_matrix(idxsurf,2));
xlabel('Temperature (degC)'); ylabel('Salinity'); zlabel('Counts');
title('10dbar');

idxmid = predictor_matrix(:,3) == 200;
figure;
histogram2(predictor_matrix(idxmid,1),predictor_matrix(idxmid,2));
xlabel('Temperature (degC)'); ylabel('Salinity'); zlabel('Counts');
title('200dbar');

idxdeep = predictor_matrix(:,3) == 1500;
figure;
histogram2(predictor_matrix(idxdeep,1),predictor_matrix(idxdeep,2));
xlabel('Temperature (degC)'); ylabel('Salinity'); zlabel('Counts');
title('1500dbar');
