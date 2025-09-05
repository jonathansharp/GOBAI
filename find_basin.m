function basin_id = find_basin(lon,lat)

% convert lon

%Coordinate bounds for each basin with Pacific/Atlantic/Indian overlap in
%Southern Ocean
x_ind = [12.0540, 12.0540, 37.8730, 152.6829, 152.6829, 125.1127, 117.6949, 114.1562, 105.9897, 101.0218, 102.4510, 98.7080, 102, 47.3085, 43.5715, 12.0540];
y_ind = [-18.0231,  -77.9140, -77.9140, -77.6793, -29.0586, -14.7138, -8.5563, -8.0516, -5.8309, 0.0230, 3.4558, 9.3104, 21, 32.6206, 12.9936, -18.0231];
x_pac = [102, 98.7080, 102.4510, 101.0218, 105.9897, 114.1562, 117.6949, 125.1127, 125.0001, 125.0001, 166.5625, 290, 298.4293, 300.5583, 302.6054, 302.6054, 304.4254, 283.7198, 280.8265, 279.2641, 277.4294, 269.3044, 262.2278, 256.4617, 240, 110.7359, 102];
y_pac = [21, 9.3104, 3.4558, 0.023, -5.8309, -8.0516, -8.5563, -14.7138, -31.77, -77.6793, -86.3625, -87.3012,  -65.0867, -64.1, -64.0130, -36.0638, -11.7340,  7.2751, 9.2345, 8.2138, 8.4485, 16.6623, 17.3664, 27.2229, 65.0065, 65.0065, 21];
x_arc= [0, 110.7359, 222.1270, 269.3044, 278.7, 279, 324,327.625,330,354.875, 362.0867, 390, 390, 0, 0];
y_arc = [59.8435, 65.0065, 65.0065, 65.7106, 71.5, 81, 81,  77.875,68,55.875,50.4563, 50, 90, 90, 59.8435];
x_med1 = [-0.8003, 15.7387, 26.9188, 26.9188, 42, 37.6129, 27.4761, -0.2668, -0.8003];
y_med1 = [42.4867, 47, 47, 38.8094, 35.7846, 28.6037, 28.1250, 32.9122, 42.4867];
x_med2 = [354.3040, 354.5015, 370.6955, 390, 390, 361.7098, 354.3040];
y_med2 = [36.8040, 35.0490, 23.8727, 23.1338, 48.1650, 44.1009, 36.8040];
x_atl1 = [269.3044, 278.7, 292, 336, 367.0867, 361.7098, 354.3040, 354.5015, 390, 390, 326.3060, 290, 285.3592, 285.3592, 304.4254, 296.1760,297.7981,297.9455,298.0438,298.4370,298.8301,298.8301,299.0267,298.9776,298.7318,298.6335,298.0929,297.2575,293.8665,292.9327,292.0972,291.2618,290.0331,288.8536,288.0673,287.2318,286.2489,285.1677,282.5630,280,279.5799,277.6502, 253.8044, 269.3044];
y_atl1 = [65.7106, 71.5, 67, 70, 50.4563, 44.1009, 36.8040, 35.0490, -8.3367, -79.5093, -87.3012, -87.3012, -87.3012, -51.7887, -11.7340, 9.7150,10.7900,11.3500,11.7710,12.1450,12.7990,13.2660,13.9210,14.9490,15.7900,16.7240,17.5650,17.9390,18.6400,18.5000,18.5000,19.0140,19.3880,19.7150,19.8550,19.9020,20.1820,20.4160,20.5560,21.9110,25.1906,30.0180, 61.5698, 65.7106];
x_atl2 = [0, 40.3434, 40.3434, 0, 0];
y_atl2 = [12, -14.7888, -78.0464, -80, 12];
x_bs = [26.9188, 26.9188, 42, 42, 26.9188];
y_bs = [47, 38.8094, 35.7846, 47, 47];
x_rs = [27.4761, 40.9015, 43.5715, 44.3975, 37.6129, 27.4761];
y_rs = [28.1250, 10.3363, 12.9936, 17.2947, 28.6037, 28.1250];
x_cs = [283.7198, 280.8265, 279.2641, 277.4294, 269.3044, 262.2278, 256.4617, 264.9318, 277.6502, 279.5799, 280,282.5630,285.1677,286.2489,287.2318,288.0673,288.8536,290.0331,291.2618,292.0972,292.9327,293.8665,297.2575,298.0929,298.6335,298.7318,298.9776,299.0267,298.8301,298.8301,298.4370,298.0438,297.9455,297.7981,296.1760,293.6700,290.1810,287.1830,285.0690,283.7198];
y_cs = [7.2751, 9.2345, 8.2138, 8.4485, 16.6623, 17.3664, 27.2229, 32.1265, 30.0180, 25.1906, 21.9110,20.5560,20.4160,20.1820,19.9020,19.8550,19.7150,19.3880,19.0140,18.5000,18.5000,18.6400,17.9390,17.5650,16.7240,15.7900,14.9490,13.9210,13.2660,12.7990,12.1450,11.7710,11.3500,10.790,9.7150,9.5750,11.0230,11.0700,10.1360,7.2751];
% x_bb = [278.7,279,324,336,292,278.7];
% y_bb = [71.5,81,81,70,67,71.5];
% add JML 6/27/2023
x_bb = [278.7,279,324,321,292,278.7];
y_bb = [71.5,81,81,68.98,67,71.5];
x_cp=[44,44,56,56,44];
y_cp=[36,50,50,36,36];
% greatlakes

x_gl=[267.5,267.5,284.5,284.5,267.5];
y_gl=[41,50,50,41,41];
%Check for points in each basin in lat data
in_ind = inpolygon(lon,lat,x_ind,y_ind);
% figure;
% plot(x_ind,y_ind,'g');plot(lon(in_ind),lat(in_ind),'.');
% % hold on
in_pac = inpolygon(lon,lat,x_pac,y_pac);
% plot(x_pac,y_pac,'g'); plot(lon(in_pac),lat(in_pac),'.');
in_arc1 = inpolygon(lon,lat,x_arc,y_arc);
in_arc2 = inpolygon(lon_360,lat,x_arc,y_arc);
in_arc=in_arc1 | in_arc2;
% plot(x_arc,y_arc,'g'); plot(lon(in_arc),lat(in_arc),'.');
in_med1 = inpolygon(lon,lat,x_med1,y_med1);
in_med2 = inpolygon(lon,lat,x_med2,y_med2);

in_med=in_med1 |in_med2;
% plot(x_med1,y_med1,'g'); plot(lon(in_med),lat(in_med),'.');
% plot(x_med2,y_med2,'g')
in_atl1 = inpolygon(lon,lat,x_atl1,y_atl1);
in_atl1_360=inpolygon(lon_360,lat,x_atl1,y_atl1);
in_atl2 = inpolygon(lon,lat,x_atl2,y_atl2);
in_atl=in_atl1|in_atl2|in_atl1_360;
% plot(x_atl1,y_atl1,'g'); plot(lon(in_atl),lat(in_atl),'.')
% plot(x_atl2,y_atl2,'g')
in_bs = inpolygon(lon,lat,x_bs,y_bs);
% plot(x_bs,y_bs,'g'); plot(lon(in_bs),lat(in_bs),'.');
in_cs = inpolygon(lon,lat,x_cs,y_cs);
% plot(x_cs,y_cs,'g'); plot(lon(in_cs),lat(in_cs),'c');
in_bb = inpolygon(lon,lat,x_bb,y_bb);
% plot(x_bb,y_bb,'g'); plot(lon(in_bb),lat(in_bb),'k');
in_rs = inpolygon(lon,lat,x_rs,y_rs);
% plot(x_rs,y_rs,'g'); plot(lon(in_rs),lat(in_rs),'.');
in_cp=inpolygon(lon,lat,x_cp,y_cp);

in_gl=inpolygon(lon,lat,x_gl,y_gl);

in_atl(in_gl)=0;

pos_sulu=[8,117.4;11,119;11,121;14,121;14.1,123.5;11.8,124.8;10,125.8;8,125.8;8,123.5;8,122.6;7,122.2;6,120.6;5.5,118.5;5,118;8,117.4];
x_sulu=pos_sulu(:,2)';
y_sulu=pos_sulu(:,1)';
% plot(pos_sulu(:,2),pos_sulu(:,1))


pos_celebes=[5,118;-6,116;-6,120;.5,120;.8,120.3;.8,121;.7,123.6;1,124.6;1.2,125.15;3.3,125.8;8,125.8;8,123.5;8,122.6;7,122.2;6,120.6;5.5,118.5;5,118;];
x_celebes=pos_celebes(:,2)';
y_celebes=pos_celebes(:,1)';
%Save coordinate pairs found in each basin to separate arrays
in_sulu=inpolygon(lon,lat,x_sulu,y_sulu);
in_celebes=inpolygon(lon,lat,x_celebes,y_celebes);

in_pac(in_sulu|in_celebes)=0;

not_in_basin=~(in_ind|in_pac|in_arc|in_med|in_atl|in_bs|in_cs|in_bb|in_rs|in_cp);

% plot(lon(not_in_basin),lat(not_in_basin),'k.');


global_basins(1).name='indian ocean';
global_basins(1).pos=in_ind;
global_basins(2).name='pacific ocean';
global_basins(2).pos=in_pac;
global_basins(3).name='artic ocean';
global_basins(3).pos=in_arc;
global_basins(4).name='meditrain sea';
global_basins(4).pos=in_med;
global_basins(5).name='atlantic ocean';
global_basins(5).pos=in_atl;
global_basins(6).name='balck sea';
global_basins(6).pos=in_bs;
global_basins(7).name='red_sea';
global_basins(7).pos=in_rs;
global_basins(8).name='gulf of mexico and caribain sea';
global_basins(8).pos=in_cs;
global_basins(9).name='bafin bay';
global_basins(9).pos=in_bb;
global_basins(10).name='caspain sea';
global_basins(10).pos=in_cp;
global_basins(11).name='atric_360';
global_basins(11).pos=in_arc;
global_basins(12).name='sulu';
global_basins(12).pos=in_sulu;
global_basins(13).name='celebes';
global_basins(13).pos=in_celebes;


