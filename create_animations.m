%% plot temp in NWA
plot_year = 2000:2020;
plot_depth = 150;
param_to_plot = 'temp';
reg = 'other';
quant = 'abs';
lat_lims = [20 45];
lon_lims = [280 320];
val_lims = [0 26];
step_val = 1;
plot_animation_lr_vs_hr(plot_year,plot_depth,param_to_plot,...
    reg,quant,lat_lims,lon_lims,val_lims,step_val,'off')

%% plot sal in NWA
plot_year = 2000:2020;
plot_depth = 150;
param_to_plot = 'sal';
reg = 'other';
quant = 'abs';
lat_lims = [20 45];
lon_lims = [280 320];
val_lims = [34.7 37.5];
step_val = 0.1;
plot_animation_lr_vs_hr(plot_year,plot_depth,param_to_plot,...
    reg,quant,lat_lims,lon_lims,val_lims,step_val,'off','left','no')

%% plot O2 in NWA
plot_year = 2000:2020;
plot_depth = 150;
param_to_plot = 'o2';
reg = 'other';
quant = 'abs';
lat_lims = [20 45];
lon_lims = [280 320];
val_lims = [140 260];
step_val = 5;
plot_animation_lr_vs_hr(plot_year,plot_depth,param_to_plot,...
    reg,quant,lat_lims,lon_lims,val_lims,step_val,'on')

%% plot NO3 in NWA
plot_year = 2000:2020;
plot_depth = 150;
param_to_plot = 'no3';
reg = 'other';
quant = 'abs';
lat_lims = [20 45];
lon_lims = [280 320];
val_lims = [0 20];
step_val = 1;
plot_animation_lr_vs_hr(plot_year,plot_depth,param_to_plot,...
    reg,quant,lat_lims,lon_lims,val_lims,step_val,'off')

%% plot temp vs. o2 in NWA
plot_year = 2004:2020;
plot_depth = 150;
param_to_plot_phys = 'temp';
param_to_plot_bgc = 'o2';
type = 'hr';
reg = 'other';
quant = 'abs';
lat_lims = [20 45];
lon_lims = [280 320];
val_lims_1 = [0 26];
step_val_1 = 1;
val_lims_2 = [100 250];
step_val_2 = 5;
plot_animation_phys_vs_bgc(plot_year,plot_depth,param_to_plot_phys,...
    param_to_plot_bgc,type,reg,quant,lat_lims,lon_lims,val_lims_1,step_val_1,...
    val_lims_2,step_val_2,'off')

%% plot temp vs. dic in NWA
plot_year = 2004:2020;
plot_depth = 150;
param_to_plot_phys = 'temp';
param_to_plot_bgc = 'dic';
type = 'hr';
reg = 'other';
quant = 'abs';
lat_lims = [20 45];
lon_lims = [280 320];
val_lims_1 = [0 26];
step_val_1 = 1;
val_lims_2 = [2050 2200];
step_val_2 = 5;
plot_animation_phys_vs_bgc(plot_year,plot_depth,param_to_plot_phys,...
    param_to_plot_bgc,type,reg,quant,lat_lims,lon_lims,val_lims_1,step_val_1,...
    val_lims_2,step_val_2,'off','left','no')

%% plot sal vs. no3 in NWA
plot_year = 2004:2020;
plot_depth = 150;
param_to_plot_phys = 'sal';
param_to_plot_bgc = 'no3';
type = 'hr';
reg = 'other';
quant = 'abs';
lat_lims = [20 45];
lon_lims = [280 320];
val_lims_1 = [34.7 37.5];
step_val_1 = 0.1;
val_lims_2 = [0 20];
step_val_2 = 1;
plot_animation_phys_vs_bgc(plot_year,plot_depth,param_to_plot_phys,...
    param_to_plot_bgc,type,reg,quant,lat_lims,lon_lims,val_lims_1,step_val_1,...
    val_lims_2,step_val_2,'off')

%% plot temp in GoA
plot_year = 2000:2020;
plot_depth = 150;
param_to_plot = 'temp';
reg = 'other';
quant = 'abs';
lat_lims = [35 60];
lon_lims = [195 240];
val_lims = [2 15];
step_val = 0.5;
plot_animation_lr_vs_hr(plot_year,plot_depth,param_to_plot,...
    reg,quant,lat_lims,lon_lims,val_lims,step_val,'off','right')

%% plot sal in GoA
plot_year = 2000:2020;
plot_depth = 150;
param_to_plot = 'sal';
reg = 'other';
quant = 'abs';
lat_lims = [35 60];
lon_lims = [195 240];
val_lims = [32.5 34.5];
step_val = 0.1;
plot_animation_lr_vs_hr(plot_year,plot_depth,param_to_plot,...
    reg,quant,lat_lims,lon_lims,val_lims,step_val,'on','right')

%% plot o2 in GoA
plot_year = 2004:2020;
plot_depth = 150;
param_to_plot = 'o2';
reg = 'other';
quant = 'abs';
lat_lims = [35 60];
lon_lims = [195 240];
val_lims = [50 300];
step_val = 10;
plot_animation_lr_vs_hr(plot_year,plot_depth,param_to_plot,...
    reg,quant,lat_lims,lon_lims,val_lims,step_val,'off','right','no')

%% plot no3 in GoA
plot_year = 2004:2020;
plot_depth = 150;
param_to_plot = 'no3';
reg = 'other';
quant = 'abs';
lat_lims = [35 60];
lon_lims = [195 240];
val_lims = [0 40];
step_val = 2;
plot_animation_lr_vs_hr(plot_year,plot_depth,param_to_plot,...
    reg,quant,lat_lims,lon_lims,val_lims,step_val,'off','right','no')

%% plot temp vs. o2 in GoA
plot_year = 2004:2020;
plot_depth = 150;
param_to_plot_phys = 'temp';
param_to_plot_bgc = 'o2';
type = 'hr';
reg = 'other';
quant = 'abs';
lat_lims = [35 60];
lon_lims = [195 240];
val_lims_1 = [2 15];
step_val_1 = 1;
val_lims_2 = [50 300];
step_val_2 = 5;
plot_animation_phys_vs_bgc(plot_year,plot_depth,param_to_plot_phys,...
    param_to_plot_bgc,type,reg,quant,lat_lims,lon_lims,val_lims_1,step_val_1,...
    val_lims_2,step_val_2,'off')

%% plot temp vs. o2 in CCS
plot_year = 2016:2025;
plot_depth = 150;
param_to_plot_phys = 'temp';
param_to_plot_bgc = 'o2';
type = 'hr';
reg = 'other';
quant = 'abs';
lat_lims = [10 50];
lon_lims = [205 250];
val_lims_1 = [2 25];
step_val_1 = 1;
val_lims_2 = [0 280];
step_val_2 = 10;
plot_animation_phys_vs_bgc(plot_year,plot_depth,param_to_plot_phys,...
    param_to_plot_bgc,type,reg,quant,lat_lims,lon_lims,val_lims_1,step_val_1,...
    val_lims_2,step_val_2,'off','left','no')

%% plot sal vs. no3 in CCS
plot_year = 2004:2020;
plot_depth = 150;
param_to_plot_phys = 'sal';
param_to_plot_bgc = 'no3';
type = 'hr';
reg = 'other';
quant = 'abs';
lat_lims = [10 50];
lon_lims = [205 250];
val_lims_1 = [32.5 35.5];
step_val_1 = 0.15;
val_lims_2 = [0 40];
step_val_2 = 2;
plot_animation_phys_vs_bgc(plot_year,plot_depth,param_to_plot_phys,...
    param_to_plot_bgc,type,reg,quant,lat_lims,lon_lims,val_lims_1,step_val_1,...
    val_lims_2,step_val_2,'off','left','no')

%% plot temp in GoM
plot_year = 2000:2020;
plot_depth = 150;
param_to_plot = 'temp';
reg = 'other';
quant = 'abs';
lat_lims = [15 33];
lon_lims = [262 290];
val_lims = [5 25];
step_val = 1;
plot_animation_lr_vs_hr(plot_year,plot_depth,param_to_plot,...
    reg,quant,lat_lims,lon_lims,val_lims,step_val,'off','right','no')

%% plot dic in GoM
plot_year = 2000:2020;
plot_depth = 150;
param_to_plot = 'dic';
reg = 'other';
quant = 'abs';
lat_lims = [15 33];
lon_lims = [262 290];
val_lims = [2000 2250];
step_val = 10;
plot_animation_lr_vs_hr(plot_year,plot_depth,param_to_plot,...
    reg,quant,lat_lims,lon_lims,val_lims,step_val,'off','left','no');
