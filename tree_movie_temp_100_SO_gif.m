
var_type='o2';


path_ERDDAP=TreeSetUp.path_ERDDAP;



tree_prefix=TreeSetUp.tree_prefix;


path_Figs=TreeSetUp.path_Figs;


end_year=TreeSetUp.end_year;

movie_name=['temp_100m_',num2str(floor(end_year)),'_SO_2015_new3.gif'];


subdir='yearly_withcycle';


path_nc_erddap=[path_ERDDAP,'netcdf\',tree_prefix,'\',subdir,'\'];

if var_type=='s'
     file_prefix='RFROMV22_SAL_';
elseif var_type=='t'
     file_prefix='RFROMV22_TEMP_';
elseif var_type=='h'
     file_prefix='RFROMV22_OHC_';
end
%%


start_year_movie=2020;
start_month_movie=1;
end_year_movie=floor(end_year);
end_year_movie=2020;
end_month_movie=12;

nyears_movie=end_year_movie-start_year_movie+1;

years_movie=repelem(start_year_movie:end_year_movie,12);
months_movie=repmat(1:12,[1 nyears_movie]);

months_movie=months_movie(start_month_movie:end-12+end_month_movie);
years_movie=years_movie(start_month_movie:end-12+end_month_movie);

nfiles=length(months_movie);





for ifile=1:nfiles

    %load in data
    year_load=years_movie(ifile);
    month_load=months_movie(ifile);

     if month_load>=10
          file_name_nc= [path_nc_erddap,file_prefix,num2str(year_load),'_',num2str(month_load),'.nc'];
       else
          file_name_nc= [path_nc_erddap,file_prefix,num2str(year_load),'_0',num2str(month_load),'.nc'];
     end

     lon_tpx=double(ncread(file_name_nc,'longitude'));
     lat_tpx=double(ncread(file_name_nc,'latitude'));

     time=double(ncread(file_name_nc,'time'));

     temp_total=ncread(file_name_nc,'ocean_temperature');

     ntime=length(time);





    for itime=1:ntime


        figure('Position', [100 100 800 900])
        clf
            set(gcf,'color','white');
            
            
            
            
            
            corrhc=double(squeeze((temp_total(:,:,13,itime))));
            date_name=datestr(time(itime)+datenum(1950,1,1));

           
            
            
            lon=lon_tpx';
            lat=lat_tpx';
        
 
             min_val=-1;
             max_val=20.;
             del_cont=1;
            
            val_lev=[min_val:1:max_val];


            
            ii=find(lon<30);
            jj=find(lon>=30);
            lon=[lon(jj),lon(ii)+360];
            corrhc=[corrhc(jj,:);corrhc(ii,:)];
            colormap(cmocean('thermal'));
            
            
            m_proj('stereographic','lat',-90,'lon',0,'radius',60)

             m_grid('xtick',12,'tickdir','in','linest','-','xaxislocation','top');

            [cs1,h1]=m_contourf(lon,lat,corrhc',[-1000,min_val:del_cont:max_val],'linestyle','none');

             m_grid('xtick',12,'tickdir','in','linest','-','xaxislocation','top');
 
            hold on

            a=gca;
            hold on
            m_coast('patch',[.8 .8 .8]);
            t2=m_text(-70,-80,date_name,'fontsize',15);
           
            caxis([min_val max_val-del_cont])
            
            hold off
           



            ja=axes('pos',[.262 .90-.0475 .51 .01]);
  
            colormap(cmocean('thermal'));
            
            [cs,h]=contourf([min_val:.01:max_val],[0 1],[1 1]'*[min_val:.01:max_val],[min_val:del_cont:max_val]);
            
            set(ja,'tickdir','out','xaxisl','top','xtick',val_lev,'ytick',[])
            caxis([min_val max_val-del_cont])
            japos=get(ja,'pos');
            set(ja,'XAxisLocation','bottom','pos',japos-[.075 .75 -.15 -.015])
            xlabel('(C^{\circ})','fontsize',12)
            

            if itime==1 & ifile==1
                exportgraphics(gcf, [path_Figs,movie_name]);
            else
                exportgraphics(gcf, [path_Figs,movie_name], Append=true);

            end
           
    end
