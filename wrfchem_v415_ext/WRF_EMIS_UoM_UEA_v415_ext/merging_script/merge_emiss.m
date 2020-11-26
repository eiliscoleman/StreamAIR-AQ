% This script will take the NAEI emission maps and replace all zero emissions
% with emissions from TNO maps. We are doing this for three domains.

for nn=1:63
%clear
%%
 if nn==1
 fname = 'wrfchemi_d01_2010-07-10_00:00:00'
 elseif nn==2
 fname = 'wrfchemi_d01_2010-07-11_00:00:00'
 elseif nn==3
 fname = 'wrfchemi_d01_2010-07-12_00:00:00'
 elseif nn==4
 fname = 'wrfchemi_d01_2010-07-13_00:00:00'
 elseif nn==5
 fname = 'wrfchemi_d01_2010-07-14_00:00:00'
 elseif nn==6
 fname = 'wrfchemi_d01_2010-07-15_00:00:00'
 elseif nn==7
 fname = 'wrfchemi_d01_2010-07-16_00:00:00'
 elseif nn==8
 fname = 'wrfchemi_d01_2010-07-17_00:00:00'
 elseif nn==9
 fname = 'wrfchemi_d01_2010-07-18_00:00:00'
 elseif nn==10
 fname = 'wrfchemi_d01_2010-07-19_00:00:00'
 elseif nn==11
 fname = 'wrfchemi_d01_2010-07-20_00:00:00'
 elseif nn==12
 fname = 'wrfchemi_d01_2010-07-21_00:00:00'
 elseif nn==13
 fname = 'wrfchemi_d01_2010-07-22_00:00:00'
 elseif nn==14
 fname = 'wrfchemi_d01_2010-07-23_00:00:00'
 elseif nn==15
 fname = 'wrfchemi_d01_2010-07-24_00:00:00'
 elseif nn==16
 fname = 'wrfchemi_d01_2010-07-25_00:00:00'
 elseif nn==17
 fname = 'wrfchemi_d01_2010-07-26_00:00:00'
 elseif nn==18
 fname = 'wrfchemi_d01_2010-07-27_00:00:00'
 elseif nn==19
 fname = 'wrfchemi_d01_2010-07-28_00:00:00'
 elseif nn==20
 fname = 'wrfchemi_d01_2010-07-29_00:00:00'
 elseif nn==21
 fname = 'wrfchemi_d01_2010-07-30_00:00:00'
 
 elseif nn==22
 fname = 'wrfchemi_d02_2010-07-10_00:00:00'
 elseif nn==23
 fname = 'wrfchemi_d02_2010-07-11_00:00:00'
 elseif nn==24
 fname = 'wrfchemi_d02_2010-07-12_00:00:00'
 elseif nn==25
 fname = 'wrfchemi_d02_2010-07-13_00:00:00'
 elseif nn==26
 fname = 'wrfchemi_d02_2010-07-14_00:00:00'
 elseif nn==27
 fname = 'wrfchemi_d02_2010-07-15_00:00:00'
 elseif nn==28
 fname = 'wrfchemi_d02_2010-07-16_00:00:00'
 elseif nn==29
 fname = 'wrfchemi_d02_2010-07-17_00:00:00'
 elseif nn==30
 fname = 'wrfchemi_d02_2010-07-18_00:00:00'
 elseif nn==31
 fname = 'wrfchemi_d02_2010-07-19_00:00:00'
 elseif nn==32
 fname = 'wrfchemi_d02_2010-07-20_00:00:00'
 elseif nn==33
 fname = 'wrfchemi_d02_2010-07-21_00:00:00'
 elseif nn==34
 fname = 'wrfchemi_d02_2010-07-22_00:00:00'
 elseif nn==35
 fname = 'wrfchemi_d02_2010-07-23_00:00:00'
 elseif nn==36
 fname = 'wrfchemi_d02_2010-07-24_00:00:00'
 elseif nn==37
 fname = 'wrfchemi_d02_2010-07-25_00:00:00'
 elseif nn==38
 fname = 'wrfchemi_d02_2010-07-26_00:00:00'
 elseif nn==39
 fname = 'wrfchemi_d02_2010-07-27_00:00:00'
 elseif nn==40
 fname = 'wrfchemi_d02_2010-07-28_00:00:00'
 elseif nn==41
 fname = 'wrfchemi_d02_2010-07-29_00:00:00'
 elseif nn==42
 fname = 'wrfchemi_d02_2010-07-30_00:00:00' 

 elseif nn==43
 fname = 'wrfchemi_d03_2010-07-10_00:00:00'
 elseif nn==44
 fname = 'wrfchemi_d03_2010-07-11_00:00:00'
 elseif nn==45
 fname = 'wrfchemi_d03_2010-07-12_00:00:00'
 elseif nn==46
 fname = 'wrfchemi_d03_2010-07-13_00:00:00'
 elseif nn==47
 fname = 'wrfchemi_d03_2010-07-14_00:00:00'
 elseif nn==48
 fname = 'wrfchemi_d03_2010-07-15_00:00:00'
 elseif nn==49
 fname = 'wrfchemi_d03_2010-07-16_00:00:00'
 elseif nn==50
 fname = 'wrfchemi_d03_2010-07-17_00:00:00'
 elseif nn==51
 fname = 'wrfchemi_d03_2010-07-18_00:00:00'
 elseif nn==52
 fname = 'wrfchemi_d03_2010-07-19_00:00:00'
 elseif nn==53
 fname = 'wrfchemi_d03_2010-07-20_00:00:00'
 elseif nn==54
 fname = 'wrfchemi_d03_2010-07-21_00:00:00'
 elseif nn==55
 fname = 'wrfchemi_d03_2010-07-22_00:00:00'
 elseif nn==56
 fname = 'wrfchemi_d03_2010-07-23_00:00:00'
 elseif nn==57
 fname = 'wrfchemi_d03_2010-07-24_00:00:00'
 elseif nn==58
 fname = 'wrfchemi_d03_2010-07-25_00:00:00'
 elseif nn==59
 fname = 'wrfchemi_d03_2010-07-26_00:00:00'
 elseif nn==60
 fname = 'wrfchemi_d03_2010-07-27_00:00:00'
 elseif nn==61
 fname = 'wrfchemi_d03_2010-07-28_00:00:00'
 elseif nn==62
 fname = 'wrfchemi_d03_2010-07-29_00:00:00'
 elseif nn==63
 fname = 'wrfchemi_d03_2010-07-30_00:00:00' 
end
%%
path_naei= '/Users/steve/emissions/3_domains/naei/summer/half_nox/';
path_tno= '/Users/steve/emissions/3_domains/tno/summer/half_nox/';


nc_naei=netcdf(strcat(path_naei,fname),'write')
nc_tno=netcdf(strcat(path_tno,fname),'nowrite')



total_naei = nc_naei{2}(:,:,:,:);

for n=3:42
   total_naei = total_naei + nc_naei{n}(:,:,:,:);
end

% d=size(nc_naei{2}) 
% d=size(total_naei)
% d=numel(total_naei)
% d=size(nc_tno{2})
%%

naei_zeros = find(total_naei==0);

%%

for n=2:42
   nc_test = nc_naei{n}(:,:,:,:);
   nc_in = nc_tno{n}(:,:,:,:);
   nc_test(naei_zeros) = nc_in(naei_zeros);
   nc_naei{n}(:,:,:,:) = nc_test;
end

%%
nc_naei{'E_OC_DOM'}(:,:,:,:) = nc_tno{'E_OC_DOM'}(:,:,:,:); % oc_dom
nc_naei{'E_OC_TRA'}(:,:,:,:) = nc_tno{'E_OC_TRA'}(:,:,:,:); % oc_tra
nc_naei{'E_BC_1'}(:,:,:,:) = nc_tno{'E_BC_1'}(:,:,:,:); % bc_1
nc_naei{'E_EC_1_25'}(:,:,:,:) = nc_tno{'E_EC_1_25'}(:,:,:,:); % EC_1_25
nc_naei{'E_OC_25_10'}(:,:,:,:) = nc_tno{'E_OC_25_10'}(:,:,:,:); % OC_25_10
nc_naei{'E_EC_25_10'}(:,:,:,:) = nc_tno{'E_EC_25_10'}(:,:,:,:); % EC_25_10
nc_oin25= nc_naei{'E_PM25'}(:,:,:,:) - nc_tno{'E_BC_1'}(:,:,:,:)- nc_tno{'E_EC_1_25'}(:,:,:,:)...
    - nc_tno{'E_OC_DOM'}(:,:,:,:)- nc_tno{'E_OC_TRA'}(:,:,:,:);

nc_oin10= nc_naei{'E_PM_10'}(:,:,:,:) - nc_tno{'E_OC_25_10'}(:,:,:,:)- nc_tno{'E_EC_25_10'}(:,:,:,:);

pm25=nc_naei{'E_PM25'}(:,:,:,:); % pm25
%pm10=nc_naei{'E_PM_10'}(:,:,:,:); % pm10

nc_oin25= min(nc_oin25,pm25*0.1) ;



nc_oin25(nc_oin25<=0.0) = 0.0;
nc_oin10(nc_oin10<=0.0) = 0.0;


clear
end
%%

