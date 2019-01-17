%Drag drop the Climate senarios and they will load 
%    Load CRU_DelT which is from Climate senarios 
%    Load CRU_RatP
% in the same floder there is file named data 
%   Load Data which is from excel
%Load the one hydrofile and change precipitaion namd into PRECIPE_R
%                           change precipitaion namd into TEMPS_R
%Chose Sen for senario run this code and save the file 

senar={'bccr_bcm2_0','cccma_cgcm3_1','cccma_cgcm3_1_t63','cnrm_cm3' ...
'csiro_mk3_0','csiro_mk3_5','gfdl_cm2_0','gfdl_cm2_1','giss_aom' ...
'giss_model_e_h','giss_model_e_r','iap_fgoals1_0_g','inmcm3_0','ipsl_cm4','miroc3_2_hires' ...
'miroc3_2_medres','mpi_echam5','mri_cgcm2_3_2a','ncar_ccsm3_0' ...
'ncar_pcm1','ukmo_hadcm3','ukmo_hadgem1'}

for(ind=1:22)
    %====================================
clearvars -except ind senar;
sennamm=senar{ind}
savepath ='C:/Works/Vietnam/Model4/Senarios/';
path='C:/Works/GCM Decadal Average/';

load([path,sennamm,'.mat'], '-regexp','CRU_DelT', 'CRU_RatP');
load([path,'Data.mat'],'-regexp','data');
load([savepath,'Vietname_Hydro_UNH.mat']);

PRECIPTS_R=PRECIPTS;
TEMPTS_R= TEMPTS;
clear('PRECIPTS','TEMPTS');



for(Sen=1:3) 
%------------------------------------------------------------------------
x=zeros(347,1200);

for(h = 1:347)
if (data(h,3)>0)     
    f = CRU_RatP(Sen,:,data(h,3),:)*data(h,4);
else
    f = CRU_RatP(1,:,1,:)*0;
end 
    g=f(:);
g=reshape(g,10,12);
x1=zeros(10,120);
for (j=1:10) 
for(i=1:120)
  bezinga=mod(i,12);
  if (bezinga==0)  
      bezinga=12; 
  end
 x1(j,i)=g(j,bezinga);       
end 
end 
x(h,:) =reshape(x1',1,120*10); 
end 

Delt=zeros(22,1200);
for(h=1:347)
Delt(data(h,2),:)=Delt(data(h,2),:) + x(h,:);
end 
clear('x','j','i','g','f','bezinga','h','x1')

PRECIPTS=Delt.*PRECIPTS_R(:,1:1200);

x=zeros(347,1200);

for(h = 1:347)
if (data(h,3)>0)     
    f = CRU_DelT(Sen,:,data(h,3),:)*data(h,4);
else
    f = CRU_DelT(1,:,1,:)*0;
end 
    g=f(:);
g=reshape(g,10,12);
x1=zeros(10,120);
for (j=1:10) 
for(i=1:120)
  bezinga=mod(i,12);
  if (bezinga==0)  
      bezinga=12; 
  end
 x1(j,i)=g(j,bezinga);       
end 
end 
x(h,:) =reshape(x1',1,120*10); 
end 

Delt=zeros(22,1200);
for(h=1:347)
Delt(data(h,2),:)=Delt(data(h,2),:) + x(h,:);
end 


TEMPTS=Delt+TEMPTS_R(:,1:1200);
TRANGETS=TRANGETS(:,1:1200);
years=1901:2000;

clear('x','j','i','g','f','bezinga','h','x1','Delt')

save([savepath,sennamm,'S',num2str(Sen,1),'.mat'],'LATm','PRECIPTS','TEMPTS','TRANGETS','areas','flowNames','obsRunoff','obsRunoffUNH','sub_basin','years')
%------------------------------------------------------------------------
end 
end 
