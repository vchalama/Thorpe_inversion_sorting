%inversion method
% Reading data files

close all
clear all
time_bulk=load('/media/Work2/filter/filter_all4/plane_data/time_bulk.txt');
th_bar=load('density_data_post.dat');
data_dir='//media/Work2/filter/filter_all4/datas_xy/'

ntt=300;
tc=170


    for i = 1:ntt
     time_t(i)=time_bulk(i,1) ;
    end 
nyy=193;
nzz=64;
grho=95.5;
ygrid(1:nyy,1:ntt)=0;
tgrid(1:nyy,1:ntt)=0;
thme(1:nyy,1:ntt)=0;
thme_s(1:nyy,1:nzz)=0;
disstotal(1:nzz,1:nyy)=0;
disstotalave(1:ntt,1:nyy)=0;
ygrid_l(1:nyy)=0;

i=1;
iii=1;
t=5;
   if ( t+tc <10 )  
        fdata = load(sprintf('%sspan3_000%1d.dat',data_dir,t+tc) );
    elseif ( 9 < t+tc  && t+tc< 100)    
        fdata = load(sprintf('%sspan3_00%2d.dat',data_dir,t+tc) );
    elseif ( 99 < t+tc && t+tc < 1000)    
        fdata = load(sprintf('%sspan3_0%3d.dat',data_dir,t+tc) );
    elseif ( 999 < t+tc && t+tc < 10000) 
        fdata = load(sprintf('%sspan3_%4d.dat',data_dir,t+tc) );
   end
 for k = 1:nzz
    for j = 1:nyy
     ygrid_l(j)=fdata(iii,2);
    iii=iii+1;
    end 
 end
 
 %%


% Initializing


for j = 1:nyy
    for kk = 1:ntt     
 	 tgrid(j,kk)=time_t(kk) ;
     ygrid(j,kk)=ygrid_l(j);
    end 
end
 ygrid18=ygrid(:,1);

%%


%%
for t=1:ntt % Reading each time step for the targeted LONG
    t+tc
    if ( t+tc <10 )  
        fdata = load(sprintf('%sspan3_000%1d.dat',data_dir,t+tc) );
    elseif ( 9 < t+tc  && t+tc< 100)    
        fdata = load(sprintf('%sspan3_00%2d.dat',data_dir,t+tc) );
    elseif ( 99 < t+tc && t+tc < 1000)    
        fdata = load(sprintf('%sspan3_0%3d.dat',data_dir,t+tc) );
    elseif ( 999 < t+tc && t+tc < 10000) 
        fdata = load(sprintf('%sspan3_%4d.dat',data_dir,t+tc) );
    end
    
    
 iii=1;    
 for k = 1:nzz
    for j = 1:nyy
     thme_s(j,k)=fdata(iii,3)+th_bar(j); % full instantaneous density 
    iii=iii+1;
    end 
 end
 
 
 for i=1:nzz  % Reading each column (LAT) at targeted LONG in each time step


      
den18=thme_s(:,i); % den18 = density of present column




sden18=den18;
I=[1:nyy]';
deviation(1:nyy)=0;
Lpsum(1:nyy)=0;
Np18(1:nyy)=0;
LTtmp(1:nyy)=0;
diss18(1:nyy)=0;
ygrid18s=ygrid18;
NT_OV=0;

jj=nyy-1;

% finding Inverted profiles
 while ( jj > 1 )

    
    if ( (den18(jj) - den18(jj-1)) > 0  )
          
             NY_TEMP_top=jj;
    while ( (den18(jj) - den18(jj-1)) > 0  )
             NY_TEMP_bot = jj-1;
             jj=jj-1;
    if ( jj == 1 )
        break
    end
    end
    
    NT_OV = NT_OV+1 ;
    NY_3D_TOT_TEMP=NY_TEMP_bot-NY_TEMP_top+1;

    
 % Dissipation estimation in each inverted profile
 
    if ( (NY_TEMP_top-NY_TEMP_bot) > 3 ) 
    clear den18temp
    clear sden18temp
    clear ygrid18temp ygrid18stemp
    clear A Aq
    
        ygrid18temp=ygrid18(NY_TEMP_bot:NY_TEMP_top);
        den18temp1=den18(NY_TEMP_bot:NY_TEMP_top);
    
    
                A(:,1)=ygrid18temp;
                maxg=max(ygrid18temp);
                ming=min(ygrid18temp);
                domaing=maxg-ming;
                leng=length(ygrid18temp);
                Aq(:,1)=ming:domaing/(leng-1):maxg;
                A(:,2)=den18temp1;
                Aq(:,2) = interp1(A(:,1),A(:,2),Aq(:,1));
                ygrid18temp=Aq(:,1);
                den18temp=Aq(:,2);
                [ sden18temp, I ]=sort(den18temp,'descend');


        for k = NY_TEMP_bot:NY_TEMP_top
          
        ygrid18stemp(k-NY_TEMP_bot+1)=ygrid18temp(I(k-NY_TEMP_bot+1));
        deviation(k)=ygrid18stemp(k-NY_TEMP_bot+1)-ygrid18temp(k-NY_TEMP_bot+1);
        dev2(k)=deviation(k)*deviation(k);
        end
        LT2=mean(dev2(NY_TEMP_bot:NY_TEMP_top));
        LT=sqrt(LT2);
        Np=-grho*(max(sden18(NY_TEMP_bot:NY_TEMP_top))-min(sden18(NY_TEMP_bot:NY_TEMP_top)));
        Np=sqrt(Np/(ygrid18(NY_TEMP_bot)-ygrid18(NY_TEMP_top)));
        disstmp=LT2*Np^3;
        diss18(NY_TEMP_bot:NY_TEMP_top)=disstmp  ;
        LTtmp(NY_TEMP_bot:NY_TEMP_top)=LT;
    end
    end
    jj=jj-1;
 end
 

LTt(i,:)=LTtmp(:);
deviationt(i,:)=deviation(:);
Lpsumt(i,:)=Lpsum(:);
disstotal(i,:)=diss18(:);

 end

 % Averaging for each LONG
 for j = 1:nyy
 disstotalave(t,j)=mean(disstotal(1:nzz,j));
 LTtave(t,j)=mean(LTt(1:nzz,j));
 end
 
end
%% plottong contours 

%     close all
%         logdisstotal=log10((disstotalave)+0.00000000001);
%        [C,h,CF] = contourf(tgrid',ygrid',logdisstotal,50);
%         set(h,'edgecolor','none');
%    %      caxis([-7,-3.5]);
%        title('Masoud','FontSize',18);
%        colorbar('location','southoutside','FontSize',10);
% %        xlim([239.0404 240.3748]);
%   %     ylim([-18.8 0]);
%        grid on
%% Saving for plotting contours
fp=fopen('DISS_THORP_MINE2.dat','w');
fprintf(fp,' title=" grid"\n ');
fprintf(fp,' VARIABLES = "t" "y" "eps" \n ');
fprintf(fp,' zone f=point, i=         %5i  j=         %5i  k=           1\n ',ntt,nyy);

for j=1:nyy
  for i=1:ntt
    fprintf(fp,'%19.10f %19.10f %19.10f  \n', tgrid(j,i),ygrid(j,i), disstotalave(i,j));
  end
end
fclose(fp);

%% Integration

for i=1:length(disstotalave(:,1))
eps_int2(i)=trapz(ygrid18,disstotalave(i,:));
end
%% Savinf the integration
fp=fopen('int_eps_cfd.dat','w');
fprintf(fp,' title=" grid"\n ');
fprintf(fp,' VARIABLES =  "time"  "int_eps" \n ');
fprintf(fp,' zone f=point, j=         %5i \n ',length(disstotalave(:,1)));
for k=1:length(disstotalave(:,1))
    fprintf(fp,'%18.10e\t%16.8e \n',tgrid(1,k),abs(eps_int2(k)));
end
fclose(fp);
