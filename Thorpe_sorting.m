%Thorpe method

% Reading data files
close all
clear all
time_bulk=load('/media/Work2/plane_data/time_bulk.txt');
th_bar=load('density_data_post.dat');
data_dir='//media/Work2/datas_xy/'

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
U_s_ave(1:ntt,1:nyy)=0;

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
 
A(:,1)=ygrid18;
maxg=max(ygrid18);
ming=min(ygrid18);
domaing=maxg-ming;
leng=length(ygrid18);
Aq(:,1)=ming:domaing/(leng-1):maxg;



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
     U_s(j,k)=fdata(iii,4); % full instantaneous velocity
    iii=iii+1;
    end 
 end
 
    
 for i=1:nzz % Reading each column (LAT) at targeted LONG in each time step


     
den18=thme_s(:,i);  %den18 = density of present column

% plotting from stretched to uniform grid
A(:,2)=den18;
Aq(:,2) = interp1(A(:,1),A(:,2),Aq(:,1));
ygrid18=Aq(:,1);
den18=Aq(:,2); 


sden18=den18;
I=[1:nyy]';
deviation(1:nyy)=0;
Lpsum(1:nyy)=0;
Np18(1:nyy)=0;
LTtmp(1:nyy)=0;
diss18(1:nyy)=0;
ygrid18s=ygrid18;


[ sden18, I ]=sort(den18,'descend'); %sden18 = sorted density of the present column

for k = 1:nyy
  ygrid18s2(k)=ygrid18(I(k)); %ygrid18s2 = corresponding height for sorted density of the present column
end

% calculating LT and disssipation
deviation(nyy)=ygrid18(nyy)-ygrid18s2(nyy);
  dev2(nyy)=deviation(nyy)*deviation(nyy);
  Lpsum(nyy)=deviation(nyy);

for k = nyy-1:-1:1
  deviation(k)=ygrid18(k)-ygrid18s2(k);
  dev2(k)=deviation(k)*deviation(k);
  Lpsum(k)=Lpsum(k+1)+deviation(k);
end


Lzeros2=find(abs(Lpsum) < 10^-12);
Lzeros(2:length(Lzeros2)+1)=Lzeros2;
Lzeros(1)=1;


for k=2:length(Lzeros)
    if (  abs(Lzeros(k-1)-Lzeros(k)) > 3 )
        kup=Lzeros(k-1);
        kbottom=Lzeros(k);
        LT2=mean(dev2(kup:kbottom));
        LT=sqrt(LT2);
        NP=((max(sden18(kup:kbottom))-min(sden18(kup:kbottom)))/(max(ygrid18(kup:kbottom))-min(ygrid18(kup:kbottom))));
        NP=sqrt(abs(NP*grho));
        disstmp=0.64*LT2*NP^3;
        diss18(kup:kbottom)=disstmp  ;
        LTtmp(kup:kbottom)=LT;
        
     end
end

LTt(i,:)=LTtmp(:);
deviationt(i,:)=deviation(:);
Lpsumt(i,:)=Lpsum(:);
disstotal(i,:)=diss18(:);

%averaging
 for j = 1:nyy
 disstotalave(t,j)=mean(disstotal(1:nzz,j));
 U_s_ave(t,j)=mean(U_s(j,1:nzz));
 end
 
 end

 
end




%% ploting for contours
fp=fopen('DISS_THORP.dat','w');
fprintf(fp,' title=" grid"\n ');
fprintf(fp,' VARIABLES = "t" "y" "eps" "U" \n ');
fprintf(fp,' zone f=point, i=         %5i  j=         %5i  k=           1\n ',ntt,nyy);

for j=1:nyy
  for i=1:ntt
    fprintf(fp,'%19.10f %19.10f %19.10f  \n', tgrid(j,i),ygrid18(j), disstotalave(i,j), U_s_ave(i,j));
  end
end
fclose(fp);

%% integrating
for i=1:length(disstotalave(:,1))
eps_int2(i)=trapz(ygrid18,disstotalave(i,:));
end

fp=fopen('int_eps_ocean.dat','w');
fprintf(fp,' title=" grid"\n ');
fprintf(fp,' VARIABLES =  "time"  "int_eps" \n ');
fprintf(fp,' zone f=point, j=         %5i \n ',length(disstotalave(:,1)));
for k=1:length(disstotalave(:,1))
    fprintf(fp,'%18.10e\t%16.8e \n',tgrid(1,k),abs(eps_int2(k)));
end
fclose(fp);
