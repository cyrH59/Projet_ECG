%% Partie 5 : 
clear all

[file,path] = uigetfile('*.mat', 'rt');
signal = load(fullfile(path, file));
data = -signal.ecg;
data2=-signal.ecg;
%for ecgnormal4.m ecgVF.m it must take off the - : 
%data = signal.ecg; data2 = signal.ecg;

Fs = signal.Fs;
N = size(data,2);
time_axis = (1:N)/Fs;
%% test : 
close all

%low pass filter
b1=zeros(1,13);
b1(1)=1;
b1(7)=-2;
b1(13)=1;
a1=[1 -2 1];

%high pass filter
b2=zeros(1,33);
a2=[1 -1];
b2(1)=-1;
b2(17)=32;
b2(18)=-32;
b2(33)=1;

%differenciator filter
b3=[1 2 0 -2 -1];
a3= [8/Fs];

% representation of the signal form
t = linspace(0,length(data)/Fs,length(data));
figure()
plot(t,data);
title("ecg")
xlabel("Time")
ylabel("Amplitude")

% low pass filter
passbas=filter(b1,a1,data);
passbas=passbas./max(passbas);
figure()
plot(t(1:3*Fs),passbas(1:3*Fs));
title("low pass");
xlabel("Frequency(Hz)")

%pass band filter by association with a low pass filter and high pass
%filter
passhaut=filter(b2,a2,data);
passhaut=passhaut./max(passhaut);
figure()
plot(t(1:Fs),passhaut(1:Fs));
title("band pass");
xlabel("signal filtered by a band pass filter")

%derivation using the differentiator filter
taille=length(passhaut);
retard=[0 0 passhaut]; %introduction of a  delay because of previous filters
passederive=filter(b3,a3, retard);
passederive=passederive./max(passederive);

figure()
plot(t(1:1*Fs),passederive(1:1*Fs));
title("signal filtered by a differentiating filter");
xlabel("Frequency(Hz)");
ylabel("Magnitude")

%square signal of signal filtered
ssq=abs(passederive).^2;
plot(t(1:1*Fs),ssq(1:1*Fs));
title("square signal");
xlabel("time(s)");
ylabel("Amplitude");

%windowing of signal filtered
N=20; 
smwi=zeros(1,length(data));
for i=101:length(smwi)
   for j=1:N
       smwi(i)=smwi(i)+(1/N)*ssq(i-j);
   end
end
   
figure()   
plot(t(1:3*Fs),smwi(1:3*Fs));
title("signal after windowing");
xlabel("Time(s)");
ylabel("Magnitude");


%QRS detection: 

seuil = 0.005;
different_seuil=zeros(2,1000);  %allow to detect different values of abscisse when the threehold is reached
different_indice=zeros(2,1000); %allow to detect different values of index of abscisse when the threehol is reached


p=0
for i=2:length(data)
    if smwi(i-1)<seuil && smwi(i)>seuil
        if j<(length(data)-10)
        j=i+1;
        end
        while smwi(j)>seuil && j<(length(data)-10)
            j=j+1;
        end
        p=p+1;
        different_seuil(1,p)=t(i+1); %it's axis value where threshold is reached for the first time
        different_seuil(2,p)=t(j); %it's axis value where threshold isn't reached starting from different_seuil(1,p)
        different_indice(1,p)=i+1; %it's axis index value where threshold is reached for the first time
        different_indice(2,p)=j; %it's axis index value where threshold isn't reached starting from different_seuil(1,p)
        
        
    end 
end

figure() 
periode=100; %number of periods

% for the file with vitricular fibrilation, it must add period because
% fibrilation is only on the second part of the signal. it can take 500
% periods. The load of charge is quite long. For the 4 normals signals, it
% took 100 periods
plot(t(1:periode*Fs),smwi(1:periode*Fs));
hold on;
for k=1:periode
scatter(different_seuil(1,k),0, 'filled','r');
scatter(different_seuil(2,k),0, 'filled','r');
end

title("signal after windowing with limits");
xlabel("Time(s)")
ylabel("Amplitude")

figure()
plot(t(1:periode*Fs),data(1:periode*Fs));
hold on;
for k=1:periode
scatter(different_seuil(1,k)-0.15,0, 'filled','r');
scatter(different_seuil(2,k)-0.115,0, 'filled','r');

end
title("ecg before being filtrated during some periods with marks that allow to detect QRS zone")
xlabel("Time(s)")
ylabel("Amplitude")


retar=0.12; % introduce of a delay linked to signal filtered. This delay is calculate with fvtool
indicedecallage=retar/0.005;
seuilecg=different_seuil-retar;
indiceecg=different_indice-indicedecallage;

figure()
plot(t(1:periode*Fs),data(1:periode*Fs));
hold on;

% initialisation of differents tabs with index values and axis values
%for R,Q,S,T,P
R=zeros(1,periode);
Rabs=zeros(1,periode);
Q=zeros(1,periode);
Qabs=zeros(1,periode);
S=zeros(1,periode);
Sabs=zeros(1,periode);
T=zeros(1,periode);
Tabs=zeros(1,periode);
Pabs=zeros(1,periode);
P=zeros(1,periode);
 
for v=1:periode
    minb=9999;
    maxb=-9999;
    
    indicemaxb=1;
    indiceminb=1;
    indiceminapresmax=1;
    minapresmaxb=9999;
    %it go on the differents intervalls et it calculate the maximum for each. 
    for m=indiceecg(1,v):indiceecg(2,v)
        if data(m)>maxb
            maxb=data(m);
            indicemaxb=m;
            
        end
    end
   % it estimate the first minimum before the previous maximum for each
   % intervall
    for z=indiceecg(1,v):indicemaxb
        if data(z)<minb
            minb=data(z);
            indiceminb=z;
        end
    end
    % it estimate the first minimum after the previous maximum for each
   % intervall
    for u=indicemaxb:indiceecg(2,v)
        if data(u)<minapresmaxb
            minapresmaxb=data(u);
            indiceminapresmax=u;
        end
        
        
    end
    % it put the different key points
     scatter(t(indicemaxb),maxb,'r','filled');
     text(t(indicemaxb),maxb,' R ','Color','red','FontSize',14);
     scatter(t(indiceminb),minb,'r','filled');
     text(t(indiceminb),minb,' Q ','Color','red','FontSize',14);
     scatter(t(indiceminapresmax),minapresmaxb,'r','filled');
     text(t(indiceminapresmax),minapresmaxb,' S ','Color','red','FontSize',14);
     
     %it put the value of index and axis values on our tabs 
     R(1,v)=maxb;
     Rabs(1,v)=indicemaxb;
     Q(1,v)=minb;
     Qabs(1,v)=indiceminb;
     S(1,v)= minapresmaxb;
     Sabs(1,v)=indiceminapresmax;
     
    
        

    


     end
 

title("ecg avant filtrage sur 1 secondes")
xlabel("Temps(s)")
ylabel("Magnitude")

%filter in order to estimate P and T point. 

b4=[1 0 0 0 0 0 -1];
a4= [1];

b5=[1 0 0 0 0 0 0 0 -1];
a5= [1 -1] ;

 

filtre1=filter(b4,a4,data2);
figure()
plot(t(1:3*Fs),data2(1:3*Fs));
title("differentiating filter");


filtre2=filter(b5,a5,filtre1);
figure()
plot(t(1:3*Fs),filtre2(1:3*Fs));
hold on;
plot(t(1:3*Fs),data(1:3*Fs),'r');
title("low pass filter for T and P and comparaison with the original signal ");
xlabel("Frequency(Hz)")
ylabel("Amplitude")
legend('original signal','signal filtered')

%Tableau permet de placer les valeurs des abscisses où le filtre s'annule
% Tab that allow to determine easily where the filter is very close to zero
Abscissezero=zeros(1,length(filtre2));


% loop that allow to detect the change of sign and the passage by zero.
for f=2:length(filtre2)
    if (filtre2(f-1)*filtre2(f)<=0)
        if abs(filtre2(f-1))>abs(filtre2(f))
            Abscissezero(1,f)=1;
        
        else
             Abscissezero(1,f-1)=1;
        end
        
   
            
        
    end
    
end

% same proccess but with index, it put in a tab the values of index where
% the filtered signal is null. 
indicezero=zeros(1,1000);
compteur=0;
for ui=1:length(data)

if Abscissezero(1,ui)==1
    compteur=compteur+1;
    indicezero(1,compteur)=ui;
end 


end

figure()
plot(t(1:periode*Fs),data(1:periode*Fs));
hold on;
retard=-4;  % delay due to the previous filtering and determined with fvtool 
threshold = -20./max(data);

for p0=1 : periode-1
    compteur = 0 ;
    compteur2= 0;
    maxy=-9999;
    absmaxy=0
    % it take all previous index 
    % it made a "if" in order to check different condition 
    % it must have the point under the threshold and under the next R and
    % over the previous S
    
    for  lk=1:length(indicezero)
        if data(indicezero(lk)+retard)> threshold && indicezero(lk)+retard < (Rabs(p0+1)-floor(0.3*(Rabs(p0+1)-Rabs(p0)))) && indicezero(lk)+retard > Sabs(p0)
            compteur = compteur +1 ;
            scatter(t(indicezero(lk)+retard),data(indicezero(lk)+retard),'r','filled');
            if compteur == 1
                % it put our point on our tab
                Tabs(p0) = indicezero(lk)+retard ;
                T(p0) = data(indicezero(lk)+retard);
            end
        
           
             
        end
        end
      
end
    
Rlength=length(R);
Tlength=length(T);
for pq=1:Rlength-1
maxy=-9999;
absmaxy=0;
% For P, it didn't find how use correctly the results of our filter. it use
% an alternative solve. it determine the max beetween 0.7RR ans R intervale
% and it check the position of this point with the others points. This
% method give good results
for hj=Rabs(pq):Rabs(pq+1)
    if data(retard+hj)>maxy && (retard+hj)>Tabs(pq) && (retard+hj)<Qabs(pq+1) && (retard+hj)>(Rabs(pq)+floor(0.7*(Rabs(pq+1)-Rabs(pq))))
        maxy=data(retard+hj);
        absmaxy=retard+hj;
    end
end
%it put our point in our tab
Pabs(pq) =absmaxy;
P(pq) = maxy;
end


%final chart
figure()
plot(t(1:Fs*periode),data(1:Fs*periode)./max(data));
hold on;
xlabel("Time(s)");
ylabel("Magnitude ");
title("Spotting of different keys points on ecgnormal4.m");

%it put all our points
for ma=1:periode
     
     scatter(t(Rabs(ma)),R(ma)./max(data),'r','filled');
     text(t(Rabs(ma)),R(ma)./max(data),' R ','Color','red','FontSize',14);
     scatter(t(Qabs(ma)),Q(ma)./max(data),'r','filled');
     text(t(Qabs(ma)),Q(ma)./max(data),' Q ','Color','red','FontSize',14);
     scatter(t(Sabs(ma)),S(ma)./max(data),'r','filled');
     text(t(Sabs(ma)),S(ma)./max(data),' S ','Color','red','FontSize',14);
     if(Pabs(ma) ~= 0)
     scatter(t(Pabs(ma)+1),P(ma)./max(data),'r','filled');
     text(t(Pabs(ma)+1),P(ma)./max(data),' P ','Color','red','FontSize',14);
     end
     if(Tabs(ma) ~= 0)
     scatter(t(Tabs(ma)),T(ma)./max(data),'r','filled');
     text(t(Tabs(ma)),T(ma)./max(data),' T ','Color','red','FontSize',14);
     end

end


%% Bradycardie Tachycardie 

BPM=zeros(1,length(Rabs));

% it go on BPM tab for each points and it calculate our BPM for each itération.

for ia=1:length(BPM)-1
    ecart=Rabs(ia+1)-Rabs(ia);
    ecartseconde=ecart/Fs;
    BPM(ia)=60/ecartseconde;
end
BPM(length(BPM))=BPM(length(BPM)-1);
BPMmoy=0;

% it create a tab and it go on this tab. For each index, it calculate the
% RR interval and it determine the value in seconds of our RR. it did the
% calcul : 60/ourtime. Then, it have a good value of pulse for each point
% of our tab. it calculate the average in order to have a good BPM. 
for ib=1:length(BPM)
    BPMmoy=BPMmoy+BPM(ib);
end

BPMmoy=BPMmoy/length(BPM);
BPMmax=max(BPM)
BPMmin=min(BPM)
    

% it made many cases to know if there is a bradycardia or tachycardia. 
disp ("The BPM is  :")
disp(BPMmoy)
if BPMmoy < 60 
    disp("The patient has bradycardia")
end
if BPMmoy > 100 
    disp(" The patient has tachycardia")
end

if BPMmoy < 100 && BPMmoy > 60
    disp("The patient has a normal BPM ");
end


% it allow to know if the patient has fibrilation ventricular. it must look
% if the key point have an anarchic position (R,S,Q,T,P). It's the case for
% the second par of the signal of ecgVF.m
if BPMmax >250
    disp("The patient has Ventricul Fibrilation")
end

alpha = 0.2 ;
compteur = 0 ;
for i = 1 : length(Rabs)-2
    if (abs(((Rabs(i+2)-Rabs(i+1))-(Rabs(i+1)-Rabs(i)))/(Rabs(i+2)-Rabs(i+1)))>alpha)
        compteur = compteur + 1 ;
    end
end

%it made a percentage in order to know what if the gravity of ectopic
%pathology. 
pourcentageectopic=compteur/length(R)*100;
disp("the patient has")
disp(pourcentageectopic)
disp("% of pulse that are ectopic")

if pourcentageectopic >10 %it can change this limit 
    disp("the patient has ectopic beats issues")
end
if pourcentageectopic<10
    disp("the patient hasn't ectopic beats issues")
end
somme_atrial = 0 ;
gamma = zeros(N-2,1);
moyennesegment = 60/BPMmoy ;
for k = 1 : N-2
    somme_atrial = 0 ;
    for n = 1 : N-k
        somme_atrial = somme_atrial + ((Rabs(n+k+1)-Rabs(n+k))/Fs-moyennesegment)*((Rabs(n+1)-Rabs(n))/Fs-moyennesegment);
    end
    gamma(k) = (1/(N-k-1))*somme_atrial ;
end

































