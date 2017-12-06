%This code is independent and produced by Lux.
%Email:lixiang19930214@gmail.com
% The input file should be with a standard format. Like:
%latitede,longtitude,date,value
%Value A,Value B, string(normally represents date), Value C
 
clear all;clc;close all;% First part of the code: Insert Data.
name='RB01_Jun_2001';
filepath=strcat(name,'.txt');
[lat,lon,date,pco2]=textread(filepath,'%f%f%s%f','headerlines',0);
 
%Random test
% order_r=randperm(6512);
% lon1=lon(order_r(1:4000));
% lat1=lat(order_r(1:4000));
% pco21=pco2(order_r(1:4000));
% lon2=lon(order_r(4001:6512));
% lat2=lat(order_r(4001:6512));
% pco22=pco2(order_r(4001:6512));
 
%Randomly add up a small displacement for each data's coordinate. In order
%to make autocorrelation matrix has solutions.
num_station=length(pco2);
for i=1:num_station
        lat(i)=lat(i)+0.001*(rand-0.5);
        lon(i)=lon(i)+0.001*(rand-0.5);
end
%Calculate each data pair distance.
for i=1:num_station
    for j=1:num_station
        d(i,j)=((lon(j)-lon(i))^2+((lat(j)-lat(i))^2))^0.5;
    end
end
%Define a resolution for correlation calculation. Normally, each continuous
%data points have a average distance approaching to 0.003. In this case, we set
%it as 0.02, in order to get a sufficient statistics population. For each
%single band. (But we could enlarge them, and it would make a difference to our results.)
%%%
%There is another method to divide its sub-zone. Like 20 gaps.
num_gap=20;
h_resolution=max(max(d)/2)/num_gap;
h_lag=h_resolution;
 
%This part is to calculate the spatial correlation with respect to the
%distance. And R is the autocorrelation function which has been normalized.
sum_total_i_j=0;
for h=0:num_gap
    h_x(h+1)=h*h_resolution;
    if (h==0)
        [x_d y_d]=find(d==0);
    else
        [x_d y_d]=find(d>=h_resolution*h&d<h_resolution*(h+1));
        
    end
    num_d_inside(h+1)=length(x_d);
    sum_total_i_j=sum_total_i_j+length(x_d);
    sum_total=0;
    for i=1:length(x_d)
        sum_total=sum_total+pco2(x_d(i))*pco2(y_d(i));
    end
    Covarance_(h+1)= sum_total/length(x_d)-mean(pco2(x_d))*mean(pco2(y_d));
    R(h+1)= Covarance_(h+1)/((length(x_d)-1)/length(x_d)*var(pco2(x_d))*...
        (length(x_d)-1)/length(x_d)*var(pco2(y_d)))^0.5;
end
 
% Curve fit, and estimate the radius of circular.
[a]=polyfit(h_x,R,4);
h_r=[0:0.0001:max(max(d)/3)]; 
% The reason why I use max(d)/3 is that prevent the fitted curve have two potiential answer.
R_fit=a(1)*h_r.^4+a(2)*h_r.^3+a(3)*h_r.^2+a(4)*h_r+a(5);
[b]=polyfit(h_x,Covarance_,4);
C_fit=b(1)*h_r.^4+b(2)*h_r.^3+b(3)*h_r.^2+b(4)*h_r+b(5);
% plot(h_r, R_fit,'r');
% hold on;
% scatter(h_x, R);
% title('Autocorrelation');
 
%This part try to find the coordinate x which R_fit equals to 0.3125*R(0) to represent Radius.
while 1
    [order]=find(R_fit<(0.3125+h_resolution)*(max(R_fit))&R_fit>(0.3125-h_resolution)*(max(R_fit)));
    if (length(order)==0)
        h_resolution=2*h_resolution;
    elseif(length(order)>2)
        h_resolution=0.5*h_resolution;
    else
        break;
    end
end
if (length(order)==1)
    radius=h_r(order);
else
    radius=mean(h_r(order));
end
 
%After we obtaining the radius and we will know the range,nugget effect
%value and drill value.
%Clear the previous data, keep data to calculate range,nugget effect
%value and drill value.
clearvars -except radius R lat lon pco2 d num_station Covarance_ C_fit h_r h_x
%C=C0+C   a represents range.
min_long=min(lon)-2*radius;max_long=max(lon)+2*radius;min_lat=min(lat)-2*radius; max_lat=max(lat)+2*radius;
aera=(max_lat-min_lat)*(max_long-min_long);
amount=0;
C=Covarance_(1);
a=2*radius;
C0=C-C_fit(1);
% After we obtaining model parameters, we start to do the Kriging
% Estimations next.
x_c=min(lon);%start point
y_c=min(lat);
x_size=round((max(lon)-min(lon))/(radius/8)); %Kriging Interpolation resolution equals to radius/8
y_size=round((max(lat)-min(lat))/(radius/8));
 
for m=1:x_size
    for n=1:y_size
        order_output=n+(m-1)*y_size;
        x_c(order_output)=min(lon)+(m-1/2)*(radius/8);
        y_c(order_output)=min(lat)+(n-1/2)*(radius/8);
        C_0=0;C_=0;d_0=0;d=0;
        %Point Kriging interpolation.
        for i=1:length(lat)%Calculation estimating point's distances among other points.
            d_0(i)= ((lon(i)-x_c(order_output))^2+((lat(i)-y_c(order_output))^2))^0.5;
        end
        [d_temp]=find(d_0<=a);%Reserve the data points which are included inside the range.
        display(length(d_temp));
        Kriging_points=length(d_temp);% A single estimation which are determined by how many data points.
        
        %Construc the Autocorrelation Matrix. Ci,j
        if (Kriging_points>0) %If the estimating point have more than one data in its range.
            d_0=d_0(d_temp);pco2_temp=pco2(d_temp);
            lon_temp=lon(d_temp);lat_temp=lat(d_temp);
            
            for i=1:Kriging_points
                C_0(i,1)=Covarance_(1)*[1.5*(1-d_0(i)/a).^2-0.5*(1-d_0(i)/a).^3]; %Autocorrelation function
                %Spherical model
            end
                C_0(Kriging_points+1,1)=1;
            
            for i=1:Kriging_points
                for j=1:Kriging_points
                    d(i,j)=((lon_temp(j)-lon_temp(i))^2+((lat_temp(j)-lat_temp(i))^2))^0.5;
                    C_(i,j)=Covarance_(1)*[(d(i,j)/a)-0.5*(1-d(i,j)/a).^3];
                end
            end
            
            for i=1:Kriging_points
                C_(i,Kriging_points+1)=1;
                C_(Kriging_points+1,i)=1;
            end
            
            %The detail of constructing the metrix will be show in P88.
            %Solve the equations.
            lamda=inv(C_)*C_0;
            pco2_e(order_output)=0;
            uncertainty_e(order_output)=Covarance_(1)+lamda(Kriging_points+1);
            
            for i=1:Kriging_points
                pco2_e(order_output)=lamda(i)*pco2_temp(i)+pco2_e(order_output);
                uncertainty_e(order_output)=uncertainty_e(order_output)-lamda(i)*C_0(i);
            end
            
        else       
%             pco2_e(order_output)=-999;
%             uncertainty_e(order_output)=-999;
        end
 
    end
end
%Draw an research region with data distribution.
order=find(pco2_e>0);
pco2_e=pco2_e(order);
uncertainty_e=abs(uncertainty_e(order));y_c=y_c(order);x_c=x_c(order);
uncertainty=sqrt(uncertainty_e);
m_proj('mercator','lat',[min_long,max_long],'long',[min_lat,max_lat]);
[x y]=m_ll2xy(x_c,y_c,'clip');
scatter(x,y,30,pco2_e,'filled');
m_gshhs_f('patch',[0.7 0.7 0.7],'edgecolor','none');
m_grid('box','fancy','tickdir','in');
colorbar;
title('Uncertainty','fontsize',14);
 
% test_num_Men=10000;
% %Using mento carlo to estimate the aera of observations.
% for i=1:test_num_Men;
%     x_test=rand*(max_long-min_long)+min_long;
%     y_test=rand*(max_lat-min_lat)+min_lat;
%     for j=1:num_station
%         d_test=(abs((lon(j)-x_test))^2+(abs(lat(j)-y_test))^2)^0.5;
%         if(d_test<=radius)
%             amount=amount+1;
%             break;
%         end
%     end
% end
%
% eff_obser=(amount/test_num_Men)*aera/(pi*radius^2);
% uncertainty_sample= var(pco2)/eff_obser;
