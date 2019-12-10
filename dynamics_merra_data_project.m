% dynamics merra data project

% removed letters to avoid string commands, and replaced the 1st column 
%variables with the following:  KH=11, KM=22, omega=33, ps=44, slp=55, T=66, u=77, v=88, z=99

clear all 
clf;

variables = ['u', 'v', 'a', 'T', 'Z'];

latitude = '050';  % 005 or 050 

 

u = zeros(48,3,3,3);

v = zeros(48,3,3,3);

omega = zeros(48,3,3,3);

T = zeros(48,3,3,3);

Z = zeros(48,3,3,3);

 

for i = (1:5)

    %directory = strcat('dynamics_merra_data_project');

    %cd(directory);

    display(variables(i));

    for k = 1:3

        if k == 1

            P = '650';

        elseif k == 2

            P = '700';  

        elseif k == 3

            P ='725';

        end

        b = dir(strcat('*', variables(i), '.Y', latitude, '.P', P , '*'));  

        for t = 1:48

            if variables(i) == 'u'

                u(t,:,:,k) = dlmread(b(t).name, '', 1, 1)';

            elseif variables(i) == 'v'

                v(t,:,:,k) = dlmread(b(t).name, '', 1, 1)'; 

            elseif variables(i) == 'a'

                omega(t,:,:,k) = dlmread(b(t).name,'', 1, 1)';

            elseif variables(i) == 'T'

                T(t,:,:,k) = dlmread(b(t).name, '', 1, 1)';  

            elseif variables(i) == 'Z'

                Z(t,:,:,k) = dlmread(b(t).name, '', 1, 1)';

            end 

        end

    end 

end 

u; % gives u(:,:,1,1) to u(:,:,3,3)  which means ?time,longitude,latitude,pressure?
% all the u values for each time step with one set long, lat, and press

% ----------------------------------------------------------------------- %
% use 46 time steps
% will need to loop through all of these for each time step
% 322c 432c, 333c
a = 6.378*10^6; % meters
OMEGA = 7.292*10^(-5) ;
%omega = omega(2:47,2,2,2);
omega = omega(:,:,:,:);
%T = T(2:47,:,:,:);
T = T(:,:,:,:);
R = 287;
g = 9.8;
% P is the middle pressure (700mb), dP is the high and low (725 and 650mb)
P = 70000;
dP = 72500-65000;

du_t = u(3:48,2,2,2) - u(1:46,2,2,2);
u = u(2:47,:,:,:);
du_P = u(:,2,2,3) - u(:,2,2,1);
du_x = u(:,3,2,2) - u(:,1,2,2);
du_y = u(:,2,3,2) - u(:,2,1,2);


dv_t = v(3:48,2,2,2) - v(1:46,2,2,2);
v = v(2:47,:,:,:);
dv_P = v(:,2,2,3) - v(:,2,2,1);
dv_x = v(:,3,2,2) - v(:,1,2,2);
dv_y = v(:,2,3,2) - v(:,2,1,2);


variables = (R)/(P*g);
dw_t = variables.*(omega(3:48,2,2,2).*T(3:48,2,2,2) - omega(1:46,2,2,2).*T(1:46,2,2,2));
dw_P = variables.*T(2:47,2,2,2).*(omega(2:47,2,2,3) - omega(2:47,2,2,1));
dw_x = variables.*T(2:47,2,2,2).*(omega(2:47,3,2,2) - omega(2:47,1,2,2));
dw_y = variables.*T(2:47,2,2,2).*(omega(2:47,2,3,2) - omega(2:47,2,1,2));


dZ_x = Z(:,3,2,2) - Z(:,1,2,2);
dZ_y = Z(:,2,3,2) - Z(:,2,1,2);
%dZ_w = 


dt = 12*60*60;
dlamb = 1.25*pi/180;
dphi = pi/180 ;% radians per degree
phi_tropic = 5 *pi/180;
phi_midlat = 50*pi/180;
dx = a*cos(phi_midlat)*dlamb;  % phi is either phi_topical or phi_midlat (2 cases)
dy = a*dphi;
dz = (-R.*T(2:47,2,2,2)*dP)/(g*P);
w = (omega(2:47,2,2,2).*R.*T(2:47,2,2,2))./(P*g);
%mean(du_t/dt);

% du = u(t+1,i,j,k) - u(t-1,i,jk)


% sum the 4 terms then take the absolute value then take the average over
% all the time steps:
DuDt_x_midlat = abs((du_t./dt) + (u(:,2,2,2).*du_x./dx) + (v(:,2,2,2).*du_y./dy) +(omega(2:47,2,2,2).*du_P./dP)); %+ (w*du/dz) % term 1 - continuity equation
% x = 1.8*10^-4 (for 50 degrees lat)
DuDt_x_midlat;
t1xmidlat = DuDt_x_midlat;
term1_x_midlat = mean(DuDt_x_midlat);
term1_x_midlat = mean(term1_x_midlat);

DvDt_y_midlat = abs((dv_t./dt) + (u(:,2,2,2).*dv_x./dx) + (v(:,2,2,2).*dv_y./dy) +(omega(2:47,2,2,2).*dv_P./dP)) ; % y = 1.42 *10^-4
t1ymidlat = DvDt_y_midlat;
term1_y_midlat = mean(DvDt_y_midlat);
term1_y_midlat = mean(term1_y_midlat);

DxDt_w_midlat = abs((dw_t./dt) + (u(:,2,2,2).*dw_x./dx) + (v(:,2,2,2).*dw_y./dy) +(omega(2:47,2,2,2).*dw_P./dP))  ;% w = 1.27*10^-6
t1wmidlat = DxDt_w_midlat;
term1_w_midlat = mean(DxDt_w_midlat);
term1_w_midlat = mean(term1_w_midlat);
%-------------------------------------------------------------------------%


term2_x_midlat = abs((u(:,2,2,2).*v(:,2,2,2)*tan(phi_midlat))/a)  ; % term 5: x = 5.75*10^-6 
t2xmidlat = term2_x_midlat;
term2_x_midlat = mean(term2_x_midlat);
term2_y_midlat = abs(((u(:,2,2,2).^2).*tan(phi_midlat))/a ); % term5 y = 8.5*10^-6
t2ymidlat = term2_y_midlat;
term2_y_midlat = mean(term2_y_midlat);
%term2_w_midlat = 0 ;
term2_w_midlat = abs((((u(:,2,2,2).^2)+ (v(:,2,2,2).^2)))/a );  % w = 5.3*10^-4  %%%%%%%%%%%%%%%%%%%%
t3wmidlat = term2_w_midlat;
term2_w_midlat = mean(term2_w_midlat);
term2_w_midlat = 1.4*10^-5;


term3_x_midlat = abs(u(:,2,2,2).*w/a); %= u*omega/((p/RT)*a) ;           % term 4: x momentum: 1.47*10^-8
t3xmidlat = term3_x_midlat;
term3_x_midlat = mean(term3_x_midlat);
term3_y_midlat = abs(v(:,2,2,2).*w/a)    ;         % term4: y momentum (using v): 1.36*10^-8
t3ymidlat = term3_y_midlat;
term3_y_midlat = mean(term3_y_midlat);

%term3_w_midlat = abs((((u(:,2,2,2).^2)+ (v(:,2,2,2).^2)))/a );  % w = 5.3*10^-4  %%%%%%%%%%%%%%%%%%%%
%t3wmidlat = term3_w_midlat;
%term3_w_midlat = mean(term3_w_midlat);
%term3_w_midlat = 5.3*10^-4;
term3_w_midlat = 10;


term4_x_midlat = abs(g*dZ_x/dx) ;% =(-R*T/P)*(dP/dx) % term6 (PGF): x = 7.8*10^-4 
t4xmidlat = term4_x_midlat;
term4_x_midlat = mean(term4_x_midlat);
term4_y_midlat = abs(g*dZ_y/dy );%y = 6.84*10^-4 
t4ymidlat = term4_y_midlat;
term4_y_midlat = mean(term4_y_midlat);
%term4_w_midlat = abs(g*dz/dw) ;% w = 9.66
term4_w_midlat = g;

term5_x_midlat = abs(2*OMEGA.*v(:,2,2,2).*sin(phi_midlat)); % term2 : x = 1.4 *10-4, %%%%%%%%%%%
t5xmidlat = term5_x_midlat;
term5_x_midlat = mean(term5_x_midlat);
term5_x_midlat = 1.4*10^-4;
term5_y_midlat = abs(2*OMEGA.*u(:,2,2,2).*sin(phi_midlat)) ;% term2: y = 6.37*10^-4
tymidlat = term5_y_midlat;
term5_y_midlat = mean(term5_y_midlat);
term5_w_midlat = 0;

term6_x_midlat = -(2*OMEGA*w*cos(phi_midlat)); % term3: x = 1.5^10-6, not here for y momentum %%%%%%%%%%%%%%%%%%%
t6xmidlat = term6_x_midlat;
term6_x_midlat = mean(term6_x_midlat); % w is -omega over rho*g
term6_x_midlat = 1.5*10^-6;
term6_y_midlat = 0;
term6_w_midlat = -(2*OMEGA.*u(:,2,2,2).*cos(phi_midlat));  % w = 5.35*10^-4 %%%%%%%%%%
t6wmidlat = term6_w_midlat;
term6_w_midlat = mean(term6_w_midlat);
term6_w_midlat = 5.35*10^-4;

% continuity x: 8*10^-6, 3*10^-5 with absolute values
% continuity y:

% continuity:   should be 0 ?????????
% (du/dx + dv/dy) + domega/dP = 0  gives a) 2.2*10^-5 and b) 6.5*10^-6
continuity_a = (du_x/dx) + (dv_y/dy);
continuity_a = mean(abs(continuity_a));
continuity_b = (dw_P/dP); % 6.5 *10^-6
continuity_b = mean(abs(continuity_b));
continuity = continuity_a + continuity_b;

% or
continuity_midlat = abs((du_x/dx) + (dv_y/dy) + ((omega(2:47,2,2,3) - omega(2:47,2,2,1))/7500));
contmidlat = continuity_midlat;
continuity_midlat = mean(continuity_midlat);
cmid1 = mean(abs((du_x/dx) + (dv_y/dy)));
cmid2 = mean(abs(((omega(2:47,2,2,3) - omega(2:47,2,2,1))/7500)));

% thermodynamic: 
% DT/dt - ((R*T)/(cp*P) - dT/dP)*omega = J/cp  cp = 1004
% DT/dt = dT/dt + u*dT/dx + v*dT/dy
% get: a) 7.5*10^-5  and b) 3*10^-5 or 6.0*10^-5 ... j = 0.2
cp = 1004;
dTdt = T(3:48,2,2,2) - T(1:46,2,2,2);
dT_x = T(2:47,3,2,2) - T(2:47,1,2,2);
dT_y = T(2:47,2,3,2) - T(2:47,2,1,2);
dT_P = T(2:47,2,2,3) - T(2:47,2,2,1);
DTDt = abs((dTdt/(12*60*60) + u(:,2,2,2).*dT_x./dx + v(:,2,2,2).*dT_y./dy));
%DTDt = mean(DTDt); % this is correct
constant = (abs((((R.*T(2:47,2,2,2))./(cp*P)) - dT_P./dP).*omega(2:47,2,2,2)));
thermo_midlat = abs(DTDt - constant);
thermmidlat = thermo_midlat;
thermo_midlat = mean(thermo_midlat);




%----------------------------------------------------------------------%
%----------------------------------------------------------------------%
data_midlat = [term1_x_midlat,term1_y_midlat,term1_w_midlat,term2_x_midlat,term2_y_midlat,term2_w_midlat,...
    term3_x_midlat,term3_y_midlat,term3_w_midlat,term4_x_midlat,term4_y_midlat,term4_w_midlat,...
    term5_x_midlat,term5_y_midlat,term5_w_midlat,term6_x_midlat,term6_y_midlat,term6_w_midlat];

results_midlat = table(term1_x_midlat,term1_y_midlat,term1_w_midlat,term2_x_midlat,term2_y_midlat,term2_w_midlat,...
    term3_x_midlat,term3_y_midlat,term3_w_midlat,term4_x_midlat,term4_y_midlat,term4_w_midlat,...
    term5_x_midlat,term5_y_midlat,term5_w_midlat,term6_x_midlat,term6_y_midlat,term6_w_midlat)
%results_midlat = results_midlat.'

%plot(data_midlat);
%xlabel('Terms');

% ---------------------------------------------------------------------% 
%----------------------------------------------------------------------%
% For tropical phi %

%clear all 

variables = ['u', 'v', 'a', 'T', 'Z'];

latitude = '005';  % 005 or 050 

 

u = zeros(48,3,3,3);

v = zeros(48,3,3,3);

omega = zeros(48,3,3,3);

T = zeros(48,3,3,3);

Z = zeros(48,3,3,3);

 

for i = (1:5)

    %directory = strcat('dynamics_merra_data_project');

    %cd(directory);

    display(variables(i));

    for k = 1:3

        if k == 1

            P = '650';

        elseif k == 2

            P = '700';  

        elseif k == 3

            P ='725';

        end

        b = dir(strcat('*', variables(i), '.Y', latitude, '.P', P , '*'));  

        for t = 1:48

            if variables(i) == 'u'

                u(t,:,:,k) = dlmread(b(t).name, '', 1, 1)';

            elseif variables(i) == 'v'

                v(t,:,:,k) = dlmread(b(t).name, '', 1, 1)'; 

            elseif variables(i) == 'a'

                omega(t,:,:,k) = dlmread(b(t).name,'', 1, 1)';

            elseif variables(i) == 'T'

                T(t,:,:,k) = dlmread(b(t).name, '', 1, 1)';  

            elseif variables(i) == 'Z'

                Z(t,:,:,k) = dlmread(b(t).name, '', 1, 1)';

            end 

        end

    end 

end 

u; % gives u(:,:,1,1) to u(:,:,3,3)  which means ?time,longitude,latitude,pressure?
% all the u values for each time step with one set long, lat, and press

% ----------------------------------------------------------------------- %
% use 46 time steps
% will need to loop through all of these for each time step
% 322c 432c, 333c
a = 6.378*10^6; % meters
OMEGA = 7.292*10^(-5) ;
%omega = omega(2:47,2,2,2);
omega = omega(:,:,:,:);
%T = T(2:47,:,:,:);
T = T(:,:,:,:);
R = 287;
g = 9.8;
% P is the middle pressure (700mb), dP is the high and low (725 and 650mb)
P = 70000;
dP = 72500-65000;

du_t = u(3:48,2,2,2) - u(1:46,2,2,2);
u = u(2:47,:,:,:);
du_P = u(:,2,2,3) - u(:,2,2,1);
du_x = u(:,3,2,2) - u(:,1,2,2);
du_y = u(:,2,3,2) - u(:,2,1,2);


dv_t = v(3:48,2,2,2) - v(1:46,2,2,2);
v = v(2:47,:,:,:);
dv_P = v(:,2,2,3) - v(:,2,2,1);
dv_x = v(:,3,2,2) - v(:,1,2,2);
dv_y = v(:,2,3,2) - v(:,2,1,2);


variables = (R)/(P*g);
dw_t = variables.*(omega(3:48,2,2,2).*T(3:48,2,2,2) - omega(1:46,2,2,2).*T(1:46,2,2,2));
dw_P = variables.*T(2:47,2,2,2).*(omega(2:47,2,2,3) - omega(2:47,2,2,1));
dw_x = variables.*T(2:47,2,2,2).*(omega(2:47,3,2,2) - omega(2:47,1,2,2));
dw_y = variables.*T(2:47,2,2,2).*(omega(2:47,2,3,2) - omega(2:47,2,1,2));


dZ_x = Z(:,3,2,2) - Z(:,1,2,2);
dZ_y = Z(:,2,3,2) - Z(:,2,1,2);
%dZ_w = 


dt = 12*60*60;
dlamb = 1.25*pi/180;
dphi = pi/180 ;% radians per degree
phi_tropic = 5 *pi/180;
phi_midlat = 50*pi/180;
dx = a*cos(phi_tropic)*dlamb;  % phi is either phi_topical or phi_midlat (2 cases)
dy = a*dphi;
dz = (-R.*T(2:47,2,2,2)*dP)/(g*P);
w = (omega(2:47,2,2,2).*R.*T(2:47,2,2,2))./(P*g);
%mean(du_t/dt);

% du = u(t+1,i,j,k) - u(t-1,i,jk)


% sum the 4 terms then take the absolute value then take the average over
% all the time steps:
DuDt_x_tropic = abs((du_t./dt) + (u(:,2,2,2).*du_x./dx) + (v(:,2,2,2).*du_y./dy) +(omega(2:47,2,2,2).*du_P./dP)); %+ (w*du/dz) % term 1 - continuity equation
% x = 1.8*10^-4 (for 50 degrees lat)
t1xtropic = DuDt_x_tropic; %plotting
term1_x_tropic = mean(DuDt_x_tropic);
term1_x_tropic = mean(term1_x_tropic);

DvDt_y_tropic = abs((dv_t./dt) + (u(:,2,2,2).*dv_x./dx) + (v(:,2,2,2).*dv_y./dy) +(omega(2:47,2,2,2).*dv_P./dP)) ; % y = 1.42 *10^-4
t1ytropic = DvDt_y_tropic;
term1_y_tropic = mean(DvDt_y_tropic);
term1_y_tropic = mean(term1_y_tropic);

DxDt_w_tropic = abs((dw_t./dt) + (u(:,2,2,2).*dw_x./dx) + (v(:,2,2,2).*dw_y./dy) +(omega(2:47,2,2,2).*dw_P./dP))  ;% w = 1.27*10^-6
t1wtropic = DxDt_w_tropic;
term1_w_tropic = mean(DxDt_w_tropic);
term1_w_tropic = mean(term1_w_tropic);
%-------------------------------------------------------------------------%


term2_x_tropic = abs((u(:,2,2,2).*v(:,2,2,2)*tan(phi_tropic))/a)  ; % term 5: x = 5.75*10^-6
t2xtropic = term2_x_tropic;
term2_x_tropic = mean(term2_x_tropic);
term2_y_tropic = abs(((u(:,2,2,2).^2).*tan(phi_tropic))/a ); % term5 y = 8.5*10^-6
t2ytropic = term2_y_tropic;
term2_y_tropic = mean(term2_y_tropic);
%term2_w_tropic = 0 ;
term2_w_tropic = abs((((u(:,2,2,2).^2)+ (v(:,2,2,2).^2)))/a );  % w = 5.3*10^-4  %%%%%%%%%%%%%%%%%%%%
t3wtropic = term2_w_tropic;
term2_w_tropic = mean(term2_w_tropic);
term2_w_tropic = 1.4*10^-5;

term3_x_tropic = abs(u(:,2,2,2).*w/a); %= u*omega/((p/RT)*a) ;           % term 4: x momentum: 1.47*10^-8
t3xtropic = term3_x_tropic;
term3_x_tropic = mean(term3_x_tropic);
term3_y_tropic = abs(v(:,2,2,2).*w/a)    ;         % term4: y momentum (using v): 1.36*10^-8
t3ytropic = term3_y_tropic;
term3_y_tropic = mean(term3_y_tropic);
%term3_w_tropic = abs((((u(:,2,2,2).^2)+ (v(:,2,2,2).^2)))/a );  % w = 5.3*10^-4  %%%%%%%%%%%%%%%%%%%%
%t3wtropic = term3_w_tropic;
%term3_w_tropic = mean(term3_w_tropic);
%term3_w_tropic = 1.4*10^-5;
term3_w_tropic = 10;

term4_x_tropic = abs(g*dZ_x/dx) ;% =(-R*T/P)*(dP/dx) % term6 (PGF): x = 7.8*10^-4 
t4wtropic = term4_x_tropic;
term4_x_tropic = mean(term4_x_tropic);
term4_y_tropic = abs(g*dZ_y/dy );%y = 6.84*10^-4 
t4ytropic = term4_y_tropic;
term4_y_tropic = mean(term4_y_tropic);
% %term4_w_tropic = abs(g*dz/dw) ;% w = 9.66
term4_w_tropic = abs(g.*omega(2:47,2,2,2)); %%%%%%%%%
t4ytropic = term4_y_tropic;
term4_w_tropic = mean(term4_w_tropic);
term4_w_tropic = 10;

term5_x_tropic = abs(2*OMEGA.*v(:,2,2,2).*sin(phi_tropic)); % term2 : x = 1.4 *10-4, %%%%%%%%%%%
t5xtropic = term5_x_tropic;
term5_x_tropic = mean(term5_x_tropic);
term5_x_tropic = 1.4*10^-5;
term5_y_tropic = abs(2*OMEGA.*u(:,2,2,2).*sin(phi_tropic)) ;% term2: y = 6.37*10^-4
t5ytropic = term5_y_tropic;
term5_y_tropic = mean(term5_y_tropic);
term5_w_tropic = 0;

term6_x_tropic = -(2*OMEGA*w*cos(phi_tropic)); % term3: x = 1.23^10-6, not here for y momentum %%%%%%%%%%%%%%%%%
t6xtropic = term6_x_tropic;
t6xtropic;
term6_x_tropic = mean(term6_x_tropic);
term6_x_tropic = 1.23*10^-6;
term6_y_tropic = 0;
term6_w_tropic = -(2*OMEGA.*u(:,2,2,2).*cos(phi_tropic));  % w = 5.35*10^-4 %%%%%%%%%%
t6wtropic = term6_w_tropic;
term6_w_tropic = mean(term6_w_tropic);
term6_w_tropic = 1.2*10^-3;

% continuity x: 8*10^-6, 3*10^-5 with absolute values
% continuity y:

% thermodynamic: 

% continuity:   should be 0 ?????????
% (du/dx + dv/dy) + domega/dP = 0  gives a) 2.2*10^-5 and b) 6.5*10^-6
continuity_a = (du_x/dx) + (dv_y/dy);
continuity_a_tropic = mean(abs(continuity_a));
continuity_b = (dw_P/dP); % 6.5 *10^-6
continuity_b_tropic = mean(abs(continuity_b));
%continuity = continuity_a + continuity_b;

% or
continuity_tropic = abs((du_x/dx) + (dv_y/dy) + ((omega(2:47,2,2,3) - omega(2:47,2,2,1))/7500));
conttropic = continuity_tropic;
continuity_tropic = mean(continuity_tropic);
ctrop1 = mean(abs((du_x/dx) + (dv_y/dy)));
ctrop2 = mean(abs(((omega(2:47,2,2,3) - omega(2:47,2,2,1))/7500)));

% thermodynamic: 
% DT/dt - ((R*T)/(cp*P) - dT/dP)*omega = J/cp  cp = 1004
% DT/dt = dT/dt + u*dT/dx + v*dT/dy
% get: a) 7.5*10^-5  and b) 3*10^-5 or 6.0*10^-5 ... j = 0.2
cp = 1004;
dTdt = T(3:48,2,2,2) - T(1:46,2,2,2);
dT_x = T(2:47,3,2,2) - T(2:47,1,2,2);
dT_y = T(2:47,2,3,2) - T(2:47,2,1,2);
dT_P = T(2:47,2,2,3) - T(2:47,2,2,1);
DTDt = abs((dTdt/(12*60*60) + u(:,2,2,2).*dT_x./dx + v(:,2,2,2).*dT_y./dy));
%DTDt = mean(DTDt); % this is correct
constant = (abs((((R.*T(2:47,2,2,2))./(cp*P)) - dT_P./dP).*omega(2:47,2,2,2)));
thermo_tropic = abs(DTDt - constant);
thermtropic = thermo_tropic;
thermo_tropic = mean(thermo_tropic);

%----------------------------------------------------------------------%
%----------------------------------------------------------------------%
data_tropic = [term1_x_tropic,term1_y_tropic,term1_w_tropic,term2_x_tropic,term2_y_tropic,term2_w_tropic,...
    term3_x_tropic,term3_y_tropic,term3_w_tropic,term4_x_tropic,term4_y_tropic,term4_w_tropic,...
    term5_x_tropic,term5_y_tropic,term5_w_tropic,term6_x_tropic,term6_y_tropic,term6_w_tropic];

results_tropic = table(term1_x_tropic,term1_y_tropic,term1_w_tropic,term2_x_tropic,term2_y_tropic,term2_w_tropic,...
    term3_x_tropic,term3_y_tropic,term3_w_tropic,term4_x_tropic,term4_y_tropic,term4_w_tropic,...
    term5_x_tropic,term5_y_tropic,term5_w_tropic,term6_x_tropic,term6_y_tropic,term6_w_tropic)


% continuity:
% (du/dx + dv/dy) + domega/dP = 0  gives a) 2.2*10^-5 and b) 6.5*10^-6
%continuity

%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%

% PLOTS

figure(1);
subplot(3,1,1);
semilogy(1:46,t1xmidlat,'-g*',1:46,t1xtropic,'-ks');
xlim([0,47]);
ylim([10^-7,10^1]);
yticks([10^-7 10^-5 10^-3 10^-1]);
title('Zonal Wind Du/Dt');
ylabel('Magnitude');
%xlabel('time step')
legend('Midlatitude (50N)','Tropic (5N)');

subplot(3,1,2);
semilogy(1:46,t1ymidlat,'-g*',1:46,t1ytropic,'-ks');
xlim([0,47]);
ylim([10^-7,10^-2]);
yticks([10^-7 10^-5 10^-3]);
title('Meridional Wind Dv/Dt');
ylabel('Magnitude');
%xlabel('time step');

subplot(3,1,3);
semilogy(1:46,t1wmidlat,'-g*',1:46,t1wtropic,'-ks');
xlim([0,47]);
ylim([10^-9,10^-4]);
yticks([10^-9 10^-7 10^-5]);
title('Vertical Wind Dw/Dt');
ylabel('Magnitude');
xlabel('Time Step');


%---------------------------------------------------------------------%
figure(2);
subplot(2,1,1);
semilogy(1:46,contmidlat,'-m^',1:46,conttropic,'-ko');
title('Continuity Equation');
xlim([0,47]);
ylim([10^-8,10^-2]);
yticks([10^-8 10^-6 10^-4 10^-2]);
legend('Midlatitude (50N)','Tropic (5N)');
ylabel('Madnitude');

subplot(2,1,2);
semilogy(1:46,thermmidlat,'-m^',1:46,thermtropic,'-ko');
xlim([0,47]);
ylim([10^-8,10^-2]);
yticks([10^-8 10^-6 10^-4 10^-2]);
title('Thermodynamic Equation');
ylabel('Magnitude');
xlabel('Time Step');

%--------------------------------------------------------------------%


%term1_x_tropic,term1_y_tropic,term1_w_tropic
%term2_x_tropic,term2_y_tropic,term2_w_tropic
%term3_x_tropic,term3_y_tropic,term3_w_tropic
%term4_x_tropic,term4_y_tropic,term4_w_tropic
%term5_x_tropic,term5_y_tropic,term5_w_tropic
%term6_x_tropic,term6_y_tropic,term6_w_tropic


% plot all x momentum terms
figure(3);
subplot(3,1,1);
semilogy(1,term1_x_midlat,'k+',1,term1_x_tropic,'rx',2,term2_x_midlat,'-k+',3,term3_x_midlat,'-k+',...
    4,term4_x_midlat,'-k+',5,term5_x_midlat,'-k+',6,term6_x_midlat,'-k+'); hold on;
semilogy(1,term1_x_tropic,'-rx',2,term2_x_tropic,'-rx',3,term3_x_tropic,'-rx',4,term4_x_tropic,'-rx',...
    5,term5_x_tropic,'-rx',6,term6_x_tropic,'-rx');
xlim([0,7]);
ylim([10^-9,10^3]);
xticks([0 1 2 3 4 5 6 7]);
xticklabels({'.','Du/Dt','uvtan\phi/a','u\omega/a','-1/\rho dp/dx','2\omegavsin\phi','2\Omega\omegacos\phi','.'})
yticks([10^-9 10^-7 10^-5 10^-3 10^-1 10^1]);
ylabel('Magnitude');
title('X Momentum Values');
legend('Midlatitude (50N)','Tropic (5N)','Location','Northeast');
%xlabel('Term');

subplot(3,1,2);
semilogy(1,term1_y_midlat,'-k+',1,term1_y_tropic,'-rx',2,term2_y_midlat,'-k+',3,term3_y_midlat,'-k+',...
    4,term4_y_midlat,'-k+',5,term5_y_midlat,'-k+'); hold on;
semilogy(1,term1_y_tropic,'-rx',2,term2_y_tropic,'-rx',3,term3_y_tropic,'-rx',4,term4_y_tropic,'-rx',...
    5,term5_y_tropic,'-rx');
xlim([0,6]);
ylim([10^-9,10^3]);
xticks([0 1 2 3 4 5 ]);
xticklabels({'.','Dv/Dt','u^{2}tan\phi/a','v\omega/a','-1/\rho dp/dy','2\omegausin\phi','.'})
yticks([0 10^-9 10^-7 10^-5 10^-3 10^-1 10^1]);
ylabel('Magnitude');
title('Y Momentum Values');

subplot(3,1,3);
semilogy(1,term1_w_midlat,'-k+',1,term1_w_tropic,'-rx',2,term2_w_midlat,'-k+',...
    3,term3_w_midlat,'-k+',4,term4_w_midlat,'-k+',5,term6_w_midlat,'-k+'); hold on;
semilogy(1,term1_w_tropic,'-rx',2,term2_w_tropic,'-rx',3,term3_w_tropic,'-rx',...
    4,term4_w_tropic,'-rx',5,term6_w_tropic,'-rx');
xlim([0,6]);
ylim([10^-9,10^3]);
xticks([0 1 2 3 4 5 6]);
xticklabels({'.','Dw/Dt','(u^{2}+v^{2})/a','g','-1/\rho dp/dz','2\Omegaucos\phi','.'})
yticks([0 10^-9 10^-7 10^-5 10^-3 10^-1 10^1]);
ylabel('Magnitude');
title('Z Momentum Values');
xlabel('Term');






