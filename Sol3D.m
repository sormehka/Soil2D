function varargout = Sol3D(varargin)
% SOL3D M-file for Sol3D.fig


% Programmer:  Sormeh Kashef Haghighi
%      SOL3D, by itself, creates a new SOL3D or raises the existing
%      singleton*.
%
%      H = SOL3D returns the handle to a new SOL3D or the handle to
%      the existing singleton*.
%
%      SOL3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SOL3D.M with the given input arguments.
%
%      SOL3D('Property','Value',...) creates a new SOL3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Sol3D_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Sol3D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Sol3D

% Last Modified by GUIDE v2.5 29-Oct-2006 17:25:42

% Inorder to be compatible with MATLAB 7.1: In the GUIDE tool menue, GUI
% options, Resize behaviour ---> Proportional  The problem before changing
% this was the resizing of 3D image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input parameters:

% v1 -- Darcy velocity in x direction  (m/s)
% n   --  effective porosity (dimensionless)
% bd  --  dry bulk density (g/cm^3)
% ax,ay,az  -- dispersivity (m)
% dm  -- molecular diffusion coefficient (m2/s) 
% d   -- decay constant (1/s)
% foc -- fraction organic carbon (%/100)
% koc -- organic carbon coefficient (cm3/g)
% t   -- snapshot time (s)
% m   -- Injected mass (kg)
% te  -- animation end time (s)
% ts  -- animation step time (s)

% Calculated parametrs:
% r -- retardation factor
% v1n -- Seepage velocity (m/s)
% dx,dy,dz -- hydrodynamic dispersion coefficient (m2/s)
% Kd -- adsorption coefficient (cm3/g)
% lx,ly,lz -- plume dimension in x,y,z 
% maxc -- maximum concentration (mg/L)



if nargin == 0  % LAUNCH GUI
    
    fig = openfig(mfilename,'reuse');
    
    % Use system color scheme for figure:
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    
    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    guidata(fig, handles);
    
    if nargout > 0
        varargout{1} = fig;
    end
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    try
        [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
    catch
        %disp(lasterr);
    end
end

clear all;

%-----CALCULATE INPUT PARAMETERS----------------------------------------------------
%logical sequence as follows: (1) editbox callback (2) find applicable input values 
%(3) calculate applicable parameters, and (4) display value

%for v: average linear velocity and hydrodynamic dispersion coeff.

function varargout = v1input_Callback(h, eventdata, handles, varargin)
v1= str2double(get(handles.v1input,'String'));
n = str2double(get(handles.ninput,'String'));
v1n=v1/n;
set (handles.v1ncalculation,'string',v1n);
ax = str2double(get(handles.axinput,'String'));
dm = str2double(get(handles.dminput,'String'));
dx=((ax*v1n)+dm);
set (handles.dxcalculation,'string',dx);


%for n: average linear velocity, hydrodynamic dispersion coeff. and retardation factor
function varargout = ninput_Callback(h, eventdata, handles, varargin)
n = str2double(get(handles.ninput,'String'));
v1= str2double(get(handles.v1input,'String'));
v1n=v1/n;
set (handles.v1ncalculation,'string',v1n);
ax = str2double(get(handles.axinput,'String'));
dm = str2double(get(handles.dminput,'String'));
dx=((ax*v1n)+dm);
set (handles.dxcalculation,'string',dx);


ay = str2double(get(handles.ayinput,'String'));
dm = str2double(get(handles.dminput,'String'));
dy=((ay*v1n)+dm);
set (handles.dycalculation,'string',dy);


az = str2double(get(handles.azinput,'String'));
dm = str2double(get(handles.dminput,'String'));
dz=((az*v1n)+dm);
set (handles.dzcalculation,'string',dz);

bd = str2double(get(handles.bdinput,'String'));
foc = str2double(get(handles.focinput,'String'));
koc = str2double(get(handles.kocinput,'String'));
r=(1+((bd*foc*koc)/n));
set (handles.rcalculation,'string',r);

%for koc: kd and retardation factor
function varargout = kocinput_Callback(h, eventdata, handles, varargin)
koc = str2double(get(handles.kocinput,'String'));
foc = str2double(get(handles.focinput,'String'));
n = str2double(get(handles.ninput,'String'));
bd = str2double(get(handles.bdinput,'String'));
kd=foc*koc;
r=(1+((bd*foc*koc)/n));
set (handles.rcalculation,'string',r);
set (handles.kdcalculation,'string',kd);

%for foc: kd and retardation factor
function varargout = focinput_Callback(h, eventdata, handles, varargin)
koc = str2double(get(handles.kocinput,'String'));
foc = str2double(get(handles.focinput,'String'));
n = str2double(get(handles.ninput,'String'));
bd = str2double(get(handles.bdinput,'String'));
kd=foc*koc;
r=(1+((bd*foc*koc)/n));
set (handles.rcalculation,'string',r);
set (handles.kdcalculation,'string',kd);

%for bd: retardation factor
function varargout = bdinput_Callback(h, eventdata, handles, varargin)
koc = str2double(get(handles.kocinput,'String'));
foc = str2double(get(handles.focinput,'String'));
n = str2double(get(handles.ninput,'String'));
bd = str2double(get(handles.bdinput,'String'));
r=(1+((bd*foc*koc)/n));
set (handles.rcalculation,'string',r);



%for a: hydrodynamic dispersion coeff.
function varargout = axinput_Callback(h, eventdata, handles, varargin)
ax = str2double(get(handles.axinput,'String'));
dm = str2double(get(handles.dminput,'String'));
v1 = str2double(get(handles.v1input,'String'));
n = str2double(get(handles.ninput,'String'));
dx=((ax*v1/n)+dm);
set (handles.dxcalculation,'string',dx);

function varargout = ayinput_Callback(h, eventdata, handles, varargin)
ay = str2double(get(handles.ayinput,'String'));
dm = str2double(get(handles.dminput,'String'));
v1 = str2double(get(handles.v1input,'String'));
n = str2double(get(handles.ninput,'String'));
dy=((ay*v1/n)+dm);
set (handles.dycalculation,'string',dy);

function varargout = azinput_Callback(h, eventdata, handles, varargin)
az = str2double(get(handles.azinput,'String'));
dm = str2double(get(handles.dminput,'String'));
v1 = str2double(get(handles.v1input,'String'));
n = str2double(get(handles.ninput,'String'));
dz=((az*v1/n)+dm);
set (handles.dzcalculation,'string',dz);

%for dm: hydrodynamic dispersion coeff.
function varargout = dminput_Callback(h, eventdata, handles, varargin)
dm = str2double(get(handles.dminput,'String'));
v1 = str2double(get(handles.v1input,'String'));
ax = str2double(get(handles.axinput,'String'));
n = str2double(get(handles.ninput,'String'));
dx=((ax*(v1/n))+dm);
set (handles.dxcalculation,'string',dx);

v1 = str2double(get(handles.v1input,'String'));
ay = str2double(get(handles.ayinput,'String'));
n = str2double(get(handles.ninput,'String'));
dy=((ay*(v1/n))+dm);
set (handles.dycalculation,'string',dy);

v1 = str2double(get(handles.v1input,'String'));
az = str2double(get(handles.azinput,'String'));
n = str2double(get(handles.ninput,'String'));
dz=((az*(v1/n))+dm);
set (handles.dzcalculation,'string',dz);

% tinput bdinput
function varargout = tinput_Callback(h, eventdata, handles, varargin)
t = str2double(get(handles.tinput,'String'));
function varargout = teinput_Callback(h, eventdata, handles, varargin)
te = str2double(get(handles.teinput,'String'));
function varargout = dinput_Callback(h, eventdata, handles, varargin)
d = str2double(get(handles.dinput,'String'));
function varargout = tsinput_Callback(h, eventdata, handles, varargin)
ts = str2double(get(handles.tsinput,'String'));



% --------PLOT SNAPSHOT------------------------------------------------------------
function varargout = plot_Callback(h, eventdata, handles, varargin)
clear figure;
% Get user input from GUI
v1 = str2double(get(handles.v1input,'String'));
n = str2double(get(handles.ninput,'String'));
bd = str2double(get(handles.bdinput,'String'));
ax= str2double(get(handles.axinput,'String'));
ay = str2double(get(handles.ayinput,'String'));
az = str2double(get(handles.azinput,'String'));
dm = str2double(get(handles.dminput,'String'));
d = str2double(get(handles.dinput,'String'));
foc = str2double(get(handles.focinput,'String'));
koc = str2double(get(handles.kocinput,'String'));
t = str2double(get(handles.tinput,'String'));
m = str2double(get(handles.minput,'String'));


% Preliminary calculations
v1n=v1/n; dx=((ax*(v1/n))+dm); dy=((ay*(v1/n))+dm); dz=((az*(v1/n))+dm); 
kd=koc*foc; r=(1+((bd*foc*koc)/n));
v1n=v1n/r; dx=dx/r; dy=dy/r; dz=dz/r; d=d/r;

% Set Grid
[x,y,z]=meshgrid(0:10:200);
[xa,ya,za]=meshgrid(0:10:200);
x0=0;y0=100;z0=0;

% calculations
if v1==0
    x0=100; 
end

m=m*1000; %conversion from kg to g 
c=m*(exp(-((x-v1n*t).^2)/(4*dx*t)-((y).^2)/(4*dy*t)-((z).^2)/(4*dz*t)-d*t))/(8*((pi*t)^1.5)*sqrt(dx*dy*dz)); 
ccol=m*(exp(-((xa-v1n*t).^2)/(4*dx*t)-((ya).^2)/(4*dy*t)-((za).^2)/(4*dz*t)-d*t))/(8*((pi*t)^1.5)*sqrt(dx*dy*dz)); 

c=m*(exp(-((x-v1n*t-x0).^2)/(4*dx*t)-((y-y0).^2)/(4*dy*t)-((z-z0).^2)/(4*dz*t)-d*t))/(8*((pi*t)^1.5)*sqrt(dx*dy*dz));
ccol=m*(exp(-((xa-v1n*t-x0).^2)/(4*dx*t)-((ya-y0).^2)/(4*dy*t)-((za-z0).^2)/(4*dz*t)-d*t))/(8*((pi*t)^1.5)*sqrt(dx*dy*dz));  
maxc=max(max(max(c)));
set (handles.maxccalculation,'string',maxc);
        

cxy=c(:,:,1);
%Create x-y plot
axes(handles.plot1)
plot3(x(:,:,1),y(:,:,1),cxy); AXIS([0 200 0 200 0 1]); 
set(handles.plot1,'XMinorTick','on')
grid on; zlabel ('C (mg/l)'); 

ccolxy=ccol(:,:,1);
% Create color plot
axes(handles.plot2);
pcolor(xa(:,:,1),ya(:,:,1),ccolxy); CAXIS([0 0.8]);
shading interp; colormap jet; xlabel ('Distance from Source (m)'); ylabel ('Distance from Source (m)');
colorbar;drawnow;pause(2);

% Dimension of plume
lx=6*(2*dx*t)^.5;
ly=6*(2*dy*t)^.5;
lz=6*(2*dz*t)^.5;
set (handles.lxcalculation,'string',lx);
set (handles.lycalculation,'string',ly);
set (handles.lzcalculation,'string',lz);    

% ------PLAY ANIMATION-----------------------------------------------------
function varargout = animate_Callback(h, eventdata, handles, varargin)

axes(handles.plot1);clear figure;
axes(handles.plot2);clear figure;

% Get user input from GUI
v1 = str2double(get(handles.v1input,'String'));
n = str2double(get(handles.ninput,'String'));
bd = str2double(get(handles.bdinput,'String'));
ax= str2double(get(handles.axinput,'String'));
ay = str2double(get(handles.ayinput,'String'));
az = str2double(get(handles.azinput,'String'));
dm = str2double(get(handles.dminput,'String'));
d = str2double(get(handles.dinput,'String'));
foc = str2double(get(handles.focinput,'String'));
koc = str2double(get(handles.kocinput,'String'));
t = str2double(get(handles.tinput,'String'));
m = str2double(get(handles.minput,'String'));
te = str2double(get(handles.teinput,'String')); %unique to animation
ts = str2double(get(handles.tsinput,'String')); %unique to animation



% Preliminary calculations

v1n=v1/n;  dx=((ax*(v1/n))+dm); dy=((ay*(v1/n))+dm); dz=((az*(v1/n))+dm); 
kd=koc*foc; r=(1+((bd*foc*koc)/n));
v1n=v1n/r; dx=dx/r; dy=dy/r; dz=dz/r; d=d/r;


% Set Grid

[x,y,z]=meshgrid(0:10:200);
[xa,ya,za]=meshgrid(0:10:200);
x0=0;y0=100;z0=0;

% calculations
if v1==0
    x0=100; 
end


m=m*1000;
for t=0:ts:te
    
    if t==0
        c=0;
        ccol=0;
    else
        set (handles.tncalculation,'string',t);
        c=m*(exp(-((x-v1n*t-x0).^2)/(4*dx*t)-((y-y0).^2)/(4*dy*t)-((z-z0).^2)/(4*dz*t)-d*t))/(8*((pi*t)^1.5)*sqrt(dx*dy*dz));
        ccol=m*(exp(-((xa-v1n*t-x0).^2)/(4*dx*t)-((ya-y0).^2)/(4*dy*t)-((za-z0).^2)/(4*dz*t)-d*t))/(8*((pi*t)^1.5)*sqrt(dx*dy*dz));  
        maxc=max(max(max(c)));
        set (handles.maxccalculation,'string',maxc);
        
        if t==ts
            mc=max(max(max(c)));
        end
        cxy=c(:,:,1);
        %Create x-y plot
        axes(handles.plot1);
        plot3(x(:,:,1),y(:,:,1),cxy); AXIS([0 200 0 200 0 mc/2]); 
        set(handles.plot1,'layer','top');
        grid on; zlabel ('C (mg/l)'); 
                
        ccolxy=ccol(:,:,1);
        % Create color plot
        axes(handles.plot2);    
        pcolor(xa(:,:,1),ya(:,:,1),ccolxy); CAXIS([0 mc/10]);
        shading interp; colormap jet; xlabel ('Distance from Source (m)'); ylabel ('Distance from Source (m)');
        colorbar;drawnow;pause(2);
               
        % Dimension of plume
        lx=6*(2*dx*t)^.5;
        ly=6*(2*dy*t)^.5;
        lz=6*(2*dz*t)^.5;
        set (handles.lxcalculation,'string',lx);
        set (handles.lycalculation,'string',ly);
        set (handles.lzcalculation,'string',lz);        
    end
end
% -------CLOSE PROGRAM-------------------------------------------------------------
function varargout = close_Callback(h, eventdata, handles, varargin)
Close all

% -------RESTORE DEFAULTS------------------------------------------------------------
function varargout = Restore_Callback(h, eventdata, handles, varargin)
v1 = 0.5; n = 0.3; bd = 1.65; 
ax= 1; ay = .2; az = .2; dm = 0.05; d = 0;
foc = 0.005; koc = 2;
t = 50; te=150; ts=5; m=5;
v1n=1.667; 
dx=1.717; dy=1.717; dz=0.05;
lx=0; ly=0; lz=0;
r=1.055; kd=0.01;

% update fields
set(handles.v1input,'String',v1); 
set(handles.ninput,'String',n); set (handles.bdinput,'String',bd);
set(handles.axinput,'String',ax); set (handles.ayinput,'String',ay); set(handles.azinput,'String',az);
set(handles.dminput,'String',dm); set(handles.dinput,'String',d);
set(handles.focinput,'String',foc); set(handles.kocinput,'String',koc);
set(handles.tinput,'String',t); set(handles.teinput,'String',te); set(handles.tsinput,'String',ts); 
set(handles.minput,'String',m); 
set (handles.dxcalculation,'string',dx);set (handles.dycalculation,'string',dy);set (handles.dzcalculation,'string',dz);
set (handles.lxcalculation,'string',lx);set (handles.lycalculation,'string',ly);set (handles.lzcalculation,'string',lz);
set (handles.v1ncalculation,'string',v1n);
set (handles.rcalculation,'string',r);set (handles.kdcalculation,'string',kd);
set (handles.maxccalculation,'string',0);

   
 
% --- Executes on button press in plots.
function varargout = plots_Callback(h, eventdata, handles, varargin)


v1 = str2double(get(handles.v1input,'String'));
n = str2double(get(handles.ninput,'String'));
bd = str2double(get(handles.bdinput,'String'));
ax= str2double(get(handles.axinput,'String'));
ay = str2double(get(handles.ayinput,'String'));
az = str2double(get(handles.azinput,'String'));
dm = str2double(get(handles.dminput,'String'));
d = str2double(get(handles.dinput,'String'));
foc = str2double(get(handles.focinput,'String'));
koc = str2double(get(handles.kocinput,'String'));
t = str2double(get(handles.tinput,'String'));
m = str2double(get(handles.minput,'String'));
te = str2double(get(handles.teinput,'String')); %unique to animation
ts = str2double(get(handles.tsinput,'String')); %unique to animation



% Preliminary calculations

v1n=v1/n;  dx=((ax*(v1/n))+dm); dy=((ay*(v1/n))+dm); dz=((az*(v1/n))+dm); 
kd=koc*foc; r=(1+((bd*foc*koc)/n));
v1n=v1n/r; dx=dx/r; dy=dy/r; dz=dz/r; d=d/r;

g=0; % for export plot
x0=0;y0=100;z0=0;
% calculations
if v1==0
    x0=100; 
end
m=m*1000;
[x,y,z]=meshgrid(0:10:200);
for l=0:ts:te
    g=g+1 % for export plot
    if l==0
        c=0;
        ccol=0;
    else
        c=m*(exp(-((x-v1n*l-x0).^2)/(4*dx*l)-((y-y0).^2)/(4*dy*l)-((z-z0).^2)/(4*dz*l)-d*l))/(8*((pi*l)^1.5)*sqrt(dx*dy*dz));
        maxc=max(max(max(c)))
        maximumc(g)=maxc; % for export plot
        
                
        % Dimension of plume
        lxx(g)=6*(2*dx*l)^.5; % for export plot
        lyy(g)=6*(2*dy*l)^.5; % for export plot 
        lzz(g)=6*(2*dz*l)^.5; % for export plot
              
    end
end
time=0:ts:te
figure; 
subplot(2,2,1);
plot(time, maximumc);
xlabel('time (s)'); ylabel('Concentration (mg/l)'); title (' Maximum concentration vs time'); 
grid on;

subplot(2,2,2);
plot(time, lxx);
xlabel('time (s)'); ylabel('plume dimension (m)'); title (' Plume dimension in X direction vs time');
grid on;

subplot(2,2,3);
plot(time, lyy);
xlabel('time (s)'); ylabel('plume dimension (m)'); title (' Plume dimension in Y direction vs time');
grid on;
hold on;
subplot(2,2,4);
plot(time, lzz);
xlabel('time (s)'); ylabel('plume dimension (m)'); title (' Plume dimension in Z direction vs time');
grid on;




% --- Executes on button press in clear.
function varargout = clear_Callback(h, eventdata, handles, varargin)
axes(handles.plot1);plot(0,0); 
axes(handles.plot2);plot(0,0); clear all


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


