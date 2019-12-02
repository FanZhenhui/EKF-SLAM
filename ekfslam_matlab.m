%% init params
clear;
clc;
%motion params
V = 0.3;
MAX_W = 120.0*pi/180.0;
DT_CONTROLS = 0.1;
sigmaV = 0.03;
sigmaW = 1.0*pi/180.0;
Q = [sigmaV^2 0; 0 sigmaW^2];
%Observation params
MAX_RANGE = 1.9;
sigmaR = 0.05;
sigmaB = (1.0*pi/180);
R = [sigmaR^2 0; 0 sigmaB^2];
%min distance can be considered as waypoint
MIN_WP = 0.03;

%% plot
lm = [1 2 2 3 2 0 0 3 1; 1 4 3 2 1 3 6 6 6; 1 2 3 4 5 6 7 8 9];
wp = [0 2 2 3 4 4 0 0; 0 0 0 0 7 7 7 0;];

ScreenSize = get(0, 'ScreenSize')*0.75;
fig = figure('Position', [0 0 ScreenSize(3), ScreenSize(4)]);
subplot(1,2,1)
plot(lm(1,:),lm(2,:),'k.')
hold on, axis equal, grid on
MAXX = max([max(lm(1,:)) max(wp(1,:))]);
MAXY= max([max(lm(2,:)) max(wp(2,:))]);
MINX = min([min(lm(1,:)) min(wp(1,:))]); 
MINY= min([min(lm(2,:)) min(wp(2,:))]);
axis([MINX-1 MAXX+1 MINY-1 MAXY+1])

Plot.xv = plot(0,0,'ro','erasemode','xor');
Plot.xt = plot(0,0,'bo','erasemode','xor');
Plot.dr = plot(0,0,'go','erasemode','xor');
Plot.pth_fusion = plot(0,0,'r','markersize',2,'erasemode','background');
Plot.pth_true = plot(0,0,'b','markersize',2,'erasemode','background');
Plot.pth_DRnoise = plot(0,0,'g','markersize',2,'erasemode','background');
Plot.obs = plot(0,0,'k','erasemode', 'xor');
xlabel('blue:true path | red:predicted path | green:noise path');
%xlabel(');

subplot(1,2,2)
hold on,  grid on
MINX_= 0;MAXX_ = 723; MINY_ = 0; MAXY_ = 0.3;
axis([MINX_ MAXX_ MINY_ MAXY_])
Plott =  plot(0,0,'b','markersize',2,'erasemode','background');
xlabel('error between true path and predicted path');
%% init states
global vtrue XX PX PaData
vtrue = zeros(3,1);
XX = zeros(3,1);
DR = zeros(3,1);
PX = zeros(3);
PosError = 0;
dt = DT_CONTROLS;
dtsum = 0;
iwp = 1;
W = 0;
counter = 0;
path_count = 0;
Assos = []; 
PaData.i=1;
PaData.PaFusion = XX;
PaData.Patrue = vtrue;
PaData.DrNoise = DR;
PaData.state(1).XX = XX;
PaData.state(1).PX = diag(PX)
PaData.PosError = PosError;
id_stop = 0;

%% slam
while iwp ~= 0 && id_stop == 0
    counter = counter +1;
    global vtrue;
    cwp = wp(:,iwp);
    d = sqrt((cwp(1) - vtrue(1))^2 + (cwp(2) - vtrue(2))^2);
    if d < MIN_WP
        iwp =iwp + 1;
        if iwp > size(wp,2)
            iwp=1;
        end
        cwp = wp(:,iwp);
    end

    %Angle between the two points
    W = AngleLimit(atan2(cwp(2) - vtrue(2), cwp(1)-vtrue(1)) - vtrue(3));
    if abs(MAX_W*dt) < W
        W = sign(W)*MAX_W*dt;
    end
    vtrue = MotionModel(vtrue, V, W, dt);
    %Control noise
    Vn = V + randn(1)*sqrt(Q(1,1));
    Wn = W + randn(1)*sqrt(Q(2,2));
    
    %calJacobian
    angleJ = XX(3);
    xUp = dt*Vn*cos(angleJ + Wn);
    yUp = dt*Vn*sin(angleJ + Wn);
    Gt = [1,0,xUp; 0,1,yUp; 0,0,1];
    Gu = [dt*cos(angleJ + Wn) xUp;
    dt*sin(angleJ + Wn) yUp;
    dt*sin(Wn)/0.1 Vn*dt*cos(Wn)/0.1;];

    PX(1:3, 1:3) = Gt*PX(1:3, 1:3)*Gt.' + Gu*Q*Gu.';
    if length(PX) > 3
        PX(1:3,4:end) = Gt*PX(1:3,4:end);
        PX(4:end,1:3) = PX(1:3,4:end)';
    end
    XX(1:3,:) = MotionModel(XX, Vn, Wn, dt);
    DR = MotionModel(DR, Vn, Wn, dt);
    
    %observation model
    z = []; ZId = [];
    Znew = true;
    for i=lm
        dist = sqrt((i(1) - vtrue(1))^2 + (i(2) - vtrue(2))^2);
        if dist <= MAX_RANGE
            %landmarks can be seen
            bearing = atan2((i(2) - vtrue(2)),(i(1) - vtrue(1)));

            l = [dist; AngleLimit(bearing - vtrue(3)); i(3)];
            z = [z l]
            %new land mark
            for y = 1:length(Assos)
                if i(3) == Assos(y)
                    ZId = [ZId y];
                    Znew = false;
                end
            end
            %for the new landmark, add id 
            if Znew
                Assos = [Assos i(3)];
                MuJ = [XX(1) + l(1)*cos(l(2) + XX(3)); XX(2) + l(1)*sin(l(2) + XX(3))];
                XX = [XX; MuJ(1); MuJ(2)];
                LenCovars = length(PX);
                PX(LenCovars + 1, LenCovars + 1) = 1; 
                PX(LenCovars + 2, LenCovars + 2) = 1;
            end
        end
        Znew = true;
    end

    %calc observe javobian
    j = 0;
    for zi = 1:size(z,2)
        zloc = z(:,zi);
        for r = 1:length(Assos) 
            if zloc(3) == Assos(r) 
                j = r;
            end
        end
        if (j==0)
            disp('Cannot find landmark');
        end
        jind = j*2 + 2;
        [zhat, H] = ObserveModel(XX, zloc, jind); 

        PHt = PX*H.';
        S = H*PHt + R;
        SInv = inv(S);
        K = PHt*SInv;
        diff = [zloc(1) - zhat(1); zloc(2) - zhat(2)];
        XLast = XX(1:3);
        XX = XX + K*[zloc(1) - zhat(1); zloc(2) - zhat(2)];
        L = K*S*K';
        PX = (eye(length(PX)) - K*H)*PX;
        j=0;
    end
    
    PosError = sqrt((XX(1) - vtrue(1))^2 + (XX(2) - vtrue(2))^2);
    
    %save data for plot
    global PaData
    CHUNK = 500;
    len = size(PaData.PaFusion,2);
    if (PaData.i == len)
    if len < CHUNK, len= CHUNK; end
    PaData.PaFusion = [PaData.PaFusion zeros(3,len)];
    PaData.Patrue = [PaData.Patrue zeros(3,len)];
    PaData.DrNoise = [PaData.DrNoise zeros(3,len)];
    pack
    end
    i = PaData.i + 1;
    PaData.i = i;
    PaData.PaFusion(:,i) = XX(1:3);
    PaData.Patrue(:,i) = vtrue;
    PaData.DrNoise(:,i) = DR;
    PaData.PosError(i) = PosError;
    %PaData.Posx(i) = i;
    PaData.state(i).x = XX;
    PaData.state(i).P = diag(PX);
    PaData.PosError(1:PaData.i)
    try
        set(Plot.xv, 'xdata', XX(1,:), 'ydata', XX(2,:))
        set(Plot.dr, 'xdata', DR(1,:), 'ydata', DR(2,:))
        set(Plot.pth_true, 'xdata', PaData.Patrue(1,1:PaData.i), 'ydata', PaData.Patrue(2,1:PaData.i))
        set(Plot.pth_DRnoise, 'xdata', PaData.DrNoise(1,1:PaData.i), 'ydata', PaData.DrNoise(2,1:PaData.i))
        set(Plot.pth_fusion, 'xdata', PaData.PaFusion(1,1:PaData.i), 'ydata', PaData.PaFusion(2,1:PaData.i))
        set(Plot.xt, 'xdata', vtrue(1,:), 'ydata', vtrue(2,:))
%         set(Plott, 'xdata', PaData.Patrue(1,1:PaData.i), 'ydata', PaData.Patrue(2,1:PaData.i))
        set(Plott, 'xdata', 1:PaData.i, 'ydata', PaData.PosError(1:PaData.i))
        drawnow
        %subplot(1,2,2)
        
    catch
        disp('whoops')
    end

    if(vtrue(2,:) < 0)
        id_stop = 1;
    end
end