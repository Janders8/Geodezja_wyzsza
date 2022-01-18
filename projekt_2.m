function main
    clear;
    % Beta Capricorni
    Rektascensja = 20.35018694 % 20h 21m 00,673
    Deklinacja = -14.78138 %−14° 46′ 52,98


    % warszawa
    Phi = 52.232222 
    Lambda = 21.008333

    % Quito
    %Phi = -0.22
    %Lambda = -78.5125


    % Melbourne
    %Phi = -37.8
    %Lambda = 144.95
    
    %godziny
    h = 0:1:24;



    
    t=hourAngle(2021,11,25,h,Lambda,Rektascensja);
    t= t';

    for i=1:25 
        if t(i) > 360 
            t(i) = t(i) - 360; 
        end 
    end

    for i = 1:25
        cosZ(i) = sind(Phi)*sind(Deklinacja) + cosd(Phi)*cosd(Deklinacja) * cosd(t(i))
        up(i) = -cosd(Deklinacja)*sind(t(i)); 
        down(i) = (cosd(Phi)*sind(Deklinacja)-sind(Phi)*cosd(Deklinacja)*cosd(t(i)));

        Az(i) = atan2d(up(i),down(i)); 

    end

    cosZ = cosZ'
    Az = Az'
    


 

% pilnowanie by kąty nie były ujemne
for i=1:25
    if Az(i) < 0 
        Az(i) = Az(i) + 360; 
    end 
end 
 
Z = acosd(cosZ);   

for i = 1:25

    px(i) = sind(Z(i))*cosd(Az(i)); 
    py(i) = sind(Z(i))*sind(Az(i)); 
    pz(i) = cosd(Z(i)); 

end

% rysowanie na kopule
%[x,y,z] = sphere(50);
%z(z < 0) = 0;
%S = surf(x,y,z); 
%axis equal
%S.FaceAlpha = 0.05; 
%hold on 
%scatter3(px,py,pz, 'r', 'filled');
% rysowanie wykresu wysokości od czasu
wysokosc = 90 - Z
plot(h,wysokosc)
title("zależność wysokości od czasu - Warszawa")
xlabel("czas - h")
ylabel("wysokość słońca nad horyzontem - stopnie")




function [t] = hourAngle(y,m,d,h, lambda, rek)
    jd = juliandate(datetime(y,m,d));
    g = GMST(jd);
    g
    UT1 = h * 1.002737909350795;
    S = UT1 * 15 + lambda + g;
    S
    t = S -rek * 15;
end
function g = GMST(JD)
    T = (JD - 2451545)/36525;
    g = 280.46061837 + 360.98564736629 * (JD - 2451545)...
        + 0.000387933 * T.^2 - T.^3/38710000;
    g = mod(g,360);
end

end