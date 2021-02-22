% BU KOD:
%       LINE & SPHERE parametrelerine  göre kesiþme noktalarýný bulur,
%       grafiklerini çizer ve sonuçlarý ekrana yazar
% 
% PARAMETRLER:
%       Line: Doðru denklemleri katsayýlarýný içeren vektör. Ýçerik: [x0 y0 z0  dx dy dz]
%               xo: x eksen sabiti, dx ise deðiþken katsaysýsý
%            Örnek: x = 1 -3*t için [xo dx] = [1 -3];
%
%       Sphere: Küre denklemi katsayýlarýný içeren vektör. Ýçerik: [xc yc zc  R]
%               xc, yc ve zc: kürenin x, y ve z eksenlerinden olan sapma miktarý
%                          R: Kürenin yarýçapý
%               Örnek: (x-1).^2 +  (y-2).^2 + (z+3).^2  = 4 için: [1 2 -3 4];
%
%       Points: Kesiþim noktalarýnýn x, y ve z deðerlerini verir
%
% ÖRNEK KOD ÇALIÞTIRMASI:
%       Points = lineAndSphereIntersept([1 2 3 1 0 -2], [1 4 0 4 ])

function Points = lineAndSphereIntersept(Line, Sphere)
%function Points = lineAndSphereIntersept
clc;
% Line=[1 2 3 1 0 -2];
% Sphere=[1 4 0 4 ];

%% SABÝTLER
tol = 1e-14;

%% KESÝÞÝM NOKTASININ BELÝRLENMESÝ

% Doðru ve Küre Merkezlerinin kýyaslanmasý
dc = Line(1:3) - Sphere(1:3);

% Denklem katsayýlarý
a = sum(Line(:, 4:6) .* Line(:, 4:6), 2);
b = 2*sum(dc.*Line(4:6), 2);
c = sum(dc.*dc, 2) - Sphere(:,4).*Sphere(:,4);
delta = b.*b - 4*a.*c;

if delta > tol % delta positif ise 2 noktada kesiþiyor demektir.
    u1 = (-b -sqrt(delta)) / 2 / a;
    u2 = (-b +sqrt(delta)) / 2 / a;
    
    Points = [Line(1:3)+u1*Line(4:6) ; Line(1:3)+u2*Line(4:6)];

elseif abs(delta) < tol % 1 noktada kesiþiyor demektir.
    u = -b/2./a;    
    Points = Line(1:3) + u*Line(4:6);
    
else % Kesiþmiyor demektir.
    Points = ones(2, 3);
    Points(:) = NaN;
end

%% GRAFIK ÇÝZÝMÝ

% Küre çizimi
[X,Y,Z] = ellipsoid(Sphere(1),Sphere(2),Sphere(3),Sphere(4),Sphere(4),Sphere(4));
surf(X,Y,Z); hold on
xlabel('X')
ylabel('Y')
zlabel('Z')

% Doðru Çizimi
L1=-50;
L2=50;
t=L1:0.1:L2;
XLine=Line(1) + Line(4).*t;
YLine=Line(2) + Line(5).*t;
ZLine=Line(3) + Line(6).*t;
plot3(XLine,YLine,ZLine,'-r','LineWidth',2); 
clear L1 L2

xlim([min(X(:))-1 max(X(:))+1])
ylim([min(Y(:))-1 max(Y(:))+1])
zlim([min(Z(:))-1 max(Z(:))+1])

% Kesiþim Noktalarýnýn Çizimi
if ~isnan(Points)
    for k=1:size(Points,1)
        plot3(Points(k,1),Points(k,2),Points(k,3),'*c','MarkerSize',14);
    end
end
hold off
legend('sphere','line','intersept')

% Sonuçlarýn ekrana yazýlmasý
fprintf(1,'\n');
fprintf(1,'Küre Denklemi: \n');
fprintf(1,'(x-%.2f).^2 + (y-%.2f).^2  + (z-%.2f).^2 = %.3f \n',Sphere(1),Sphere(2),Sphere(3),Sphere(4));

fprintf(1,'\n');
fprintf(1,'Doðru Denklemi: \n');
fprintf(1,'x = %.2f  %+.2f*t \n',Line(1),Line(4));
fprintf(1,'y = %.2f  %+.2f*t \n',Line(2),Line(5));
fprintf(1,'z = %.2f  %+.2f*t \n',Line(3),Line(6));
fprintf(1,'\n');

fprintf(1,'Kesiþim Noktalarý: \n');
if all(isnan(Points))
    fprintf(1,'Kesiþen Nokta Yok!');
else    
    for k=1:size(Points,1)
        fprintf(1,'%d. Nokta [x, y, z] = [%6.3f  %6.3f  %6.3f ]\n',k, Points(k,1),Points(k,2),Points(k,3));
    end
end
fprintf(1,'\n');

