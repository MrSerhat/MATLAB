% BU KOD:
%       LINE & SPHERE parametrelerine  g�re kesi�me noktalar�n� bulur,
%       grafiklerini �izer ve sonu�lar� ekrana yazar
% 
% PARAMETRLER:
%       Line: Do�ru denklemleri katsay�lar�n� i�eren vekt�r. ��erik: [x0 y0 z0  dx dy dz]
%               xo: x eksen sabiti, dx ise de�i�ken katsays�s�
%            �rnek: x = 1 -3*t i�in [xo dx] = [1 -3];
%
%       Sphere: K�re denklemi katsay�lar�n� i�eren vekt�r. ��erik: [xc yc zc  R]
%               xc, yc ve zc: k�renin x, y ve z eksenlerinden olan sapma miktar�
%                          R: K�renin yar��ap�
%               �rnek: (x-1).^2 +  (y-2).^2 + (z+3).^2  = 4 i�in: [1 2 -3 4];
%
%       Points: Kesi�im noktalar�n�n x, y ve z de�erlerini verir
%
% �RNEK KOD �ALI�TIRMASI:
%       Points = lineAndSphereIntersept([1 2 3 1 0 -2], [1 4 0 4 ])

function Points = lineAndSphereIntersept(Line, Sphere)
%function Points = lineAndSphereIntersept
clc;
% Line=[1 2 3 1 0 -2];
% Sphere=[1 4 0 4 ];

%% SAB�TLER
tol = 1e-14;

%% KES���M NOKTASININ BEL�RLENMES�

% Do�ru ve K�re Merkezlerinin k�yaslanmas�
dc = Line(1:3) - Sphere(1:3);

% Denklem katsay�lar�
a = sum(Line(:, 4:6) .* Line(:, 4:6), 2);
b = 2*sum(dc.*Line(4:6), 2);
c = sum(dc.*dc, 2) - Sphere(:,4).*Sphere(:,4);
delta = b.*b - 4*a.*c;

if delta > tol % delta positif ise 2 noktada kesi�iyor demektir.
    u1 = (-b -sqrt(delta)) / 2 / a;
    u2 = (-b +sqrt(delta)) / 2 / a;
    
    Points = [Line(1:3)+u1*Line(4:6) ; Line(1:3)+u2*Line(4:6)];

elseif abs(delta) < tol % 1 noktada kesi�iyor demektir.
    u = -b/2./a;    
    Points = Line(1:3) + u*Line(4:6);
    
else % Kesi�miyor demektir.
    Points = ones(2, 3);
    Points(:) = NaN;
end

%% GRAFIK ��Z�M�

% K�re �izimi
[X,Y,Z] = ellipsoid(Sphere(1),Sphere(2),Sphere(3),Sphere(4),Sphere(4),Sphere(4));
surf(X,Y,Z); hold on
xlabel('X')
ylabel('Y')
zlabel('Z')

% Do�ru �izimi
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

% Kesi�im Noktalar�n�n �izimi
if ~isnan(Points)
    for k=1:size(Points,1)
        plot3(Points(k,1),Points(k,2),Points(k,3),'*c','MarkerSize',14);
    end
end
hold off
legend('sphere','line','intersept')

% Sonu�lar�n ekrana yaz�lmas�
fprintf(1,'\n');
fprintf(1,'K�re Denklemi: \n');
fprintf(1,'(x-%.2f).^2 + (y-%.2f).^2  + (z-%.2f).^2 = %.3f \n',Sphere(1),Sphere(2),Sphere(3),Sphere(4));

fprintf(1,'\n');
fprintf(1,'Do�ru Denklemi: \n');
fprintf(1,'x = %.2f  %+.2f*t \n',Line(1),Line(4));
fprintf(1,'y = %.2f  %+.2f*t \n',Line(2),Line(5));
fprintf(1,'z = %.2f  %+.2f*t \n',Line(3),Line(6));
fprintf(1,'\n');

fprintf(1,'Kesi�im Noktalar�: \n');
if all(isnan(Points))
    fprintf(1,'Kesi�en Nokta Yok!');
else    
    for k=1:size(Points,1)
        fprintf(1,'%d. Nokta [x, y, z] = [%6.3f  %6.3f  %6.3f ]\n',k, Points(k,1),Points(k,2),Points(k,3));
    end
end
fprintf(1,'\n');

