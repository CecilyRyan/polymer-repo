function [Sigs, sigsd,sigsp,Sigm,d,p]= OwensWendt(A1, A2, Temp, name,RT)

% Owens Wendt Theory for Surface Enegry Calculation 
% Using water and Diiodomethane
% change liquids as necessary


TM= Temp; % Melt temp concerned about mixing at
RT= RT;   %Room Temperature
% Insert Liquid data here

sigl=[72.8   % water
      50.8]; %DIMETH
             %known surface tension values for liquid at room temp

% Break into 2 component model
% Liquid on PTFE.. Polar component= 0 Dispersive= 18 mJ/m^2
%thetaPTFE= [113.7
           % 100.7 ]; %DIIMETH
                     % Liquid contact angle on PTFE


% Computing the Dispersive and Polar components of each
sigld=[21.8  % water%(sigl.^2.*(cosd(thetaPTFE)+1).^2)/72;
       50.4]; % DIMETH 
             % Dispersive component
siglp=[51.0
         .4] ;%sigl-sigld Polar component


    

% Insert contact angles for Polymer in question
theta1=(A1 ) ;   %water
theta2=(A2 ) ; % glycerol






theta= [mean(theta1)
        mean(theta2)];

Y=zeros(1,length(theta));
X=zeros(1,length(theta));

for i=1:length(theta)
     Y(i)= sigl(i)*(cosd(theta(i))+1)/(2*sqrt(sigld(i)));
     X(i)= sqrt(siglp(i))/sqrt(sigld(i));
end

P = polyfit(X,Y,1); % x = x data, y = y data, 1 = order of the polynomial.
fit = polyval(P,X);

plot(X,Y,'o',X,fit,'-')
legend('data','Linear Fit')

% using the slope and the slope intercept to fing the surface tension
sigsp=(P(1))^2; % slope     ie: polar
sigsd=(P(2))^2; % Intercept ie: Dispersive

Sigs= sigsp+sigsd;
fprintf('Surface Tension for ' )
fprintf(name)
fprintf('\n\n')
fprintf('          AT ROOM TEMPERATURE\n')
fprintf('Overall surface Energy = %3.3f mj/m^2  \n', Sigs);
fprintf('Dipsersive component   = %3.3f  mj/m^2 \n',sigsd);
fprintf('Polar component        = %3.3f  mj/m^2 \n\n',sigsp);


% Ratio Method of determining the Melt polar and Dispersive components

dgdt=.06;% accepted value in Ohmega paper % (11/9)* Go/Tc*(1-(20/Tc))^(2/9);
Sigm= Sigs-(dgdt*(TM-RT));
d=(sigsd/Sigs)*Sigm;
p=(sigsp/Sigs)*Sigm;


fprintf('          At %3.0f degrees Celcius\n',TM);
fprintf('Overall surface Energy = %3.3f mj/m^2  \n', Sigm);
fprintf('Dipsersive component   = %3.3f  mj/m^2 \n',d);
fprintf('Polar component        = %3.3f  mj/m^2 \n\n',p);
end
