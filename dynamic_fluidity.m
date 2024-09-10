function Fluidity = dynamic_fluidity(TS)

% function: Fluidity = dynamic_fluidity(TS)
%
% This function compute dynamic fluidity as in Aguilera 2024 of a
% multidimensional time serie
%
% INPUT: TS, your time serie in observation x variables format (e.g. for
% EEG in time x channels format)
%
% OUTPUT: Fluidity, time serie with a value of dynamic fluidity per time
% point.
%
% This function is based on work published in Faranda, D.,Messori, G. & Yiou, P. Dynamical proxies of North Atlantic predictability and extremes. Sci. Rep. 7, 41278 (2017)
%
% Matthieu Aguilera, Funsy team, Université de Strasbourg, 2024

quanti=0.98;
STEP = 1;

L = length(TS);

u = 0;

Fluidity = [];
for j = 1:STEP:L

    u=u+1;

    distance=pdist2(TS(j,:),TS(setdiff(1:STEP:length(TS), j),:));
    logdista=-log(distance);

    Fluidity(u)=extremal_Sueveges(logdista,quanti);


end
end

function [theta]=extremal_Sueveges(Y,p)
u=quantile(Y, p);
q=1-p;
Li=find(Y>u);
Ti=diff(Li);
Si=Ti-1;
Nc=length(find(Si>0));
N=length(Ti);

theta=(sum(q.*Si)+N+Nc- sqrt( (sum(q.*Si) +N+Nc).^2  -8*Nc*sum(q.*Si))  )./(2*sum(q.*Si));
end