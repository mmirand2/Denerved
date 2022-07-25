function Elh = ElastanceBasic(t,T,Ts,Tr,Ed,Es)

%will have to change like the rat function
% global tpars ppars
for i = 1:length(t)
    if t(i)<=Ts
        Elh(i) = Ed + ((Es-Ed)/2).*(1-cos(pi.*t(i)/Ts));
    elseif (t(i)<=Ts+Tr) && (t(i)>Ts)
        Elh(i) = Ed + ((Es-Ed)./2)*(cos((pi./Tr).*(t(i)-Ts)) + 1);
    else
        Elh(i)=Ed;
    end
end