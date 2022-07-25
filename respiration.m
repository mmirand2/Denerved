t= timeS;
A=0;
B= 0;

Tinsp= 1.5;
T=8/3*(Tinsp);


tresp= mod(t,T);


for tt=1:length(t)
if tresp(tt)< Tinsp
    x(tt)= (A/2)*( cos((2*pi*tresp(tt))/Tinsp)-1 )+B;
else
    x(tt)= B;
end

end

figure
plot(t,x)


j=0.0025 + (0.0025* x)*10;
plot(t,j)