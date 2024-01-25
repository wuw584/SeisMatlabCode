sv = zeros(1,90);
for ip = 0:90
    rf = recoff(1,0.1,1,5,3,2,sin(deg2rad(ip))/5);
    disp(rf);
    disp(ip);
    sv(ip+1) = rf(9);
end
disp(sv)
    plot(ip+1,sv, 'k:', ip, rf(12), 'k--' ,ip,rf(11) , 'k-')
    legend('P-wave trans potential', 'SV-wave re potential' , 'P-wave re potential');


    xlabel('incident angle / ^o');
    ylabel('reflectance of P-wave');

