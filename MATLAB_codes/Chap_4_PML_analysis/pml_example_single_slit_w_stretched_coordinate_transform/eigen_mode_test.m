sorted_Peig_value=zeros(size(Peigvalue,1),size(Peigvalue,2),size(Peigvalue,3));

for k=3:-1:1
    [q1 q2]=sort(imag(Peigvalue(1,:,k)),'descend');
    sorted_Peig_value(1,:,k)=Peigvalue(1,q2,k);    
%     temp=Peigvalue(1,q2,k);    
%     sorted_Peig_value(1,:,k)=temp;
end

% [ei_1 pl_1]=get_ei_pl(Peigvalue(1,:,1),k0);
% [ei_2 pl_2]=get_ei_pl(Peigvalue(1,:,2),k0);
% [ei_3 pl_3]=get_ei_pl(Peigvalue(1,:,3),k0);
% [ei_1 pl_1]=get_ei_pl(sorted_Peig_value(1,:,1),k0);
% [ei_2 pl_2]=get_ei_pl(sorted_Peig_value(1,:,2),k0);
% [ei_3 pl_3]=get_ei_pl(sorted_Peig_value(1,:,3),k0);


figure(11);subplot(1,2,1);bar(imag(Peigvalue(1,:,1)/k0));axis([1-100 2*L+100 -1.4 1.4]);subplot(1,2,2);bar(1./(-real(Peigvalue(1,:,1))));axis([1-100 2*L+100 0 100000e-6]);
figure(12);subplot(1,2,1);bar(imag(Peigvalue(1,:,2)/k0));axis([1-100 2*L+100 -1.4 1.4]);subplot(1,2,2);bar(1./(-real(Peigvalue(1,:,2))));axis([1-100 2*L+100 0 10e-6]);
figure(13);subplot(1,2,1);bar(imag(Peigvalue(1,:,3)/k0));axis([1-100 2*L+100 -1.4 1.4]);subplot(1,2,2);bar(1./(-real(Peigvalue(1,:,3))));axis([1-100 2*L+100 0 10e-6]);
figure(21);subplot(1,2,1);bar(imag(sorted_Peig_value(1,:,1)/k0));axis([1-100 2*L+100 -1.4 1.4]);subplot(1,2,2);bar(1./(-real(sorted_Peig_value(1,:,1))));axis([1-100 2*L+100 0 100000e-6]);
figure(22);subplot(1,2,1);bar(imag(sorted_Peig_value(1,:,2)/k0));axis([1-100 2*L+100 -1.4 1.4]);subplot(1,2,2);bar(1./(-real(sorted_Peig_value(1,:,2))));axis([1-100 2*L+100 0 10e-6]);
figure(23);subplot(1,2,1);bar(imag(sorted_Peig_value(1,:,3)/k0));axis([1-100 2*L+100 -1.4 1.4]);subplot(1,2,2);bar(1./(-real(sorted_Peig_value(1,:,3))));axis([1-100 2*L+100 0 10e-6]);
% figure(11);subplot(1,2,1);bar(ei_1);axis([1-100 2*L+100 -1.4 1.4]);subplot(1,2,2);bar(pl_1);axis([1-100 2*L+100 0 10]);
% figure(12);subplot(1,2,1);bar(ei_2);axis([1-100 2*L+100 -1.4 1.4]);subplot(1,2,2);bar(pl_2);axis([1-100 2*L+100 0 10]);
% figure(13);subplot(1,2,1);bar(ei_3);axis([1-100 2*L+100 -1.4 1.4]);subplot(1,2,2);bar(pl_3);axis([1-100 2*L+100 0 10]);





%%
em=real(eprm);
ed=epra;

n_spp=real(sqrt(em*ed/(em+ed)));

fzero(@(b)sqrt(b^2-em)/em+sqrt(b^2-ed)/ed*tanh(sqrt(b^2-ed)*400*nm/2*k0),n_spp)
fzero(@(b)sqrt(b^2-em)/em+sqrt(b^2-ed)/ed*coth(sqrt(b^2-ed)*400*nm/2*k0),n_spp)

fzero(@(b)sqrt(b^2-em)/em+sqrt(b^2-ed)/ed*tanh(sqrt(b^2-ed)*100*nm/2*k0),n_spp)
fzero(@(b)sqrt(b^2-em)/em+sqrt(b^2-ed)/ed*coth(sqrt(b^2-ed)*100*nm/2*k0),n_spp)

bb=linspace(0.01,4,401);
y1=sqrt(bb.^2-em)/em+sqrt(bb.^2-ed)/ed.*tanh(sqrt(bb.^2-ed)*400*nm/2*k0);
y2=sqrt(bb.^2-em)/em+sqrt(bb.^2-ed)/ed.*coth(sqrt(bb.^2-ed)*400*nm/2*k0);
y3=sqrt(bb.^2-em)/em+sqrt(bb.^2-ed)/ed.*tanh(sqrt(bb.^2-ed)*100*nm/2*k0);
y4=sqrt(bb.^2-em)/em+sqrt(bb.^2-ed)/ed.*coth(sqrt(bb.^2-ed)*100*nm/2*k0);
figure(4);
subplot(4,1,1);plot(bb,y1);
subplot(4,1,2);plot(bb,y2);
subplot(4,1,3);plot(bb,y3);
subplot(4,1,4);plot(bb,y4);

