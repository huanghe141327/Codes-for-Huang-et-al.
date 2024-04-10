function [out1,out2,out3,out4,out5,out6]=ncal_loop(ind,OE_VA,OE_VA1,u1,u2,u3,nlutbt,save_cost,save_re,save_cot,save_h,modisbt,sa,sy)
%OE_VA,OE_VA1：prior values
%u1,u2,u3: COT,CER,CTH lookup table interval
%nlutbt: look-up table for CRTM-BTS
%save_cot,save_re,save_cot,save_h
%modisbt
%sa sy
out1=[];
out2=[];
out3=[];
out4=[];
out5=[];
out6=[];
for j=1:length(ind)
    num=ind(j);
    out6=[out6,num];
    s_bt=modisbt(num,:);%MODIS bt
    try   
   
    for loop=1:300  
      %  disp(['start'])
    i_re=OE_VA(num,1);
    i_cot=OE_VA(num,2);
    i_level=OE_VA(num,3); 

    i_bt=interpn(u1,u2,u3,nlutbt,i_re(1,1),i_cot(1,1),i_level(1,1),'linear');
temp1=interpn(u1,u2,u3,nlutbt,i_re(1,1)+0.01,i_cot(1,1),i_level(1,1),'linear');
temp2=interpn(u1,u2,u3,nlutbt,i_re(1,1),i_cot(1,1)+0.01,i_level(1,1),'linear');
temp3=interpn(u1,u2,u3,nlutbt,i_re(1,1),i_cot(1,1),i_level(1,1)+0.01,'linear');


    i_btt=[i_bt(1,7:9),i_bt(1,11:16)];
    
costy= (i_btt-s_bt)*inv(sy)*(i_btt-s_bt)'...
      +(OE_VA1(num,:)-OE_VA(num,:))*inv(sa) *(OE_VA1(num,:)-OE_VA(num,:))';
    
    save_cost(num,loop)=costy;
    save_re(num,loop)=i_re;
    save_cot(num,loop)=i_cot;
    save_h(num,loop)=i_level;
%gradient cal.

i_btx1=[temp1(1,7:9),temp1(1,11:16)];
i_btx2=[temp2(1,7:9),temp2(1,11:16)];
i_btx3=[temp3(1,7:9),temp3(1,11:16)];

 costy1= (i_btx1-s_bt)*inv(sy)*(i_btx1-s_bt)'...
        +(OE_VA1(num,:)-OE_VA(num,:))*inv(sa) *(OE_VA1(num,:)-OE_VA(num,:))';
deltax1=(costy1-costy)/0.01;

 costy2= (i_btx2-s_bt)*inv(sy)*(i_btx2-s_bt)'...
        +(OE_VA1(num,:)-OE_VA(num,:))*inv(sa) *(OE_VA1(num,:)-OE_VA(num,:))';
deltax2=(costy2-costy)/0.01;

 costy3= (i_btx3-s_bt)*inv(sy)*(i_btx3-s_bt)'...
        +(OE_VA1(num,:)-OE_VA(num,:))*inv(sa) *(OE_VA1(num,:)-OE_VA(num,:))';
deltax3=(costy3-costy)/0.01;

deltaa=[deltax1,deltax2,deltax3];

%iterations
if loop<=200

deltaa=0.05*deltaa;
else

deltaa=0.01*deltaa;

end

OE_VA1(num,:)=OE_VA(num,:);
OE_VA(num,1:3)=OE_VA(num,1:3)-deltaa;
temp1=OE_VA(num,1);
temp11=OE_VA1(num,1);
temp1(temp1>89)=temp11(temp1>89)*0.9+8.9;
temp1(temp1<5)=temp11(temp1<5)*0.9+0.5;
OE_VA(num,1)=temp1;

temp2=OE_VA(num,2);
temp22=OE_VA1(num,2);
temp2(temp2>50)=temp22(temp2>50)*0.9+0.5;
temp2(temp2<0.01)=temp22(temp2<0.01)*0.9+0.001;
OE_VA(num,2)=temp2;

temp3=OE_VA(num,3);
temp33=OE_VA1(num,3);
temp3(temp3>max(u3))=temp33(temp3>max(u3))*0.9+0.1*max(u3);
temp3(temp3<0.01)=0.01;
OE_VA(num,3)=temp3;
%disp(['OK',num2str(loop)])
    end

    catch
        disp(['error',num2str(j)])
    end
out1=[out1;OE_VA(num,:)];
out2=[out2;save_cost(num,:)];
out3=[out3;save_re(num,:)];
out4=[out4;save_cot(num,:)];
out5=[out5;save_h(num,:)];
end

end
