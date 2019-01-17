monthly_averages=zeros(22,100);
summ=zeros(57,22);

for(j= 1:57)
 for(i=1:22)
  indd=((j-1)*22 + i);
  x=Allout(indd,:);
  x1=reshape(x',12,100);
  monavg=sum(x1); %one basin monthly average 
  monthly_averages(i,:)=monavg;
 end 
 summ(j,:)=mean(monthly_averages');
end


boxplot(summ)
hold on 
plot(summ(1,:))
