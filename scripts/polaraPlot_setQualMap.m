%
% setting legend for quality
%
k=abs(kvar); sg=sign(kvar);
hold on;
j=1;h1(1)=plot(c{i}(1,1),sg*c{i}(1,k),conr{j},'MarkerFaceColor',conr{j}(2:2));
j=2;h1(2)=plot(c{i}(1,1),sg*c{i}(1,k),conr{j},'MarkerFaceColor',conr{j}(2:2));
j=4;h1(3)=plot(c{i}(1,1),sg*c{i}(1,k),conr{j},'MarkerFaceColor',conr{j}(2:2));
j=6;h1(4)=plot(c{i}(1,1),sg*c{i}(1,k),conr{j},'MarkerFaceColor',conr{j}(2:2));
j=8;h1(5)=plot(c{i}(1,1),sg*c{i}(1,k),conr{j},'MarkerFaceColor',conr{j}(2:2));
j=12;h1(6)=plot(c{i}(1,1),sg*c{i}(1,k),conr{j},'MarkerFaceColor',conr{j}(2:2));


for j=1:nR
   plot(c{i}(j,1),sg*c{i}(j,k),conr{cq(j)+2},'MarkerFaceColor',...
       conr{cq(j)+2}(2:2),'MarkerSize',10)
end


%
qst{1}='CQ: N/A';
qst{2}='CQ: not conveged';
qst{3}='CQ: Coefficients convegence (X100)';
qst{4}='CQ: Coefficients convegence (X10)';
qst{5}='CQ: Coefficients convegence (X1)';
qst{6}='CQ: Residue convegence';
hlg=legend(h1,qst(:),'FontSize',10);
hold off;

