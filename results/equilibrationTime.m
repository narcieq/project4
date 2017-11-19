% input file [T <E> Cv <M> x <|M|>]
% intersted y index [2 6 3 5]
% lattice [40 60 80 100]
% T=[2:0.05:2.3]

clc;
clear;

m=4;% file number
n=4;% interested values <E> <|M|> Cv  x 
lattice=[40 60 80 100];
yindex=[2 6 3 5]; %<E> -2~0  <|M|> 0~1 Cv 0~3  x 0~137
Tstep=0.05;
%taskE_40_1000000.txt
ymax=[-1 1 3 20];
ymin=[-2 0.2 0 0];
A=['r' 'g' 'b' 'c' 'm' 'y' 'k'];

%data import
for j=1:n
    for i=1:m
    str1='taskE_';
    str2=num2str(lattice(i));
    str3='_1000000.txt';
    data=importdata([str1 str2 str3]);
    subplot(2,2,j);
    plot(data(:,1),data(:,yindex(j)),A(i));
    hold on;
    grid on;
    end
    legend('Lattice 40' ,'Lattice 60','Lattice 80','Lattice 100')
    xlabel('T/(kT/J)','FontName','Times New Roman','FontSize',12) 
    xlim([2 2.3]);
    ylim([ymin(j) ymax(j)]);  
    set(gca,'XTick',[2:Tstep:2.3])
      switch(j)
        case 1
            title('Expectation value of energy <E>');
            ylabel('<E>','FontName','Times New Roman','FontSize',12)
        case 2
            title('Expectation value of absolute magnetic moment <|M|>');
             ylabel('<|M|>','FontName','Times New Roman','FontSize',12)
        case 3
            xlim([2.23 2.3]);
            set(gca,'XTick',[2.23:0.01/2:2.3])
            title('Specific heat Cv');
            ylabel('Cv','FontName','Times New Roman','FontSize',12)
          case 4
              title('susceptibility ¦Ö ');      
              ylabel('¦Ö','FontName','Times New Roman','FontSize',12)
      end
     hold on
end
 


hold off
