a=importdata('Region_frequency.txt');
T = sum(a.data);
tot = T(1,2);
C=a.data(:,2);
% b=convertallfiletocell(a);
% c=sortrows(b,3);
% for i=1:length(c)
%     C(i,1)=(str2double(c{i,3}));
%     B{i,1}=c{i,1};
% end

%barh(C,'stacked');

set(gca,'YTickLabel',a.textdata);
figure1 = figure('XVisual','');

% Create axes
axes1 = axes('Parent',figure1,...
    'YTickLabel',a.textdata,...
    'YTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]);
box(axes1,'on');
hold(axes1,'all');

% Create bar
bar(C,'Horizontal','on','BarLayout','stacked');
