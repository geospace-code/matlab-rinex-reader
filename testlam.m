close all;
clear all;

qahat = mkqahat2;
qzhat = decorrel(qahat);
[L,D] = ldldecom(qzhat);

figure;
axis ([-1 1 -1 1]);
hold on;

profile on;

for i = -1:0.01:1;
  for j = -1:0.01:1;
    
    afloat = [i;j];
    incr   = afloat - rem(afloat,1);
    afloat = rem(afloat,1);
    
    Chi2   = chistart (D,L,afloat,2);
    afixed = lsearch (afloat,L,D,Chi2,1);
    afixed = round(afixed + repmat(incr,1,1));

    if isequal(afixed,[0;0]);
      plot (i,j,'.');
    end;
  
  end;
end;

profile report;

points = mkpullin (qzhat);
h = plot ([points(:,1)' points(1,1)], [points(:,2)' points(1,2)],'r-');
set (h,'LineWidth',3);
