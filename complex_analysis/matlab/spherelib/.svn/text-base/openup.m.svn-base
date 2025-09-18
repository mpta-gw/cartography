function openup(ah,xory,perc)
% OPENUP(ah,par,perc)
%
% Opens axes in which a single data set is plotted, to the 
%
% INPUT:
% 
% ah     Axis handle
% xory   1 right
%        2 top
%        3 left
%        4 bottom
%        5 left and right
%        6 top and bottom
% perc   Percentage of the data range [10]
%
% Last modified by fjsimons-at-alum.mit.edu, 11-29-2003

defval('ah',gca);
defval('xory',2);
defval('perc',10);

switch xory
 case {1,3,5}
  wat='xlim';
  wit='XData';
 case {2,4,6}
  wat='ylim';
  wit='YData';
end

for index=1:length(ah)
  %  set(ah(index),wat,get(ah(index),wat)...
  %		+[0 range(get(ah(index),wat)*perc/100)])
  % The 'indeks' is so STEM also works
  switch xory
    case {1,2}
     set(ah(index),wat,minmax(get(indeks(getkids(ah(index)),1),wit))+...
		   [0 range(get(indeks(getkids(ah(index)),1),wit))*perc/100])
    case {3,4}
     set(ah(index),wat,minmax(get(indeks(getkids(ah(index)),1),wit))+...
		   [-range(get(indeks(getkids(ah(index)),1),wit))*perc/100 0])
   case {5,6}
     oldr=minmax(get(indeks(getkids(ah(index)),1),wit));
     if oldr(1)==oldr(2)
       oldr(1)=oldr(1)-oldr(1)/2;
       oldr(2)=oldr(2)+oldr(2)/2;
       if oldr(1)==0
	 oldr=[-1 1];
       end
     end
     set(ah(index),wat,oldr+...
		   [-range(get(indeks(getkids(ah(index)),1),wit))*perc/100 ...
		    range(get(indeks(getkids(ah(index)),1),wit))*perc/100])
  end
end

