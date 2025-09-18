function [ h ] = smartfigure( hin )
% smartfigure is a wrapper routine for to MATLAB function figure
% is has several calling syntaxes:
%
% A) call with 1 argument hin
%    if hin is a string
%    --> this simply sets the environmental variable 'SMARTFIGURE' to hin
%    if hin is a number
%    --> this sets the environmental variable 'SMARTFIGURE' to num2str(hin)
% B) call with no argument
%    the behavior depend on the value of the environmental variable
%    'SMARTFIGURE':
%
%      Value  | Behavior 
% ------------------------------------------------------------------------
%       ''    | identical to 'NEXT' - behaviour just as figure
%             |
%      'NEXT' | h=figure; is called and h returned
%             | i.e. the hande of a new figure is returned
%             |
%      'CURR' | h=gcf; is called and h returned
%             | i.e. the hande of the current figure is returned
%             |
%      'SAME' | h=figure(gcf); is called and h returned
%             | i.e. the hande of the current figure is returned
%             |
%      'num'  | where num is an integer number:
%             | h=figure(num); is called, h is returned and
%             | num2str(num+1) is stored in the env. var. 'SMARTFIGURE'
%             |
% 'n ORDER l' | where n is a number, and l is a list of numbers
%             | h=figure(n); is called, h is returned, 
%             | n is found in the list l, m = next list entry,
%             | or m = first list entry if
%             |   - n is not in the list or
%             |   - n is the last list entry.
%             | 'm ORDER l' is stored in the env. var. 'SMARTFIGURE'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             | Additional feature:
% 'Position [xmin,ymin,xwidth,ywidth]':
%             | sets the environmental variable 'SMARTFIGUREPOSITION'
%             | which then is checked on a no-argument call, and the new
%             | figure is placed to the position [xmin,ymin,xwidth,ywidth]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin>0
        if isstr(hin)
            pind=strfind(upper(hin),'POSITION');
            if pind
                setenv('SMARTFIGUREPOSITION',hin(pind+8:end));
            else
                setenv('SMARTFIGURE',hin);
            end
        else
            setenv('SMARTFIGURE',num2str(hin));
        end
    else
        val =getenv('SMARTFIGURE');
        uval=upper(val);
        oind=strfind(uval,'ORDER');
        if strcmp(val,'')
            hout=figure;
        elseif strcmp(uval,'NEXT')
            hout=figure;
        elseif strcmp(val,'CURR')
            hout=gcf;
        elseif strcmp(uval,'SAME')
            hout=figure(gcf);
        elseif oind
            valstr =uval(1:oind-1);
            v      =str2num(valstr); %#ok<*ST2NM>
            liststr=uval(oind+5:end);
            list   =str2num(liststr);
            hout=figure(v);
            mind=find(list==v, 1);
            if isempty(mind)
                m=list(1);
            else
                try, %#ok<NOCOM>
                    m=list(mind+1);
                catch
                    m=list(1);
                end
            end
            valout=[num2str(m),' ORDER ',liststr];
            setenv('SMARTFIGURE',valout);
        else
            v=str2num(val);
            hout=figure(v);
            setenv('SMARTFIGURE',num2str(v+1));
        end
        valpos =getenv('SMARTFIGUREPOSITION');
        if not(strcmp(valpos,''))
            set(hout,'Position',str2num(valpos));
        end
    end
    if nargout>0
        h=hout;
    end
end

