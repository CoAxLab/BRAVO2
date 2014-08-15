function Data=dload(filename)
% DLOAD: loads a column-oriented ascii-data file into memory
% synopsis
%		Data=dload(filename)
% Description
%		the file filename has to be tab or space delimited all-numeric or character datafile
%		first row has to be a header file with valid variable names.
%       if there is an output argument, dload gives it to a struture with
%       field names as the columns
%       Otherwise it assigns them to the workspace.
%       If transforming to numericals is unsuccessful, it leaves them as
%       strings
fid=fopen(filename,'r');
if (fid==-1)
    fprintf(['Error: Could not find ' filename '\n']);
    Header=[];Data=[];
    return;    
end;
Header=fgetl(fid);
H={};
Data=[];
%eval(['global ' Header]);
Head=Header;
while length(Head)>0
    [r,Head]=strtok(Head);
    if (~isempty(r))
        H={H{1:end} r};
    end;
end;
%Header=H;
try 
    A=textread(filename,'%f','headerlines',1);    
catch   
    A=textread(filename,'%s','headerlines',1);    
end;
Indx=[1:length(H):length(A)]';
% assign into variable names
for col=1:length(H)
    if (max(Indx)>length(A))
        warning(['File format does not have rows*columns entries (e.g. empty cells):' filename]);
        Data=[];
        return;
    end;
    if (iscell(A))
        V={A{Indx}}';
        d=str2num(char(V));
        if (~isempty(d))
            V=d;
        end;
    else
        V=A(Indx);
    end;
    name=strrep(H{col},'-','_');
    if (nargout==0)        % no output variable: just assign to global scope as single variables
        assignin('caller',name,V);
    else 
        Data=setfield(Data,name,V); 
    end;
    Indx=Indx+1;
end;
fclose(fid);


