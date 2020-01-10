function write_csv_with_names(filename,data,columns,rows)
% Write a csv file with name filname from data
%   filename : string ( e.g. data.csv )
%   data : matlab array ( e.g. [[1,2;3,4;5,6]])
%   columns : cell array of names, or matlab array ( e.g. {'x','y'} or [1,2] )
%   rows : cell array of names, or matlab array ( e.g. {'x','y','z'} or [1,2,3] )
% Warning : not necessarily fast.
n_cols=size(data,2);
n_rows=size(data,1);

if nargin<4
  rows=make_empty_names(n_rows);
  if nargin<3
    columns=make_empty_names(n_cols);
  end
end

if numel(columns) ~= n_cols
    columns=make_empty_names(n_cols);
end
if numel(rows) ~= n_rows
    rows=make_empty_names(n_rows);
end

rows=check_convert_input(rows);
columns=check_convert_input(columns);

% column names to line
line = columns{1};
for i = 2:n_cols
    line = [line,',',columns{i}];
end

% Create the file and write columns names
fid = fopen(filename,'w');
fprintf(fid,'%s\r\n',line);


% now write lines
for i=1:n_rows
  line=rows{i};
  for j=1:n_cols
    line = [line,',',num2str(data(i,j))];
  end
  fprintf(fid,'%s\n',line);
end

fclose(fid);

end

function list=check_convert_input(list)
% Making sure columns and rows are cell arrays
  if ~iscell(list)
    n_items=numel(list);
    new_list=cell(n_items,1);
    for i=1:n_items
      new_list{i}=num2str(list(i));
    end
    list=new_list;
  end
end

function names=make_empty_names(n)
  % Just making a n x 1 cell array of empty strings
  names=cell(n,1);
  for i=1:n
    names{i}='';
  end
end
