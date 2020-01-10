function [M,w]=normalize_rows_ND(M)
d=size(M,2);
w=sqrt(sum(M.^2,2));
M(:,:)=M./(w*ones(1,d));
end
