function [M,w]=normalize_rows_ND(M)
w=sqrt(sum(M.^2,2));
M(:,:)=M./(w*ones(1,3));
end