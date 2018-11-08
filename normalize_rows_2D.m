function [M,w]=normalize_rows_2D(M)
w=sqrt(sum(M.^2,2));
M(:,:)=M./(w*ones(1,2));
end
