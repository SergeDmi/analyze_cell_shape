function [normsP,normsF]=compute_normals(points,faces)
    % Marginal gain
	np=size(points,1);
	nf=size(faces,1);
	
	normsP=zeros(np,3);
	normsF=zeros(np,3);
	%surfF=zeros(nf,3);
	%surfP=zeros(np,3);
	cntfP=zeros(np,1);
	cntfF=zeros(nf,1);
    %A=zeros(1,3);
    %B=zeros(1,3);
    %C=zeros(1,3);
    dir=zeros(1,3);
	col3=ones(3,1);

	
    for f=1:nf
        %A(:)=points(faces(f,1),:);
        %B(:)=points(faces(f,2),:);
        %C(:)=points(faces(f,3),:);
        %dir(:)=cross((B-A),(C-A));
		dir(:)=cross((points(faces(f,2),:)-points(faces(f,1),:)),(points(faces(f,3),:)-points(faces(f,1),:)));
	    %surfP(faces(f,:),:)=surfP(faces(f,:),:)+
		normsP(faces(f,:),:)=normsP(faces(f,:),:)+col3*dir;
		cntfP(faces(f,:))=cntfP(faces(f,:))+1;
		normsF(f,:)=dir;
		cntfF(f)=cntfF(f)+1;
		
		
		% Link A to B
		
		
	end
    [normsP(:,:),surfP]=normalize_rows_3D(normsP);
	[normsF(:,:),surfF]=normalize_rows_3D(normsF);
	
	
end