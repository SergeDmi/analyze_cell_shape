function [experiments]=load_exps(root,folders)
% We scan over many folders (experiments)
% for which we get several sub-experiment
% for which there is pre and post
% for which there are several cells
% ...
% Let's go !

%% Ok now we have to declare the root and folders
%root='/home/dmitrief/Sync/ownCloud/work/yeast_wall/pombe_3D_files/';
%folders={'180214/','180220/','180222/'};
nfolders=numel(folders);
count=0;
clear experiments;
c_pre=0;
for N=1:nfolders
	fold=[root,folders{N}];
	files=dir(fold);
	% Get a logical vector that tells which is a directory.
	dirFlags=[files.isdir];
	folds=files(dirFlags);
	nf=numel(folds);
	% We go over all GFP_* folders
	for n=1:nf
		%disp(folds(n).name);
		foldname=folds(n).name;
		if ~strcmp(foldname(1),'.')
			fname=[fold,foldname];
			prename=[fname,'/pre/'];
			postname=[fname,'/post/'];
			if exist(prename,'file')==7
				exps=dir(prename);
				expFlags=[exps.isdir];
				exps=exps(expFlags);
				ne=numel(exps);
				% Now we go over all cells in pre/post
				for e=1:ne					
					expname_pre =[prename, exps(e).name,'/'];
					fname_pre =dir([expname_pre,'*.ply']);
					if numel(fname_pre)>0
						c_pre=c_pre+1;
						plyname=fname_pre(1).name;
						posts=dir(postname);
						pFlags=[posts.isdir];
						posts=posts(pFlags);
						np=numel(posts);
						for p=1:np
							expname_post=[postname,posts(p).name,'/'];
							fname_post =dir([expname_post,'*.ply']);
							if numel(fname_post)>0
								if strcmp(fname_post(1).name(2:5),plyname(2:5))
									count=count+1;
									experiments(count).plyname=plyname;
									experiments(count).prename =[expname_pre,plyname];
									experiments(count).postname=[expname_post,fname_post(1).name];
								end
							end
						end
						%plyname_pre=[expname_pre,fname_pre(1).name];
					end
				end
			end
		end
	end
end

Nexp=count;
%% Now importing ply files into coordinates, and doing the analysis already
for n=1:Nexp
	experiments(n).pre_pombe=import_ply(experiments(n).prename);
	experiments(n).post_pombe=import_ply(experiments(n).postname);
    % we shouldn't do analysis here !
	%experiments(n).pre_analysis=analyze_pombe(experiments(n).pre_pombe);
	%experiments(n).post_analysis=analyze_pombe(experiments(n).post_pombe);
	%disp([num2str(n) '/' num2str(Nexp)]);
end

%res=summary_experiments(experiments);
end