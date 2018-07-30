function val = check_pair(experiment, options)
% Display a pair and wait for user to decide wether to keep or not
% S. Dmitrieff  2018
% Institut Jacques Monod
% www.biophysics.fr

if nargin < 1
    error('First argument should be an experiment');
end

if nargin < 2
    options=pombe_default_options();
end
scale=1.1;
% By default we don't keep
val=0;
h=options.thickness;
%% Here we plot
p1=experiment.pre_pombe.points*options.pixel_size;
p2=experiment.post_pombe.points*options.pixel_size;
offset=(max(p1(:,2))-min(p2(:,2)))*1.3;
pp=[p1;p2];


figure
subplot(4,2,1)
gl1=logical((p1(:,3)<h).*(p1(:,3)>-h));
scatter(p1(gl1,1),p1(gl1,2),5,'k');
title('top view')
axis equal
axis(scale*[min(pp(:,1)) max(pp(:,1)) min(pp(:,2)) max(pp(:,2))])
subplot(4,2,2)
gl2=logical((p2(:,3)<h).*(p2(:,3)>-h));
scatter(p2(gl2,1),p2(gl2,2),5,'r');
axis equal
axis(scale*[min(pp(:,1)) max(pp(:,1)) min(pp(:,2)) max(pp(:,2))])

subplot(4,2,3)
gl1=logical((p1(:,2)<h).*(p1(:,2)>-h));
scatter(p1(gl1,1),p1(gl1,3),5,'k');
title('side view')
axis equal
axis(scale*[min(pp(:,1)) max(pp(:,1)) min(pp(:,3)) max(pp(:,3))])
subplot(4,2,4)
gl2=logical((p2(:,2)<h).*(p2(:,2)>-h));
scatter(p2(gl2,1),p2(gl2,3),5,'r');
axis equal
axis(scale*[min(pp(:,1)) max(pp(:,1)) min(pp(:,3)) max(pp(:,3))])

subplot(4,2,[5:8])
scatter3(p1(:,1),p1(:,2),p1(:,3),5,'k');
hold all
hImage=scatter3(p2(:,1),p2(:,2)+offset,p2(:,3),5,'r');
cc=min(numel(experiment(1).prename)-1,20);
title(['...' experiment(1).prename(end-cc:end)], 'Interpreter', 'none');
axis equal
view([1,-0.7,-0.8])
%rotate3d on

%% Here we wait for user action

hFig   = gcf;
set(gcf, 'Position', [100, 100, 1000, 1000])
set(hFig,  'KeyPressFcn',   {@callback_key_down});
set(hFig,  'DeleteFcn',     {@callback_quit});


%% prepare the dialog




    function callback_key_down(hObject, eventData)
        %keyboard action
        T = [];
        if eventData.Character == ' '
            %reset_view;
            %set_selection([], []);
			val=1;
			callback_quit();
        elseif eventData.Character == 'q'
            callback_quit();
        end


	end


    function callback_quit(hObject, eventData)
        %save_with_confirmation;
        %set(hFig, 'pointer', pointer_saved);
        if ishandle(hFig)
            delete(hFig);
		end

    end



%% Callbacks
if nargout > 0
    waitfor(hImage);
end




end
