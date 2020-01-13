function [checked,experiment] = check_states(experiment, options)
% Display a pair and wait for user to decide wether to keep or not
% The user can also select a state (by typing its number)
% S. Dmitrieff  2018
% Institut Jacques Monod
% www.biophysics.fr

% @TODO : create a config file for plotting

if nargin < 1
    error('First argument should be an experiment');
end

if nargin < 2
    options=analysis_default_options();
end

pixel_size=1.0;
if isfield(options,'pixel_size')
  pixel_size=options.pixel_size;
end


scale=1.1;
% By default we don't keep
checked=0;
h=options.thickness;
%% Here we plot
n_states=numel(experiment.states);
selected_state=n_states;

if isfield(options,'add_offset')
  add_offset=options.add_offset;
else
  p1=experiment.states(1).shape.points*pixel_size;
  add_offset=(max(p1(:,2))-min(p1(:,2)))/3.0;
end

default_view=[1,0.5,0.5];
current_view=default_view;

finished=0;

while ~finished

  do_plot_of_states();





  %% Here we wait for user action

  hFig   = gcf;
  set(gcf, 'Position', [100, 100, 1000, 1000])
  set(hFig,  'KeyPressFcn',   {@callback_key_down});
  set(hFig,  'DeleteFcn',     {@callback_close});


  %% Callbacks
  if nargout > 0
      waitfor(hImage);
      finished
  end

end

%% prepare the dialog

  function do_plot_of_states()
    figure
    lasty=-add_offset;

    for s=1:n_states
      pixel_size
      points=experiment.states(s).shape.points*pixel_size;
      offset=lasty+add_offset-min(points(:,2));
      lasty=offset+(max(points(:,2))-min(points(:,2)));

      % Subplot of the top view
      subplot(4,n_states,s)
      gl1=logical((points(:,3)<h).*(points(:,3)>-h));
      scatter(points(gl1,1),points(gl1,2),5,'k');
      title('top view')
      axis equal
      axis(scale*[min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))])

      % Subplot of the side view
      subplot(4,n_states,n_states+s)
      gl1=logical((points(:,2)<h).*(points(:,2)>-h));
      scatter(points(gl1,1),points(gl1,3),5,'k');
      title('side view')
      axis equal
      axis(scale*[min(points(:,1)) max(points(:,1)) min(points(:,3)) max(points(:,3))])

      % subplot of all cells
      subplot(4,n_states,[(2*n_states+1):4*n_states]);
      hold all
      hImage=scatter3(points(:,1),points(:,2)+offset,points(:,3),1,'k');

    end

      cc=min(numel(experiment.states(1).name)-1,20);
      title(['...' experiment.states(1).name(end-cc:end)], 'Interpreter', 'none');
      axis equal;
      view(current_view)
  end


    function callback_key_down(hObject, eventData)
        %keyboard actions
        if eventData.Character == ' '
            %reset_view;
            %set_selection([], []);
	         checked=1;
	         callback_quit();
        elseif eventData.Character == 'q'
            callback_quit();
        elseif eventData.Character == 'r'
              reset_view();
        elseif eventData.Character == 's'
            select_state(input('Please input the state to select (and then select the plot window again) '))
        elseif is_numeric(eventData.Character)
          select_state(str2double(eventData.Character));
        elseif numel(strfind(eventData.Key,'arrow'))>0 || numel(strfind(eventData.Key,'page'))>0
          move_view(eventData.Key)
        else
          try_rotation(eventData.Character)
        end


	end

  function reset_view()
    current_view=default_view;
    view(current_view);
  end

    function move_view(arrow)
      intens=0.1*norm(current_view);
      switch arrow
        case 'rightarrow'
            current_view=current_view+[1,0.0,0.0]*intens;
        case 'leftarrow'
            current_view=current_view-[1,0.0,0.0]*intens;
        case 'uparrow'
            current_view=current_view+[0.0,1,0.0]*intens;
        case 'downarrow'
            current_view=current_view-[0.0,1,0.0]*intens;
        case 'pageup'
            current_view=current_view+[0.0,0.0,1]*intens;
        case 'pagedown'
            current_view=current_view-[0.0,0.0,1]*intens;
      end
      view(current_view)
    end

    function callback_quit(hObject, eventData)
        %save_with_confirmation;
        %set(hFig, 'pointer', pointer_saved);
        close_figure()
        finished=1;

    end

    function callback_close(hObject, eventData)
        %save_with_confirmation;
        %set(hFig, 'pointer', pointer_saved);
        close_figure()

    end


    function select_state(s)
      if s>0 && s<=n_states
        selected_state=s;
        disp(['Selected state : ' num2str(s)])
      else
        disp('Invalid state')
      end
    end

    function try_rotation(command)
      changed=0;
      if command=='f';
        changed=1;
        experiment.states(selected_state).shape.points=-experiment.states(selected_state).shape.points;
      end
      if changed
        close_figure();
      end
    end

    function close_figure()
      if ishandle(hFig)
          delete(hFig);
      end
    end

end

function val=is_numeric(string)
  val=~isnan(str2double(string));
end
