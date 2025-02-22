function int_ts = interp_ts(ts, badInds, PLOT)
% (spline) interpolate ts around the points listed in badInds.
% catie 3/19/07


int_ts = ts;

if isempty(badInds)
  return;
end

k=1;
while (k <= length(badInds))
  
  % if first point in timeseries, 
  if badInds(k)==1
    % find next good one
    curr_ind = k;
    while (curr_ind < length(badInds))
      if badInds(curr_ind+1) == badInds(curr_ind)+1
	curr_ind= curr_ind+1;
      else
	break;
      end
    end
   
    k = curr_ind+1;
    value = ts(k);
    int_ts(1:k) = value;
 
  
  % if last point in timeseries, duplicate the prev value.
  elseif badInds(k)==length(ts)
    int_ts(badInds(k)) = ts(badInds(k)-1);
    break;
  else
  
    % if somewhere in the middle:
    xlower = badInds(k)-1;
    curr_ind = k;
    ISEND = false;
    while(curr_ind < length(badInds))
      
      if badInds(curr_ind+1) == badInds(curr_ind)+1
	curr_ind = curr_ind+1;
	% if we happen to be at the end
	if badInds(curr_ind) == length(ts)
	  ISEND = true;
	end
      else
	break;
      end
    end
    k = curr_ind+1; % start here within the vec BadInds next time
    xupper = badInds(curr_ind)+1;
    
    ylower = ts(xlower);
    if ISEND
      yupper = ylower; % if at end, duplicate values from first good.
    else
      yupper = ts(xupper);
    end
    
    dist = xupper - xlower; 
    x = [xlower xupper];
    y = [ylower, yupper];
    xi = [xlower+1 : xupper-1];
    meanval = mean(y);
    yi = ones(size(xi))*meanval;
    %yi = interp1(x,y,xi,'spline');
    int_ts(xi) = yi;
  end
end

if PLOT
  plot(ts,'g'); hold on; plot(int_ts,'r');
  legend('original','interpolated');
end
