function SW2 = SHIN(TIME1,SW1,TIME2,lat,lon,tz,REF1,REF2)
% Interpolates the shortwave radiative flux (SW) at TIME1
% (SWhortwave INterpolation) using transmissivity. Transmissivity is
% calculated using the average cos(SZA) during the interval in TIME1 and
% interpolated to a fine (instantaneous) time step. The instantaneous
% transmissivity is then converted into a radiative flux, which is averaged
% to the desired time step (TIME2).
%
% SYNTAX:
%	SW2 = SHIN(TIME1,SW1,TIME2,lat,lon,tz,REF1,REF2)
%
% INPUTS:
%	TIME1	= Nx7 matrix - time_builder format time
%	SW		= Nx1 vector - downwelling solar surface radiative flux [Wm^-2]
%	lat		= 1x1 scalar - degrees north
%	lon		= 1x1 scalar - degrees west
%	tz		= 1x1 sclar - # of time zones West of UTC
%	REF1	= string - argument describing how the data is referenced to the time stamp
%	REF2	= string - ''
%			'END' - averaged data referenced to the end of the interval
%			'MID' - averaged data referenced to the middle of the interval
%			'BEG' - averaged data referenced to the beginning of the
%			interval
%
% OUTPUTS:
%	SW2	= Nx1 vector - SW radiative flux interpolated to TIME2 [Wm^-2]
%
% DEPENDENCIES:
%	SolarGeometry.m
%	Get_dt.m
%	TimeIndex.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHECKS/INITIAL FORMATTING %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(TIME1) ~= length(SW1)
	error('Time and shortwave vectors must be the same length')
end
if size(TIME1,2) ~=7 && size(TIME1,2) ~= 1
	error('TIME variable must either be a time_builder format matrix or a vector of serial dates')
end
if size(TIME1,2) == 7
	TIME1 = TIME1(:,7);
end
if size(TIME2,2) ~=7 && size(TIME2,2) ~= 1
	error('TIME variable must either be a time_builder format matrix or a vector of serial dates')
end
if size(TIME2,2) == 7
	TIME2 = TIME2(:,7);
end

% Time needs to be referenced to the beginning of the averaging period for
% both time series (TIME1 and TIME2)
dt1 = Get_dt(TIME1);								% Find the time step (serial format)
dt2 = Get_dt(TIME2);								% Find the time step (serial format)
if strcmp(REF1,'END')								% Referenced to end of time step
	TIME1 = TIME1 - dt1;							% Reference to the beginning of the averaging period	
elseif strcmp(REF1,'MID')							% Referenced to the middle of time step
	TIME1 = TIME1 - dt1./2;							% Referenced to the beginning of the averaging period			
elseif ~strcmp(REF1,'BEG') && ~strcmp(REF1,'MID') && ~strcmp(REF1,'END')
	error('Unrecognized REF1 option')
end
if strcmp(REF2,'END')								% Referenced to end of time step
	TIME2 = TIME2 - dt2;							% Reference to the beginning of the averaging period	
elseif strcmp(REF2,'MID')							% Referenced to the middle of time step
	TIME2 = TIME2 - dt2./2;							% Referenced to the beginning of the averaging period			
elseif ~strcmp(REF2,'BEG') && ~strcmp(REF2,'MID') && ~strcmp(REF2,'END')
	error('Unrecognized REF2 option')
end

%%%%%%%%%%
%% CODE %%
%%%%%%%%%%
% Transmissivity of the original data (TIME1)
TAU1 = Trans_AVGcosSZA(TIME1,SW1,lat,lon,tz,'BEG');	% Transmissivity @ TIME1 using average SZA

% Transmissivity @ a fine time step (t_fine) relative to TIME2
t_fine = time_builder(TIME2(1),TIME2(end)+dt2,1/12);

% Transmissivity interpolated to the fine time step
if strcmp(REF1,'END')
	TAU_fine = interp1(TIME1+dt1,TAU1,t_fine(:,7));
	% Indices outside of range of TIME1
	r1 = find(t_fine(:,7) < TIME1(1)+dt1);
	r2 = find(t_fine(:,7) > TIME1(end)+dt1);
elseif strcmp(REF1,'MID')
	TAU_fine = interp1(TIME1+dt1/2,TAU1,t_fine(:,7));
	% Indices outside of range of TIME1
	r1 = find(t_fine(:,7) < TIME1(1)+dt1/2);
	r2 = find(t_fine(:,7) > TIME1(end)+dt1/2);
end

if length(r1) > 2*max(dt1,dt2)*24*12 || length(r2) > 2*max(dt1,dt2)*24*12
	error('The range of TIME2 is greater than the range of TIME1 by more than a time step.')
end

% Fill in NaNs due to any mismatch between TIME1 and TIME2
TAU_fine(isnan(TAU_fine)) = interp1(t_fine(~isnan(TAU_fine),7),...
	TAU_fine(~isnan(TAU_fine)),t_fine(isnan(TAU_fine),7),'linear','extrap');

% Convert to a radiative flux @ a fine time step (t_fine)
EL_fine = SolarGeometry_v2(t_fine,lat,lon,tz);		% Solar geometry for fine time step
SW_fine = TAU_fine.*sind(EL_fine).*1365;			% Radiative flux @ fine time step

% Average to the outgoing time (TIME2)
TIND = TimeIndex(t_fine,dt2);						% Indices of one time step
SW2 = NaN(size(TIND,1),1);							% Pre-allocate w/ NaNs
for n = 1:size(TIND,1)
	IND = TIND(n,1):TIND(n,2);						% Interval index
	SW2(n,1) = mean(SW_fine(IND));					% Avg. SW flux over period
end
