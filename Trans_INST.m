function TAU = Trans_INST(SW,EL)
% Calculates transmittance of downwelling shortwave irradiance using the
% instantaneous SZA at the time stamp.
%
% SYNTAX:
%	TAU = Trans_INST(SW,EL)
% 
% INPUTS:
%	SW	= Nx1 vector - Downwelling solar surface radiative flux [Wm^-2]
%	EL	= Nx1 vector - elevation angle [degrees]
%
% OUTPUTS:
%	TAU = Nx1 vector - Transmissivity of the flux [0-1]

%%%%%%%%%%%%
%% CHECKS %%
%%%%%%%%%%%%
if length(EL) ~= length(SW)
	error('Time and elevation angle vectors must be the same length')
end

%%%%%%%%%%
%% CODE %%
%%%%%%%%%%
TOA = 1365.*sind(EL);						% Quick and dirty TOA
TAU = SW ./ TOA;							% Transmittance
TAU(isnan(TAU)) = 0;						% Correct for physicality
TAU(TAU > 1) = 1;							% Correct for physicality
