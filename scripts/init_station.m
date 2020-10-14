function station = init_station(par)
station = containers.Map; % station monitor
station('num_occupied_pole') = 0; 
station('FLEX_list') = [];
station('ASAP_list') = [];
station('num_empty_pole') = par.station.num_poles;
station('D_init') = 0;
station('pow_cap') = par.station.num_poles * 6.6 ; % this value is arbitrary for now
station('cost_dc') = 18.86; % this value is arbitrary for now
end