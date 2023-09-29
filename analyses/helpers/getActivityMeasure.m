function act = getActivityMeasure(this_p,act_struct)

if strcmp(this_p.activityMeasure,"casc_50c_beh")
    act = act_struct.casc_50c_beh;
elseif strcmp(this_p.activityMeasure,"casc_50g_beh")
    act = act_struct.casc_50g_beh;
elseif strcmp(this_p.activityMeasure,"casc_100c_beh")
    act = act_struct.casc_100c_beh;
elseif strcmp(this_p.activityMeasure,"casc_100g_beh")
    act = act_struct.casc_100g_beh;
elseif strcmp(this_p.activityMeasure,"casc_200g_beh")
    act = act_struct.casc_200g_beh;
elseif strcmp(this_p.activityMeasure,"spks_beh")
    act = act_struct.spks_beh;
elseif strcmp(this_p.activityMeasure,"dFF_beh")
    act = act_struct.dFF_beh;
elseif strcmp(this_p.activityMeasure,"F_beh")
    act = act_struct.F_beh;
elseif strcmp(this_p.activityMeasure,"Fneu_beh")
    act = act_struct.Fneu_beh;
elseif strcmp(this_p.activityMeasure,"dFFn_beh")
    act = act_struct.dFFn_beh;
elseif strcmp(this_p.activityMeasure,"dFFnn_beh")
    act = act_struct.dFFnn_beh;
elseif strcmp(this_p.activityMeasure,"dFFns_beh")
    act = act_struct.dFFns_beh;
elseif strcmp(this_p.activityMeasure,"spksn_beh")
    act = act_struct.spksn_beh;
elseif strcmp(this_p.activityMeasure,"spksnn_beh")
    act = act_struct.spksnn_beh;
end

end