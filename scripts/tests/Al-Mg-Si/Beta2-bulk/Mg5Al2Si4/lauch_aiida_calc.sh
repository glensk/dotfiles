#!/usr/bin/env sh
aiida_launch_workflow_alalloy.py \
    --calc_method "vc-relax" \
    --code_node "1" \
    --structure_group_name "NN_relaxed_Mg5Al2Si4_n2p2_v2ag" \
    --workchain_group_name "NN_relaxed_Mg5Al2Si4_n2p2_v2ag_calc" \
    --base_parameter_node "9b370584-3f56-471c-a724-dbaadf022ec5" \
    --pseudo_familyname "SSSP_v1.1_eff" \
    --kptper_recipang 80 \
    --nume2bnd_ratio 0.75 \
    --max_wallclock_seconds 21600 \
    --max_active_calculations 300 \
    --sleep_interval 600

# wallclock_secs adapted to daint (this time is fine up to 264 atoms)
# one day secs allowed on daint: 86400
# --code_node "pw-v6.3" \ not working anymore since import from daniel
