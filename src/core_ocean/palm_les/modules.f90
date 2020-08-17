!> @file modules.f90
!------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1997-2018 Leibniz Universitaet Hannover
!------------------------------------------------------------------------------!
!
! Current revisions:
! ------------------
!
! 2018-11-15 cbegeman
! Add bubble property variables
!
! Former revisions:
! -----------------
!
! 2018-10-25 cbegeman
! Add dirichlet bottom boundary conditions for salinity
!
! $Id: modules.f90 3083 2018-06-19 14:03:12Z gronemeier $
! set dt_3d = 0.01
!
! 3065 2018-06-12 07:03:02Z Giersch
! Variables concerning stretching introduced or revised
!
! 3045 2018-05-28 07:55:41Z Giersch
! z_max_do2d removed
!
! 3026 2018-05-22 10:30:53Z schwenkel
! Changed the name specific humidity to mixing ratio, since we are computing
! mixing ratios.
!
! 3014 2018-05-09 08:42:38Z maronga
! Added default values of u_max, v_max, and w_max to avoid floating invalid
! during spinup
!
! 3004 2018-04-27 12:33:25Z Giersch
! precipitation_rate removed
!
! 3003 2018-04-23 10:22:58Z Giersch
! The inversion height is defined as a global variable now which belongs to the
! module statistics
!
! 2968 2018-04-13 11:52:24Z suehring
! +topo_min_level
!
! 2964 2018-04-12 16:04:03Z raasch
! *_time_count variables are all initialized with zero now
!
! 2918 2018-03-21 15:52:14Z gronemeier
! -l_grid, -l_wall
!
! 2906 2018-03-19 08:56:40Z Giersch
! Module control_parameters has been extended with ENVIRONMENT variables
! read/write_svf
!
! 2894 2018-03-15 09:17:58Z Giersch
! _prerun flags were removed, Control paramters restart_string and length have
! been added
!
! 2881 2018-03-13 16:24:40Z suehring
! Added flag for switching on/off calculation of soil moisture
!
! 2797 2018-02-08 13:24:35Z suehring
! +ghf_av
!
! 2776 2018-01-31 10:44:42Z Giersch
! Variable synthetic_turbulence_generator has been abbreviated and _prerun flags
! for skipping module related restart data has beed introduced
!
! 2765 2018-01-22 11:34:58Z maronga
! Set initial value for time_since_reference_point
!
! 2746 2018-01-15 12:06:04Z suehring
! +plant_canopy
!
! 2742 2018-01-12 14:59:47Z suehring
! +tsurf_av
!
! 2735 2018-01-11 12:01:27Z suehring
! +r_a_av
!
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Implementation of uv exposure model (FK)
! + turbulence closure variables (control_parameters)
! + arrays for prognostic equation of disspiation (arrays_3d)
! + km_av, kh_av (TG)
! Implementation of chemistry module (FK)
! -lod
! +topo_distinct (MS)
!
! 2669 2017-12-06 16:03:27Z raasch
! CONTIGUOUS-attribut added to 3d pointer arrays,
! coupling_char extended to LEN=8
!
! 2575 2017-10-24 09:57:58Z maronga
! Renamed phi -> latitude, moved longitude from radiation model to modules
!
! 2563 2017-10-19 15:36:10Z Giersch
! Variable wind_turbine was added to control_parameters
!
! 2550 2017-10-16 17:12:01Z boeske
! complex_terrain namelist parameter added
!
! 2508 2017-10-02 08:57:09Z suehring
! Change default value for pt/q/s/sa_vertical_gradient_level
!
! 2499 2017-09-22 16:47:58Z kanani
! Default changed to fft_method = 'temperton-algorithm'
!
! 2408 2017-09-05 15:47:53Z gronemeier
! Changed default value of mg_cycles from -1 to 4.
!
! 2375 2017-08-29 14:10:28Z schwenkel
! Moved mass_of_solute, molecular_weight_of_solute, molecular_weight_of_water,
! vanthoff back from particle attributes because they can now also be used in
! bulk microphysics.
! Added aerosol_bulk, aerosol_nacl, aerosol_c3h4o4, aerosol_nh4no3
!
! 2372 2017-08-25 12:37:32Z sward
! y_shift namelist parameter added
!
! 2339 2017-08-07 13:55:26Z gronemeier
! corrected timestamp in header
!
! 2338 2017-08-07 12:15:38Z gronemeier
! moved 1d-model varaibles to own module model_1d_mod
!
! 2337 2017-08-07 08:59:53Z gronemeier
! -old_dt_1d
! +l1d_diss
!
! 2326 2017-08-01 07:23:24Z gronemeier
! Updated variable descriptions
!
! 2320 2017-07-21 12:47:43Z suehring
! -ptnudge, qnudge, tnudge, td_lsa_lpt, td_lsa_q, td_sub_lpt, td_sub_q, ug_vert,
!  vg_vert, unudge, vnudge, wsubs_vert, shf_surf, p_surf, pt_surf, q_surt,
!  qsws_surf, tmp_tnudge, timenudge, time_surf, time_vert
!
! 2300 2017-06-29 13:31:14Z raasch
! default value for host changed to '????', default value for loop_optimization
! changed to 'cache', default value for termination_time_needed set to 35.0
!
! 2298 2017-06-29 09:28:18Z raasch
! missing variable descriptions have been added,
! type of write_binary changed from CHARACTER to LOGICAL
! -plot_precision, plot_3d_precision, return_addres, return_username,
! avs_data_file, exchange_mg, sendrecvcound_yxd, sendrecv_in_background,
! port_name, profile_number, cross_ts_numbers, cross_ts_number_count,
! dots_crossindex, dots_index, cross_ts_uymax, cross_ts_uymax_computed,
! cross_ts_uymin, cross_ts_uymin_computed
!
! 2296 2017-06-28 07:53:56Z maronga
! Added parameters for model spinup
!
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison'
! includes two more prognostic equations for cloud drop concentration (nc)
! and cloud water content (qc).
!
! 2277 2017-06-12 10:47:51Z kanani
! Added doxygen comments for variables/parameters,
! removed unused variables dissipation_control, do2d_xy_n, do2d_xz_n, do2d_yz_n,
! do3d_avs_n, lptnudge, lqnudge, lunudge, lvnudge, lwnudge, skip_do_avs,
! sums_up_fraction_l.
!
! 2259 2017-06-08 09:09:11Z gronemeier
! Implemented synthetic turbulence generator
!
! 2256 2017-06-07 13:58:08Z suehring
! Change default value of zeta_min to -20
! Increase dimension for wall_heatflux, etc.
!
! 2233 2017-05-30 18:08:54Z suehring
!
! 2232 2017-05-30 17:47:52Z suehring
! Renamed wall_flags_0 and wall_flags_00 into advc_flags_1 and advc_flags_2,
! respectively. Moreover, introduced further flag array wall_flags_0.
!
! Adjustments for new topography concept:
!   -fwxm, fwxp, fwym, fwyp, fxm, fxp, fym, fyp, rif_wall, wall_e_x, wall_e_y,
!   -wall_v, wall_u, wall_w_x, wall_w_y, wall_qflux, wall_sflux, wall_nrflux,
!   -wall_qrflux
!
! Adjustments for new surface concept:
!   +land_surface
!   -z0, z0h, z0q, us, ts, qs, qsws, nrs, nrsws, qrs, qrsws, ssws, ss, saswsb
!   -nzb_diff_u, nzb_diff_v, nzt_diff
!   -uswst, vswst, tswst, sswst, saswst, qswst, qrswst, nrswst, qswst_remote
!
! Generic tunnel setup:
!   +tunnel_height, tunnel_length, tunnel_width_x, tunnel_width_y,
!   +tunnel_wall_depth
!
! Topography input via netcdf
!   +lod
!
! 2200 2017-04-11 11:37:51Z suehring
! -monotonic_adjustment
!
! 2174 2017-03-13 08:18:57Z maronga
! Changed default values for most_method to 'newton'
!
! 2118 2017-01-17 16:38:49Z raasch
! -acc_rank, background_communication, i_left, i_right, j_south, j_north,
!  num_acc_per_node, on_device
!
! 2107 2017-01-09 12:21:49Z kanani
! Preparation for doxygen comments (Giersch)
!
! 2050 2016-11-08 15:00:55Z gronemeier
! Implement turbulent outflow condition
!
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean and rho_av to rho_ocean_av
!
! 2011 2016-09-19 17:29:57Z kanani
! +urban_surface, +lsf_exception, +varnamelength
!
! 2007 2016-08-24 15:47:17Z kanani
! Increased DIMENSION of data_output, data_output_user, do2d, do3d
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1992 2016-08-12 15:14:59Z suehring
! +constant_top_scalarflux, top_scalarflux
! default of bc_s_t adjusted
!
! 1968 2016-07-18 12:01:49Z suehring
! Changed dimension for MPI-datatypes type_x_int and type_y_int
!
! 1960 2016-07-12 16:34:24Z suehring
! Separate humidity and passive scalar
! +bc_s_t_val, diss_s_s, diss_l_s, flux_s_s, flux_l_s, s, sp, s1, s2, s3, ssws_av,
!  s_init, s_surf, sums_wsss_ws_l, ss, ssws, sswst, ts_m, wall_sflux
! +constant_scalarflux, ibc_s_b, ibc_s_t, s_vertical_gradient_level_ind
!
! Unused variables removed
! -gamma_x, gamma_y, gamma_z, var_x, var_y, var_z
!
! Change initial values (in order to unify gradient calculation):
! pt_vertical_gradient_level, sa_vertical_gradient_level
!
! 1957 2016-07-07 10:43:48Z suehring
! +fl_max, num_leg, num_var_fl, num_var_fl_user, var_fl_max, virtual_flight
!
! 1918 2016-05-27 14:35:57Z raasch
! default timestep switched from -1.0 to +1.0 in order to avoid wrong sign of
! initially calculated divergence
!
! 1906 2016-05-24 14:38:08Z suehring
! default value of mg_switch_to_pe0_level changed to -1
!
! 1849 2016-04-08 11:33:18Z hoffmann
! bfactor, mass_of_solute, molecular_weight_of_solute, molecular_weight_of_water,
! vanthoff moved to mod_particle_attributes.
! dt_micro and several cloud_parameters moved to microphysics_mod.
! 1d-microphysics profiles moved to microphysics_mod.
!
! 1845 2016-04-08 08:29:13Z raasch
! -nzb_2d
!
! 1833 2016-04-07 14:23:03Z raasch
! spectra parameter moved to spectra module
!
! 1831 2016-04-07 13:15:51Z hoffmann
! curvature_solution_effects removed
! turbulence renamed collision_turbulence, drizzle renamed
! cloud_water_sedimentation
!
! 1822 2016-04-07 07:49:42Z hoffmann
! icloud_scheme removed. microphysics_sat_adjust, microphysics_kessler,
! microphysics_seifert added.
!
! 1815 2016-04-06 13:49:59Z raasch
! cpp-directive for decalpha removed
!
! 1808 2016-04-05 19:44:00Z raasch
! MPI module used by default on all machines
!
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
!
! 1788 2016-03-10 11:01:04Z maronga
! Added roughness length for moisture (z0q)
!
! 1786 2016-03-08 05:49:27Z raasch
! module spectrum moved to new separate module
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf variables moved to the netcdf-interface module
!
! 1779 2016-03-03 08:01:28Z raasch
! coupling_char extended to LEN=3
!
! 1764 2016-02-28 12:45:19Z raasch
! some reformatting
!
! 1762 2016-02-25 12:31:13Z hellstea
! +nest_* variables, size of volume_flow arrays increased by one element
!
! 1738 2015-12-18 13:56:05Z raasch
! +mean_surface_level_height
!
! 1695 2015-10-27 10:03:11Z maronga
! Removed rif (forgotten in last revision)
!
! 1693 2015-10-27 08:35:45Z maronga
! Renamed zp -> z_mo
!
! 1691 2015-10-26 16:17:44Z maronga
! Renamed Obukhov length. Added ol, removed rif. Increased number of profiles
! (pr_palm). Changed default values for dissipation_1d to 'detering' and
! (mixing_length_1d to 'blackadar'. Added most_method. rif_min and rif_max
! renamed to zeta_min and zeta_max and new values assigned.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1677 2015-10-02 13:25:23Z boeske
! +ngp_yz_int, type_xz_int, type_yz_int
!
! 1666 2015-09-23 07:31:10Z raasch
! +user_interface_current_revision, user_interface_required_revision in
! control_parameters
!
! 1639 2015-08-31 14:46:48Z knoop
! Bugfix: string 'unknown' extended to match LEN=13
!
! 1575 2015-03-27 09:56:27Z raasch
! +ngp_xz
!
! 1560 2015-03-06 10:48:54Z keck
! +recycling_yshift
!
! 1557 2015-03-05 16:43:04Z suehring
! +monotonic_adjustment
!
! 1551 2015-03-03 14:18:16Z maronga
! Increased pr_palm to 120. Increased length of dots_unit and dots_label to 13
! digits. Increased length of domask, do2d, and do3d to 20 digits.
!
! 1496 2014-12-02 17:25:50Z maronga
! Renamed "radiation" -> "cloud_top_radiation"
!
! 1484 2014-10-21 10:53:05Z kanani
! Changes due to new module structure of the plant canopy model:
!   canopy-model related parameters/variables moved to module
!   plant_canopy_model_mod
!
! 1468 2014-09-24 14:06:57Z maronga
! Adapted for use on up to 6-digit processor cores.
! Increased identifier string length for user-defined quantities to 20.
!
! 1450 2014-08-21 07:31:51Z heinze
! ntnudge from 100 to 1000 increased to allow longer simulations
!
! 1431 2014-07-15 14:47:17Z suehring
! +var_d
!
! 1429 2014-07-15 12:53:45Z knoop
! +ensemble_member_nr to prepare the random_generator for ensemble runs
!
! 1382 2014-04-30 12:15:41Z boeske
! Renamed variables which store large scale forcing tendencies
! pt_lsa -> td_lsa_lpt, pt_subs -> td_sub_lpt,
! q_lsa  -> td_lsa_q,   q_subs  -> td_sub_q
!
! 1365 2014-04-22 15:03:56Z boeske
! Usage of large scale forcing enabled:
! increase pr_palm from 90 to 100 to allow for more standard profiles
! + ngp_sums_ls, pt_lsa, pt_subs, q_lsa, q_subs, tmp_tnudge, sums_ls_l,
! use_subsidence_tendencies
!
! 1361 2014-04-16 15:17:48Z hoffmann
! tend_* removed
! call_microphysics_at_all_substeps added
! default of drizzle set to true
!
! 1359 2014-04-11 17:15:14Z hoffmann
! particle_attributes moved to mod_particle_attributes.f90
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1327 2014-03-21 11:00:16Z raasch
! REAL constants defined as wp-kind
! -avs_output, data_output_format, do3d_compress, iso2d_output, netcdf_output
!
! 1320 2014-03-20 08:40:49Z raasch
! ONLY-attribute added to USE-statements,
! kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! old module precision_kind is removed,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1318 2014-03-17 13:35:16Z raasch
! module cpulog moved to new separate module-file
! interface for cpu_log removed
!
! 1314 2014-03-14 18:25:17Z suehring
! + log_z_z0, number_of_sublayers, z0_av_global
! 1308 2014-03-13 14:58:42Z fricke
! +ntdim_2d_xy, ntdim_2d_xz, ntdim_2d_yz, ntdim_3d
!
! 1257 2013-11-08 15:18:40Z raasch
! set default values for grid indices of maximum velocity components
! u|v|w_max_ijk
!
! 1241 2013-10-30 11:36:58Z heinze
! Usage of nudging enabled
! +nudging, ntnudge, ptnudge, qnudge, tnudge, unudge, vnudge, wnudge
! increase pr_palm from 80 to 90 to allow for more standard profiles
!
! Enable prescribed time depenend surface fluxes and geostrophic wind read in
! from external file LSF_DATA
! +large_scale_forcing, lsf_surf, lsf_vert, nlsf, time_surf, shf_surf, qsws_surf,
!  pt_surf, q_surf, p_surf, time_vert, ug_vert, vg_vert, wsubs_vert
!
! 1221 2013-09-10 08:59:13Z raasch
! wall_flags_0 changed to 32bit int, +wall_flags_00,
! +rflags_s_inner, rflags_invers
!
! 1216 2013-08-26 09:31:42Z raasch
! +transpose_compute_overlap,
! several variables are now defined in the serial (non-parallel) case also
!
! 1212 2013-08-15 08:46:27Z raasch
! +tri
!
! 1179 2013-06-14 05:57:58Z raasch
! +reference_state, ref_state, use_initial_profile_as_reference, vpt_reference,
! use_reference renamed use_single_reference_value
!
! 1159 2013-05-21 11:58:22Z fricke
! -bc_lr_dirneu, bc_lr_neudir, bc_ns_dirneu, bc_ns_neudir
! +use_cmax
!
! 1128 2013-04-12 06:19:32Z raasch
! +background_communication, i_left, i_right, j_north, j_south, req, req_count,
! send_receive, sendrecv_in_background, wait_stat
!
! 1115 2013-03-26 18:16:16Z hoffmann
! unused variables removed
!
! 1113 2013-03-10 02:48:14Z raasch
! +on_device
!
! 1111 2013-03-08 23:54:10Z raasch
! +tric, nr_timesteps_this_run
!
! 1106 2013-03-04 05:31:38Z raasch
! array_kind renamed precision_kind, pdims defined in serial code
! bugfix: default value assigned to coupling_start_time
!
! 1095 2013-02-03 02:21:01Z raasch
! FORTRAN error in r1092 removed
!
! 1092 2013-02-02 11:24:22Z raasch
! character length in some derived type changed for better alignment
!
! 1065 2012-11-22 17:42:36Z hoffmann
! + c_sedimentation, limiter_sedimentation, turbulence, a_1, a_2, a_3, b_1, b_2,
! + b_3, c_1, c_2, c_3, beta_cc
!
! bottom boundary condition of qr, nr changed from Dirichlet to Neumann
!
! 1053 2012-11-13 17:11:03Z hoffmann
! necessary expansions according to the two new prognostic equations (nr, qr)
! of the two-moment cloud physics scheme:
! +*_init, flux_l_*, diss_l_*, flux_s_*, diss_s_*, *sws, *swst, tend_*, *, *_p
! +t*_m, *_1, *_2, *_3, *_av, bc_*_b, bc_*_t, ibc_*_b, ibc_*_t, bc_*_t_val,
! +*_vertical_gradient, *_surface_initial_change, *_vertical_gradient_level,
! +*_vertical_gradient_level_ind, *_surface, constant_waterflux_*,
! +cloud_scheme, icloud_scheme, surface_waterflux_*, sums_ws*s_ws_l, wall_*flux
!
! constants for the two-moment scheme:
! +a_vent, a_term, b_vent, b_term, c_evap, c_term, cof, eps_sb, k_cc, k_cr, k_rr,
! +k_br, kappa_rr, kin_vis_air, mu_constant_value, nc, pirho_l, dpirho_l, rho_1,
! +schmidt, schmidt_p_1d3, stp, x0, xmin, xmax, dt_precipitation, w_precipitation
!
! steering parameters for the two_moment scheme:
! +mu_constant, ventilation_effect
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1031 2012-10-19 14:35:30Z raasch
! +output_format_netcdf
!
! 1015 2012-09-27 09:23:24Z raasch
! +acc_rank, num_acc_per_node,
! -adjust_mixing_length
!
! 1010 2012-09-20 07:59:54Z raasch
! pointer free version can be generated with cpp switch __nopointer
!
! 1003 2012-09-14 14:35:53Z raasch
! -grid_matching, nxa, nya, etc., nnx_pe, nny_pe, spl_*
!
! 1001 2012-09-13 14:08:46Z raasch
! -asselin_filter_factor, cut_spline_overshoot, dt_changed, last_dt_change,
! last_dt_change_1d, long_filter_factor, overshoot_limit_*, ups_limit_*
! several pointer/target arrays converted to normal ones
!
! 996 2012-09-07 10:41:47Z raasch
! -use_prior_plot1d_parameters
!
! 978 2012-08-09 08:28:32Z fricke
! +c_u_m, c_u_m_l, c_v_m, c_v_m_l, c_w_m, c_w_m_l,
! +bc_lr_dirneu, bc_lr_neudir, bc_ns_dirneu, bc_ns_neudir
! -km_damp_x, km_damp_y, km_damp_max, outflow_damping_width
! +z0h, z0h_av, z0h_factor, z0h1d
! +ptdf_x, ptdf_y, pt_damping_width, pt_damping_factor
!
! 964 2012-07-26 09:14:24Z raasch
! -cross_linecolors, cross_linestyles, cross_normalized_x, cross_normx_factor,
! cross_normalized_y, cross_normy_factor, cross_pnc_local,
! cross_profile_numbers, cross_profile_number_counter, cross_uxmax,
! cross_uxmax_computed, cross_uxmax_normalized,
! cross_uxmax_normalized_computed, cross_uxmin, cross_uxmin_computed,
! cross_uxmin_normalized, cross_uxmin_normalized_computed, cross_uymax,
! cross_uymin, cross_xtext, dopr_crossindex, dopr_label, linecolors, linestyles,
! nz_do1d, profil_output, z_max_do1d, z_max_do1d_normalized
!
! 951 2012-07-19 14:22:52Z hoffmann
! changing profile_columns and profile_rows
!
! 940 2012-07-09 14:31:00Z raasch
! +neutral
!
! 927 2012-06-06 19:15:04Z raasch
! +masking_method
!
! 880 2012-04-13 06:28:59Z raasch
! gathered_size, subdomain_size moved to control_parameters
!
! 866 2012-03-28 06:44:41Z raasch
! interface for global_min_max changed
!
! 861 2012-03-26 14:18:34Z suehring
! +wall_flags_0
! -boundary_flags
! +nzb_max
! +adv_sca_1, +adv_mom_1
!
! 849 2012-03-15 10:35:09Z raasch
! +deleted_particles, deleted_tails, tr.._count_sum, tr.._count_recv_sum in
! particle_attributes,
! +de_dx, de_dy, de_dz in arrays_3d,
! first_call_advec_particles renamed first_call_lpm
!
! 828 2012-02-21 12:00:36Z raasch
! +dissipation_classes, radius_classes, use_kernel_tables,
! particle feature color renamed class
!
! 825 2012-02-19 03:03:44Z raasch
! +bfactor, curvature_solution_effects, eps_ros, molecular_weight_of_water,
! vanthoff, -b_cond in cloud_parameters,
! dopts_num increased to 29, particle attributes speed_x|y|z_sgs renamed
! rvar1|2|3
! wang_collision_kernel and turbulence_effects_on_collision replaced by
! collision_kernel, hall_kernel, palm_kernel, wang_kernel
!
! 809 2012-01-30 13:32:58Z marongas
! Bugfix: replaced .AND. and .NOT. with && and ! in the preprocessor directives
!
! 807 2012-01-25 11:53:51Z maronga
! New cpp directive "__check" implemented which is used by check_namelist_files.
! New parameter check_restart has been defined which is needed by
! check_namelist_files only.
!
! 805 2012-01-17 15:53:28Z franke
! Bugfix collective_wait must be out of parallel branch for runs in serial mode
!
! 801 2012-01-10 17:30:36Z suehring
! Dimesion of sums_wsus_ws_l, ! sums_wsvs_ws_l, sums_us2_ws_l, sums_vs2_ws_l,
! sums_ws2_ws_l, sums_wspts_ws_l, sums_wsqs_ws_l, sums_wssas_ws_l increased.
! for thread-safe summation in advec_ws.
!
! RCS Log replace by Id keyword, revision history cleaned up
!
! Revision 1.95  2007/02/11 13:18:30  raasch
! version 3.1b (last under RCS control)
!
! Revision 1.1  1997/07/24 11:21:26  raasch
! Initial revision
!
!
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of global variables
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of all arrays defined on the computational grid.
!------------------------------------------------------------------------------!
 MODULE arrays_3d

    USE kinds

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ddzu                   !< 1/dzu
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ddzu_pres              !< modified ddzu for pressure solver
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  dd2zu                  !< 1/(dzu(k)+dzu(k+1))
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  dzu                    !< vertical grid size (u-grid)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ddzw                   !< 1/dzw
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  dzw                    !< vertical grid size (w-grid)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  hyp                    !< hydrostatic pressure
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_init                !< initial profile of potential temperature
                                                                   !< (or total water content with active cloud physics)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rdf                    !< rayleigh damping factor for velocity components
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  rdf_sc                 !< rayleigh damping factor for scalar quantities
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  sa_init                !< initial profile of salinity (ocean)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  u_init                 !< initial profile of horizontal velocity component u
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  v_init                 !< initial profile of horizontal velocity component v
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  u_stk                  !< Stokes dirft in x-direction
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  v_stk                  !< Stokes dirft in y-direction
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  u_stk_zw               !< Stokes dirft in x-direction, at w-levels
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  v_stk_zw               !< Stokes dirft in y-direction, at w-levels
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  zu                     !< vertical grid coordinate of u-grid (in m)
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  zw                     !< vertical grid coordinate of w-grid (in m)
#ifdef __SP
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  u_LS_forcing           !< large scale forcing for u
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  v_LS_forcing           !< large scale forcing for v
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  pt_LS_forcing          !< large scale forcing for pt
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  sa_LS_forcing          !< large scale forcing for sa
#endif

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ddzuw                 !< 1/(dzu*dzw)

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  d           !< divergence
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  kh          !< eddy diffusivity for heat
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  km          !< eddy diffusivity for momentum
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  sgs_diss    !< sub-grid-scale dissipation
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tend        !< tendency field (time integration)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  tric        !< coefficients of the tridiagonal matrix for solution of the Poisson equation in Fourier space

#if defined( __nopointer )
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  e          !< subgrid-scale turbulence kinetic energy (sgs tke)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  e_p        !< prognostic value of sgs tke
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  p          !< perturbation pressure
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  prho       !< potential density
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  pt         !< potential temperature
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  pt_p       !< prognostic value of potential temperature
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  rho_ocean  !< density of ocean
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  alpha_T    !< drhodT
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  beta_S     !< drhodS
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  solar3d    !< 3d solar tendency
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  sa         !< ocean salinity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  sa_p       !< prognostic value of ocean salinity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  tdiss_m    !< weighted tendency of diss for previous sub-timestep (Runge-Kutta)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  te_m       !< weighted tendency of e for previous sub-timestep (Runge-Kutta)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  tpt_m      !< weighted tendency of pt for previous sub-timestep (Runge-Kutta)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  tsa_m      !< weighted tendency of sa for previous sub-timestep (Runge-Kutta)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  tu_m       !< weighted tendency of u for previous sub-timestep (Runge-Kutta)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  tv_m       !< weighted tendency of v for previous sub-timestep (Runge-Kutta)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  tw_m       !< weighted tendency of w for previous sub-timestep (Runge-Kutta)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  u          !< horizontal velocity component u (x-direction)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  u_p        !< prognostic value of u
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  v          !< horizontal velocity component v (y-direction)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  v_p        !< prognostic value of v
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  w          !< vertical velocity component w (z-direction)
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  w_p        !< prognostic value of w
#else
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  e_1     !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  e_2     !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  e_3     !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  p       !< pointer: perturbation pressure
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  prho_1  !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  pt_1    !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  pt_2    !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  pt_3    !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  rho_1   !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  solar3d_1
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  alpha_T_1
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  beta_S_1
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  sa_1    !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  sa_2    !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  sa_3    !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  u_1     !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  u_2     !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  u_3     !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  v_1     !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  v_2     !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  v_3     !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  w_1     !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  w_2     !< pointer for swapping of timelevels for respective quantity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  w_3     !< pointer for swapping of timelevels for respective quantity

    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  e          !< pointer: subgrid-scale turbulence kinetic energy (sgs tke)
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  e_p        !< pointer: prognostic value of sgs tke
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  prho       !< pointer: potential density
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  pt         !< pointer: potential temperature
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  pt_p       !< pointer: prognostic value of potential temperature
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  rho_ocean  !< pointer: density of ocean
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  alpha_T    !< pointer: thermal expansion coefficent
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  beta_S     !< pointer: haline contraction coefficient
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  solar3d    !< pointer: 3d solar tendency
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  sa         !< pointer: ocean salinity
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  sa_p       !< pointer: prognostic value of ocean salinity
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  te_m       !< pointer: weighted tendency of e for previous sub-timestep (Runge-Kutta)
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  tpt_m      !< pointer: weighted tendency of pt for previous sub-timestep (Runge-Kutta)
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  tsa_m      !< pointer: weighted tendency of sa for previous sub-timestep (Runge-Kutta)
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  tu_m       !< pointer: weighted tendency of u for previous sub-timestep (Runge-Kutta)
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  tv_m       !< pointer: weighted tendency of v for previous sub-timestep (Runge-Kutta)
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  tw_m       !< pointer: weighted tendency of w for previous sub-timestep (Runge-Kutta)
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  u          !< pointer: horizontal velocity component u (x-direction)
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  u_p        !< pointer: prognostic value of u
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  v          !< pointer: horizontal velocity component v (y-direction)
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  v_p        !< pointer: prognostic value of v
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  w          !< pointer: vertical velocity component w (z-direction)
    REAL(wp), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::  w_p        !< pointer: prognostic value of w
#endif

    REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::  tri    !<  array to hold the tridiagonal matrix for solution of the Poisson equation in Fourier space (4th dimension for threads)

!    SAVE

 END MODULE arrays_3d


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of variables needed for time-averaging of 2d/3d data.
!------------------------------------------------------------------------------!
 MODULE averaging

    USE kinds

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  ol_av                  !< avg. Obukhov length
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  shf_av                 !< avg. surface heat flux
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  shf_sol_av             !< avg. surface solar flux

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  e_av          !< avg. subgrid-scale tke
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  kh_av         !< avg. eddy diffusivity for heat
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  km_av         !< avg. eddy diffusivity for momentum
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  p_av          !< avg. perturbation pressure
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  pt_av         !< avg. potential temperature
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  rho_ocean_av  !< avg. ocean density
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  alpha_T_av       !< avg. thermal expansion coefficient
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  beta_S_av        !< avg. haline contraction coefficient
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  solar3d_av    !< avg. 3d solar tendency
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  sa_av         !< avg. salinity
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  u_av          !< avg. horizontal velocity component u
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  v_av          !< avg. horizontal velocity component v
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, TARGET ::  w_av          !< avg. vertical velocity component

 END MODULE averaging


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of general constants.
!------------------------------------------------------------------------------!
 MODULE constants

    USE kinds

    REAL(wp) ::  pi = 3.141592654_wp  !< PI
    REAL(wp) ::  adv_mom_1            !< 1/4 - constant used in 5th-order advection scheme for momentum advection (1st-order part)
    REAL(wp) ::  adv_mom_3            !< 1/24 - constant used in 5th-order advection scheme for momentum advection (3rd-order part)
    REAL(wp) ::  adv_mom_5            !< 1/120 - constant used in 5th-order advection scheme for momentum advection (5th-order part)
    REAL(wp) ::  adv_sca_1            !< 1/2 - constant used in 5th-order advection scheme for scalar advection (1st-order part)
    REAL(wp) ::  adv_sca_3            !< 1/12 - constant used in 5th-order advection scheme for scalar advection (3rd-order part)
    REAL(wp) ::  adv_sca_5            !< 1/60 - constant used in 5th-order advection scheme for scalar advection (5th-order part)

!    SAVE

 END MODULE constants


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of parameters for program control
!------------------------------------------------------------------------------!
 MODULE control_parameters

    USE kinds

    TYPE file_status
       LOGICAL ::  opened         !< file is currently open
       LOGICAL ::  opened_before  !< file is currently closed, but has been openend before
    END TYPE file_status

    INTEGER(iwp), PARAMETER ::  varnamelength = 10        !< length of output variable names

    TYPE(file_status), DIMENSION(200) ::                &  !< indicates if file is open or if it has been opened before
                             openfile = file_status(.FALSE.,.FALSE.)

    CHARACTER (LEN=1)    ::  timestep_reason = ' '                        !< 'A'dvection or 'D'iffusion criterion, written to RUN_CONTROL file
                                                                          !< '_NV': vertically nested atmosphere PE, '_N##': PE of nested domain ##
    CHARACTER (LEN=8)    ::  most_method = 'newton'                       !< namelist parameter
    CHARACTER (LEN=8)    ::  run_date                                     !< date of simulation run, printed to HEADER file
    CHARACTER (LEN=8)    ::  run_time                                     !< time of simulation run, printed to HEADER file
    CHARACTER (LEN=9)    ::  simulated_time_chr                           !< simulated time, printed to RUN_CONTROL file
    CHARACTER (LEN=12)   ::  version = ' '                                !< PALM version number
    CHARACTER (LEN=12)   ::  revision = ' '                               !< PALM revision number
    CHARACTER (LEN=12)   ::  user_interface_current_revision = ' '        !< revision number of the currently used user-interface (must match user_interface_required_revision)
    CHARACTER (LEN=12)   ::  user_interface_required_revision = ' '       !< required user-interface revision number
    CHARACTER (LEN=20)   ::  bc_e_b = 'neumann'                           !< namelist parameter
    CHARACTER (LEN=20)   ::  bc_p_b = 'neumann'                           !< namelist parameter
    CHARACTER (LEN=20)   ::  bc_p_t = 'neumann'                         !< namelist parameter
    CHARACTER (LEN=20)   ::  bc_pt_b = 'neumann'                        !< namelist parameter
    CHARACTER (LEN=20)   ::  bc_pt_t = 'neumann'                 !< namelist parameter
    CHARACTER (LEN=20)   ::  bc_sa_t = 'neumann'                          !< namelist parameter
    CHARACTER (LEN=20)   ::  bc_sa_b = 'neumann'                          !< namelist parameter
    CHARACTER (LEN=20)   ::  bc_uv_b = 'neumann'                        !< namelist parameter
    CHARACTER (LEN=20)   ::  bc_uv_t = 'neumann'                        !< namelist parameter
    CHARACTER (LEN=20)   ::  random_generator = 'random-parallel'         !< namelist parameter
    CHARACTER (LEN=20)   ::  turbulence_closure = 'Moeng_Wyngaard'        !< namelist parameter
    CHARACTER (LEN=64)   ::  host = '????'                                !< hostname on which PALM is running, ENVPAR namelist parameter provided by mrun
    CHARACTER (LEN=80)   ::  log_message                                  !< user-defined message for debugging (sse data_log.f90)
    CHARACTER (LEN=80)   ::  run_identifier                               !< run identifier as given by mrun option -d, ENVPAR namelist parameter provided by mrun
    CHARACTER (LEN=100)  ::  initializing_actions = 'set_constant_profiles'                   !< namelist parameter
    CHARACTER (LEN=100)  ::  restart_string = ' '                         !< for storing strings in case of writing/reading restart data
    CHARACTER (LEN=210)  ::  run_description_header                       !< string containing diverse run informations as run identifier, coupling mode, host, ensemble number, run date and time
    CHARACTER (LEN=1000) ::  message_string = ' '                         !< dynamic string for error message output

    CHARACTER (LEN=varnamelength), DIMENSION(500) ::  data_output = ' '
    CHARACTER (LEN=varnamelength), DIMENSION(500) ::  data_output_user = ' '  !< namelist parameter
    CHARACTER (LEN=varnamelength), DIMENSION(500) ::  doav = ' '              !< label array for multi-dimensional,
                                                                              !< averaged output quantities

    CHARACTER (LEN=varnamelength), DIMENSION(300) ::  data_output_pr = ' '

    CHARACTER (LEN=varnamelength), DIMENSION(200) ::  data_output_pr_user = ' '  !< namelist parameter

    CHARACTER (LEN=varnamelength), DIMENSION(0:1,500) ::  do3d = ' '  !< label array for 3d output quantities

    INTEGER(iwp) ::  abort_mode = 1                    !< abort condition (nested runs)
    INTEGER(iwp) ::  average_count_pr = 0              !< number of samples in vertical-profile output
    INTEGER(iwp) ::  average_count_3d = 0              !< number of samples in 3d output
    INTEGER(iwp) ::  current_timestep_number = 0       !< current timestep number, printed to RUN_CONTROL file
    INTEGER(iwp) ::  disturbance_level_ind_b           !< lowest grid index where flow disturbance is applied
    INTEGER(iwp) ::  disturbance_level_ind_t           !< highest grid index where flow disturbance is applied
    INTEGER(iwp) ::  doav_n = 0                        !< number of 2d/3d output quantities subject to time averaging
    INTEGER(iwp) ::  dopr_n = 0                        !< number of profile output quantities subject to time averaging
    INTEGER(iwp) ::  dopr_time_count = 0               !< number of output intervals for profile output
    INTEGER(iwp) ::  dp_level_ind_b = 0                !< lowest grid index for external pressure gradient forcing
    INTEGER(iwp) ::  ensemble_member_nr = 0            !< namelist parameter
    INTEGER(iwp) ::  grid_level                        !< current grid level handled in the multigrid solver
    INTEGER(iwp) ::  ibc_e_b                           !< integer flag for bc_e_b
    INTEGER(iwp) ::  ibc_p_b                           !< integer flag for bc_p_b
    INTEGER(iwp) ::  ibc_p_t                           !< integer flag for bc_p_t
    INTEGER(iwp) ::  ibc_pt_b                          !< integer flag for bc_pt_b
    INTEGER(iwp) ::  ibc_pt_t                          !< integer flag for bc_pt_t
    INTEGER(iwp) ::  ibc_sa_t                          !< integer flag for bc_sa_t
    INTEGER(iwp) ::  ibc_sa_b                          !< integer flag for bc_sa_b
    INTEGER(iwp) ::  ibc_uv_b                          !< integer flag for bc_uv_b
    INTEGER(iwp) ::  ibc_uv_t                          !< integer flag for bc_uv_t
    INTEGER(iwp) ::  intermediate_timestep_count       !< number of current Runge-Kutta substep
    INTEGER(iwp) ::  intermediate_timestep_count_max   !< maximum number of Runge-Kutta substeps
    INTEGER(iwp) ::  io_group = 0                      !< I/O group to which the PE belongs (= #PE / io_blocks)
    INTEGER(iwp) ::  io_blocks = 1                     !< number of blocks for which I/O is done in sequence (total number of PEs / maximum_parallel_io_streams)
    INTEGER(iwp) ::  iran = -1234567                   !< integer random number used for flow disturbances
    INTEGER(iwp) ::  length = 0                        !< integer that specifies the length of a string in case of writing/reading restart data
    INTEGER(iwp) ::  maximum_grid_level                !< number of grid levels that the multigrid solver is using
    INTEGER(iwp) ::  maximum_parallel_io_streams = -1  !< maximum number of parallel io streams that the underlying parallel file system allows, set with mrun option -w, ENVPAR namelist parameter, provided by mrun
    INTEGER(iwp) ::  mid                               !< masked output running index
    INTEGER(iwp) ::  nr_timesteps_this_run = 0         !< number of timesteps (cpu time measurements)
    INTEGER(iwp) ::  number_stretch_level_start        !< number of user-specified start levels for stretching
    INTEGER(iwp) ::  number_stretch_level_end          !< number of user-specified end levels for stretching
    INTEGER(iwp) ::  nz_do3d = -9999                   !< namelist parameter
    INTEGER(iwp) ::  runnr = 0                         !< number of run in job chain
    INTEGER(iwp) ::  timestep_count = 0                !< number of timesteps carried out since the beginning of the initial run

    INTEGER(iwp) ::  dist_nxl                               !< left boundary of disturbance region
    INTEGER(iwp) ::  dist_nxr                               !< right boundary of disturbance region
    INTEGER(iwp) ::  dist_nyn                               !< north boundary of disturbance region
    INTEGER(iwp) ::  dist_nys                               !< south boundary of disturbance region
    INTEGER(iwp) ::  do3d_no(0:1) = 0                            !< number of 3d output quantities
    INTEGER(iwp) ::  do3d_time_count(0:1) = 0                    !< number of output intervals for 3d data
    INTEGER(iwp) ::  dz_stretch_level_end_index(9)               !< vertical grid level index until which the vertical grid spacing is stretched
    INTEGER(iwp) ::  dz_stretch_level_start_index(9)             !< vertical grid level index above which the vertical grid spacing is stretched
    INTEGER(iwp) ::  pt_vertical_gradient_level_ind(10) = -9999  !< grid index values of pt_vertical_gradient_level(s)
    INTEGER(iwp) ::  sa_vertical_gradient_level_ind(10) = -9999  !< grid index values of sa_vertical_gradient_level(s)
    INTEGER(iwp) ::  stokes_drift_method = -9999                 !< method to compute Stokes drift
                                                                 !<  1: exponential profile
                                                                 !<  2: use empirical wave spectrum of Donelan et al., 1985
#ifdef __SP
    INTEGER(wp) ::  mode_LS = 0                         !< large scale forcing mode
#endif

    INTEGER(iwp), DIMENSION(0:1) ::  ntdim_3d     !< number of output intervals for 3d data

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  grid_level_count  !< internal switch for steering the multigrid v- and w-cycles

    LOGICAL ::  poisfft_initialized = .FALSE.
    LOGICAL ::  call_psolver_at_all_substeps = .TRUE.            !< namelist parameter
    LOGICAL ::  constant_heatflux = .TRUE.                       !< heat flux at all surfaces constant?
    LOGICAL ::  constant_top_heatflux = .TRUE.                   !< heat flux at domain top constant?
    LOGICAL ::  constant_top_salinityflux = .TRUE.               !< salinity flux at ocean domain top?
    LOGICAL ::  create_disturbances = .TRUE.                     !< namelist parameter
    LOGICAL ::  data_output_2d_on_each_pe = .TRUE.               !< namelist parameter
    LOGICAL ::  disturbance_created = .FALSE.                    !< flow disturbance imposed?
    LOGICAL ::  do3d_at_begin = .FALSE.                          !< namelist parameter
    LOGICAL ::  do_sum = .FALSE.                                 !< contribute to time average of profile data?
    LOGICAL ::  dp_external = .FALSE.                            !< namelist parameter
    LOGICAL ::  dp_smooth = .FALSE.                              !< namelist parameter
    LOGICAL ::  dt_fixed = .FALSE.                               !< fixed timestep (namelist parameter dt set)?
    LOGICAL ::  dt_3d_reached                                    !< internal timestep for particle advection
    LOGICAL ::  dt_3d_reached_l                                  !< internal timestep for particle advection
    LOGICAL ::  force_print_header = .FALSE.                     !< namelist parameter
    LOGICAL ::  les_mw = .FALSE.                                 !< use Moeng-Wyngaard turbulence closure for LES mode
    LOGICAL ::  linear_eqnOfState = .FALSE.                      !< namelist parmaeter for linear equation of state in ocean
    REAL(wp) :: rho_ref = 1000.0_wp                              !< reference density for linear eos
    REAL(wp) :: alpha_const = 2.0E-4                             !< fixed alpha_T value
    REAL(wp) :: beta_const = 8.0E-4                              !< fixed beta_S value
    REAL(wp) :: pt_ref = 15.0_wp                                 !< potential temperature reference falue
    REAL(wp) :: sa_ref = 35.0_wp                                 !< salinity reerence value for fixed linear density equation
    LOGICAL ::  idealized_diurnal = .FALSE.                      !< flag for diurnal cycle
    LOGICAL ::  run_control_header = .FALSE.                     !< onetime output of RUN_CONTROL header
    LOGICAL ::  scalar_rayleigh_damping = .TRUE.                 !< namelist parameter
    LOGICAL ::  stokes_force = .FALSE.                           !< switch for use of Stokes forces
    LOGICAL ::  stop_dt = .FALSE.                                !< internal switch to stop the time stepping
    LOGICAL ::  synchronous_exchange = .FALSE.                   !< namelist parameter
    LOGICAL ::  syn_turb_gen = .FALSE.                           !< flag for synthetic turbulence generator module
    LOGICAL ::  terminate_run = .FALSE.                          !< terminate run (cpu-time limit, restarts)?
    LOGICAL ::  transpose_compute_overlap = .FALSE.              !< namelist parameter
    LOGICAL ::  use_initial_profile_as_reference = .FALSE.       !< use of initial profiles as reference state?
    LOGICAL ::  use_prescribed_profile_data = .FALSE.            !< use of prescribed wind profiles?
                                                                 !< (namelist parameters u_profile, v_profile)
    LOGICAL ::  wall_adjustment = .TRUE.                         !< namelist parameter
    LOGICAL ::  write_binary = .FALSE.                           !< ENVPAR namelist parameter to steer restart I/O (ENVPAR is created by palmrun)

    REAL(wp) :: ideal_solar_heatflux = 0.0_wp                    !< maximum daytime solar heatflux (same sign convention as top_heatflux)
    REAL(wp) :: ideal_solar_division = 0.67_wp                   !< value for breakdown of double exponential
    REAL(wp) :: ideal_solar_efolding1 = 1.0_wp/1.0_wp            !< efolding depth for IR in solar (m^-1)
    REAL(wp) :: ideal_solar_efolding2 = 1.0_wp/17.0_wp           !< efolding depth for blue in solar (m^-1)
    REAL(wp) ::  averaging_interval = 0.0_wp                   !< namelist parameter
    REAL(wp) ::  averaging_interval_pr = 9999999.9_wp          !< namelist parameter
    REAL(wp) ::  bc_pt_t_val                                   !< vertical gradient of pt near domain top
    REAL(wp) ::  cfl_factor = -1.0_wp                          !< namelist parameter
    REAL(wp) ::  disturbance_amplitude = 0.25_wp               !< namelist parameter
    REAL(wp) ::  disturbance_energy_limit = 0.01_wp            !< namelist parameter
    REAL(wp) ::  disturbance_level_b = -9999999.9_wp           !< namelist parameter
    REAL(wp) ::  disturbance_level_t = -9999999.9_wp           !< namelist parameter
    REAL(wp) ::  dp_level_b = 0.0_wp                           !< namelist parameter
    REAL(wp) ::  dt = -1.0_wp                                  !< namelist parameter
    REAL(wp) ::  dt_averaging_input = 0.0_wp                   !< namelist parameter
    REAL(wp) ::  dt_averaging_input_pr = 9999999.9_wp          !< namelist parameter
    REAL(wp) ::  dt_data_output = 3600.0_wp                 !< namelist parameter
    REAL(wp) ::  dt_data_output_av = 3600.0_wp              !< namelist parameter
    REAL(wp) ::  dt_disturb = 20.0_wp                     !< namelist parameter
    REAL(wp) ::  dt_dopr = 3600.0_wp                        !< namelist parameter
    REAL(wp) ::  dt_dopr_listing = 9999999.9_wp                !< namelist parameter
    REAL(wp) ::  dt_do3d = 9999999.9_wp                        !< namelist parameter
    REAL(wp) ::  dt_max = 20.0_wp                              !< namelist parameter
    REAL(wp) ::  dt_restart = 9999999.9_wp                     !< namelist parameter
    REAL(wp) ::  dt_run_control = 0.0_wp                      !< namelist parameter
    REAL(wp) ::  dt_3d = 0.01_wp                               !< time step
#ifdef __SP
    REAL(wp) ::  dt_LS = 60.0_wp                               !< large scale time step
    REAL(wp) ::  factor_LS_tracer = 0.0_wp                     !< large scale frocing factor for tracers
    REAL(wp) ::  factor_LS_vel = 0.0_wp                        !< large scale frocing factor for velocity
#endif

    REAL(wp) ::  dz_max = 1000.0_wp                            !< namelist parameter
    REAL(wp) ::  dz_stretch_factor = 1.08_wp                   !< namelist parameter
    REAL(wp) ::  dz_stretch_level = -9999999.9_wp              !< namelist parameter
    REAL(wp) ::  e_init = 0.0_wp                               !< namelist parameter
    REAL(wp) ::  e_min = 0.0_wp                                !< namelist parameter
    REAL(wp) ::  end_time = 10800.0_wp                             !< namelist parameter
    REAL(wp) ::  f = 0.0_wp                                    !< Coriolis parameter
    REAL(wp) ::  fs = 0.0_wp                                   !< Coriolis parameter
    REAL(wp) ::  g = 9.81_wp                                   !< gravitational acceleration
    REAL(wp) ::  kappa = 0.4_wp                                !< von Karman constant
    REAL(wp) ::  latitude = 0.0_wp                             !< namelist parameter
    REAL(wp) ::  longitude = 0.0_wp                            !< namelist parameter
    REAL(wp) ::  maximum_cpu_time_allowed = 0.0_wp             !< given wall time for run
    REAL(wp) ::  old_dt = 1.0E-10_wp                           !< length of previous timestep
    REAL(wp) ::  omega = 7.29212E-5_wp                         !< namelist parameter
    REAL(wp) ::  particle_maximum_age = 9999999.9_wp           !< namelist parameter
    REAL(wp) ::  prandtl_number = 1.0_wp                       !< namelist parameter
    REAL(wp) ::  prho_reference                                !< reference state of potential density
    REAL(wp) ::  pt_surface = 293.0_wp                         !< namelist parameter
    REAL(wp) ::  pt_surface_initial_change = 0.0_wp            !< namelist parameter
    REAL(wp) ::  rayleigh_damping_factor = -1.0_wp             !< namelist parameter
    REAL(wp) ::  rayleigh_damping_height = -1.0_wp             !< namelist parameter
    REAL(wp) ::  restart_time = 9999999.9_wp                   !< namelist parameter
    REAL(wp) ::  rho_surface                                   !< surface value of density
    REAL(wp) ::  sa_surface = 35.0_wp                          !< namelist parameter
    REAL(wp) ::  simulated_time = 0.0_wp                       !< elapsed simulated time
    REAL(wp) ::  simulated_time_at_begin                       !< elapsed simulated time of previous run (job chain)
    REAL(wp) ::  skip_time_data_output = 0.0_wp                !< namelist parameter
    REAL(wp) ::  skip_time_data_output_av = 9999999.9_wp       !< namelist parameter
    REAL(wp) ::  skip_time_dopr = 9999999.9_wp                 !< namelist parameter
    REAL(wp) ::  skip_time_do3d = 9999999.9_wp                 !< namelist parameter
    REAL(wp) ::  surface_pressure = 1013.25_wp                 !< namelist parameter
    REAL(wp) ::  termination_time_needed = 35.0_wp             !< namelist parameter
    REAL(wp) ::  time_disturb = 0.0_wp                         !< time since last flow disturbance
    REAL(wp) ::  time_dopr = 0.0_wp                            !< time since last profile output
    REAL(wp) ::  time_dopr_av = 0.0_wp                         !< time since last averaged profile output
    REAL(wp) ::  time_dopr_listing = 0.0_wp                    !< time since last profile output (ASCII) on file
    REAL(wp) ::  time_dosp = 0.0_wp                            !< time since last spectra output
    REAL(wp) ::  time_dosp_av = 0.0_wp                         !< time since last averaged spectra output
    REAL(wp) ::  time_do3d = 0.0_wp                            !< time since last 3d output
    REAL(wp) ::  time_do_av = 0.0_wp                           !< time since last averaged-data output
    REAL(wp) ::  time_do_sla = 0.0_wp                          !< time since last
    REAL(wp) ::  time_restart = 9999999.9_wp                   !< time at which run shall be terminated and restarted
    REAL(wp) ::  time_run_control = 0.0_wp                     !< time since last RUN_CONTROL output
    REAL(wp) ::  time_since_reference_point = 0.0_wp           !< time after atmosphere-ocean coupling has been activated, or time after spinup phase of LSM has been finished
    REAL(wp) ::  top_heatflux = 5e-5_wp                   !< namelist parameter
    REAL(wp) ::  top_momentumflux_u = 0.0_wp             !< namelist parameter
    REAL(wp) ::  top_momentumflux_v = 0.0_wp             !< namelist parameter
    REAL(wp) ::  top_salinityflux = 0.0_wp               !< namelist parameter
    REAL(wp) ::  wall_adjustment_factor = 1.8_wp               !< adjustment factor for mixing length l

    REAL(wp) ::  dpdxy(1:2) = 0.0_wp                               !< namelist parameter
    REAL(wp) ::  dz(10) = -1.0_wp                                  !< namelist parameter
    REAL(wp) ::  dzconst = 2.5_wp
    REAL(wp) ::  dz_stretch_level_start(9) = -9999999.9_wp         !< namelist parameter
    REAL(wp) ::  dz_stretch_level_end(9) = 9999999.9_wp            !< namelist parameter
    REAL(wp) ::  dz_stretch_factor_array(9) = 1.08_wp              !< namelist parameter
    REAL(wp) ::  pt_vertical_gradient(10) = 0.0_wp                 !< namelist parameter
    REAL(wp) ::  pt_vertical_gradient_level(10) = -999999.9_wp     !< namelist parameter
    REAL(wp) ::  sa_vertical_gradient(10) = 0.0_wp                 !< namelist parameter
    REAL(wp) ::  sa_vertical_gradient_level(10) = -999999.9_wp     !< namelist parameter
    REAL(wp) ::  tsc(10) = (/ 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, &    !< array used for controlling time-integration at different substeps
                 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp /)
    REAL(wp) ::  u_profile(100) = 9999999.9_wp                     !< namelist parameter
    REAL(wp) ::  uv_heights(100) = 9999999.9_wp                    !< namelist parameter
    REAL(wp) ::  v_profile(100) = 9999999.9_wp                     !< namelist parameter
    REAL(wp) ::  u0_stk = -9999999.9_wp                            !< surface Stokes drift in m/s, x-component
    REAL(wp) ::  v0_stk = -9999999.9_wp                            !< surface Stokes drift in m/s, y-component
    REAL(wp) ::  d_stk = -9999999.9_wp                             !< Stokes drift exponential decay depth scale in m
    REAL(wp) ::  wind_speed = -9999999.9_wp                        !< 10-meter wind speed used to compute empirical wave spectra (m/s)
    REAL(wp) ::  wind_dir = -9999999.9_wp                          !< wind direction used to compute empirical wave spectra, degree counter-clockwise from x-direction
    REAL(wp) ::  wave_age = -9999999.9_wp                          !< wave age c_p/(U_{10}\cos\theta) used to compute empirical wave spectra

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  dp_smooth_factor  !< smoothing factor for external pressure gradient forcing

!    SAVE

 END MODULE control_parameters


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of grid spacings.
!------------------------------------------------------------------------------!
 MODULE grid_variables

    USE kinds

    REAL(wp) ::  ddx          !< 1/dx
    REAL(wp) ::  ddx2         !< 1/dx2
    REAL(wp) ::  dx = 2.5_wp  !< horizontal grid size (along x-direction)
    REAL(wp) ::  dx2          !< dx*dx
    REAL(wp) ::  ddy          !< 1/dy
    REAL(wp) ::  ddy2         !< 1/dy2
    REAL(wp) ::  dy = 2.5_wp  !< horizontal grid size (along y-direction)
    REAL(wp) ::  dy2          !< dy*dy

!    SAVE

 END MODULE grid_variables


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of array bounds, number of gridpoints, and wall flag arrays.
!------------------------------------------------------------------------------!
 MODULE indices

    USE kinds

    INTEGER(iwp) ::  nbgp = 3       !< number of boundary ghost points
    INTEGER(iwp) ::  nnx            !< number of subdomain grid points in x-direction
    INTEGER(iwp) ::  nx = 31         !< nx+1 = total number of grid points in x-direction
    INTEGER(iwp) ::  nxl            !< left-most grid index of subdomain (excluding ghost points)
    INTEGER(iwp) ::  nxlg           !< left-most grid index of subdomain (including ghost points)
    INTEGER(iwp) ::  nxlu           !< =nxl+1 (at left domain boundary with inflow from left), else =nxl (used for u-velocity component)
    INTEGER(iwp) ::  nxr            !< right-most grid index of subdomain (excluding ghost points)
    INTEGER(iwp) ::  nxrg           !< right-most grid index of subdomain (including ghost points)
    INTEGER(iwp) ::  nx_on_file     !< nx of previous run in job chain
    INTEGER(iwp) ::  nny            !< number of subdomain grid points in y-direction
    INTEGER(iwp) ::  ny = 31         !< ny+1 = total number of grid points in y-direction
    INTEGER(iwp) ::  nyn            !< north-most grid index of subdomain (excluding ghost points)
    INTEGER(iwp) ::  nyng           !< north-most grid index of subdomain (including ghost points)
    INTEGER(iwp) ::  nys            !< south-most grid index of subdomain (excluding ghost points)
    INTEGER(iwp) ::  nysg           !< south-most grid index of subdomain (including ghost points)
    INTEGER(iwp) ::  nysv           !< =nys+1 (at south domain boundary with inflow from south), else =nys (used for v-velocity component)
    INTEGER(iwp) ::  ny_on_file     !< ny of previous run in job chain
    INTEGER(iwp) ::  nnz            !< number of subdomain grid points in z-direction
    INTEGER(iwp) ::  nz = 256         !< total number of grid points in z-direction
    INTEGER(iwp) ::  nzb            !< bottom grid index of computational domain
    INTEGER(iwp) ::  nzt            !< nzt+1 = top grid index of computational domain

    INTEGER(idp) ::  ngp_3d        !< number of grid points of the total domain

    INTEGER(iwp) ::  ngp_2dh  !< number of grid points of a horizontal cross section through the total domain


    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  advc_flags_1            !< flags used to degrade order of advection scheme
    INTEGER(iwp), DIMENSION(:,:,:), ALLOCATABLE ::  advc_flags_2            !< flags used to degrade order of advection scheme

!    SAVE

 END MODULE indices


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Interfaces for special subroutines which use optional parameters.
!------------------------------------------------------------------------------!
 MODULE interfaces

    INTERFACE

!------------------------------------------------------------------------------!
! Description:
! ------------
!> @todo Missing subroutine description.
!------------------------------------------------------------------------------!
       SUBROUTINE global_min_max ( i1, i2, j1, j2, k1, k2, array, mode, offset, &
                                   result, result_ijk, result1, result1_ijk )

          USE kinds

          CHARACTER (LEN=*), INTENT(IN) ::  mode                      !< mode of global min/max function: can be 'min', 'max', 'minmax', 'abs', or 'absoff'
          INTEGER(iwp), INTENT(IN)      ::  i1                        !< internal index of min/max function
          INTEGER(iwp), INTENT(IN)      ::  i2                        !< internal index of min/max function
          INTEGER(iwp), INTENT(IN)      ::  j1                        !< internal index of min/max function
          INTEGER(iwp), INTENT(IN)      ::  j2                        !< internal index of min/max function
          INTEGER(iwp), INTENT(IN)      ::  k1                        !< internal index of min/max function
          INTEGER(iwp), INTENT(IN)      ::  k2                        !< internal index of min/max function
          INTEGER(iwp)                  ::  result_ijk(3)             !< grid index result of min/max function
          INTEGER(iwp), OPTIONAL        ::  result1_ijk(3)            !< optional grid index result of min/max function
          REAL(wp)                      ::  offset                    !< min/max function calculates absolute value with respect to an offset
          REAL(wp)                      ::  result                    !< result of min/max function
          REAL(wp), OPTIONAL            ::  result1                   !< optional result of min/max function
          REAL(wp), INTENT(IN)          ::  array(i1:i2,j1:j2,k1:k2)  !< input array of min/max function

       END SUBROUTINE global_min_max

    END INTERFACE

!    SAVE

 END MODULE interfaces


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of variables which define processor topology and the exchange of
!> ghost point layers. This module must be placed in all routines containing
!> MPI-calls.
!------------------------------------------------------------------------------!
 MODULE pegrid

    USE kinds

#if defined( __parallel )
#if defined( __mpifh )
    INCLUDE "mpif.h"
#else
    USE MPI
#endif
#endif
    CHARACTER(LEN=2) ::  send_receive = 'al'     !<
    CHARACTER(LEN=7) ::  myid_char = ''          !< character string containing processor id number

    INTEGER(iwp) ::  comm1dx                     !< communicator for domain decomposition along x
    INTEGER(iwp) ::  comm1dy                     !< communicator for domain decomposition along y
    INTEGER(iwp) ::  comm2d                      !< standard 2d (xy) communicator used in PALM for the process group the PE belongs to
    INTEGER(iwp) ::  comm_inter                  !< intercommunicator that connects atmosphere/ocean process groups
    INTEGER(iwp) ::  comm_palm                   !< internal communicator used during the MPI setup at the beginning of a run
    INTEGER(iwp) ::  ierr                        !< standard error parameter in MPI calls
    INTEGER(iwp) ::  myid = 0                    !< id number of processor element
    INTEGER(iwp) ::  myidx = 0                   !< id number of processor elements with same position along x-direction
    INTEGER(iwp) ::  myidy = 0                   !< id number of processor elements with same position along y-direction
    INTEGER(iwp) ::  ndim = 2                    !< dimension of the virtual PE grid
    INTEGER(iwp) ::  ngp_y                       !< number of subdomain grid points along y including ghost points
    INTEGER(iwp) ::  npex = -1                   !< number of processor elements in x-direction
    INTEGER(iwp) ::  npey = -1                   !< number of processor elements in y-direction
    INTEGER(iwp) ::  numprocs = 1                !< total number of appointed processor elements
    INTEGER(iwp) ::  numprocs_previous_run = -1  !< total number of appointed processor elements in previous run (job chain)
    INTEGER(iwp) ::  pleft                       !< MPI-address of the processor left of the current one
    INTEGER(iwp) ::  pnorth                      !< MPI-address of the processor north of the current one
    INTEGER(iwp) ::  pright                      !< MPI-address of the processor right of the current one
    INTEGER(iwp) ::  psouth                      !< MPI-address of the processor south of the current one
    INTEGER(iwp) ::  req_count = 0               !< MPI return variable - checks if Send-Receive operation is already finished
    INTEGER(iwp) ::  run_id = -9999              !< id number of run
    INTEGER(iwp) ::  sendrecvcount_xy            !< number of subdomain gridpoints to be exchanged in direct transpositions (y --> x, or x --> y) or second (2d) transposition x --> y
    INTEGER(iwp) ::  sendrecvcount_yz            !< number of subdomain gridpoints to be exchanged in third (2d) transposition y --> z
    INTEGER(iwp) ::  sendrecvcount_zx            !< number of subdomain gridpoints to be exchanged in first (2d) transposition z --> x
    INTEGER(iwp) ::  sendrecvcount_zyd           !< number of subdomain gridpoints to be exchanged in direct transpositions z --> y (used for calculating spectra)
    INTEGER(iwp) ::  tasks_per_node = -9999      !< MPI tasks per compute node
    INTEGER(iwp) ::  threads_per_task = 1        !< number of OPENMP threads per MPI task
    INTEGER(iwp) ::  type_x                      !< derived MPI datatype for 2-D ghost-point exchange - north / south
    INTEGER(iwp) ::  type_xy                     !< derived MPI datatype for 2-D ghost-point exchange - north / south
    INTEGER(iwp) ::  type_y                      !< derived MPI datatype for 2-D exchange in atmosphere-ocean coupler

    INTEGER(iwp) ::  pdims(2) = 1  !< number of processors along x-y dimension
    INTEGER(iwp) ::  req(100)      !< MPI return variable indicating if send-receive operation is finished

    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  hor_index_bounds               !< horizontal index bounds
    INTEGER(iwp), DIMENSION(:,:), ALLOCATABLE ::  hor_index_bounds_previous_run  !< horizontal index bounds of previous run

    LOGICAL ::  collective_wait = .FALSE.          !< switch to set an explicit MPI barrier in front of all collective MPI calls

#if defined( __parallel )
    INTEGER(iwp) ::  ibuf(12)                 !< internal buffer for calculating MPI settings
    INTEGER(iwp) ::  pcoord(2)                !< PE coordinates along x and y
    INTEGER(iwp) ::  status(MPI_STATUS_SIZE)  !< MPI status variable used in various MPI calls

    INTEGER(iwp), DIMENSION(MPI_STATUS_SIZE,100) ::  wait_stat  !< MPI status variable used in various MPI calls

    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ngp_xz      !< number of ghost points in xz-plane on different multigrid level
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ngp_xz_int  !< number of ghost points in xz-plane on different multigrid level
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ngp_yz      !< number of ghost points in yz-plane on different multigrid level
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  ngp_yz_int  !< number of ghost points in yz-plane on different multigrid level
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  type_x_int  !< derived MPI datatype for 2-D integer ghost-point exchange - north / south
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  type_xz     !< derived MPI datatype for 3-D integer ghost-point exchange - north / south
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  type_xz_int !< derived MPI datatype for 3-D integer ghost-point exchange - north / south
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  type_y_int  !< derived MPI datatype for 2-D integer ghost-point exchange - left / right
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  type_yz     !< derived MPI datatype for 3-D integer ghost-point exchange - left / right
    INTEGER(iwp), DIMENSION(:), ALLOCATABLE ::  type_yz_int !< derived MPI datatype for 3-D integer ghost-point exchange - left / right

    LOGICAL ::  reorder = .TRUE.           !< switch to allow MPI the reorder of ranking (e.g. row-major or column-major)

    LOGICAL, DIMENSION(2) ::  cyclic = (/ .TRUE. , .TRUE. /)  !< boundary conditions of the virtual PE grid
    LOGICAL, DIMENSION(2) ::  remain_dims                     !< internal array used to determine sub-topologies for transpositions
#endif

!    SAVE

 END MODULE pegrid


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of variables which control PROFIL-output.
!------------------------------------------------------------------------------!
 MODULE profil_parameter

    USE kinds

    INTEGER(iwp), PARAMETER ::  crmax = 50  !< maximum number of coordinate systems for profile output

    CHARACTER (LEN=100), DIMENSION(crmax) ::  cross_profiles = &  !< quantities to be plotted into one coordinate system, respectively
                           (/ ' u v                           ', &
                              ' pt                            ', &
                              ' w"pt" w*pt*                   ', &
                              ' w"u" w*u* w"v" w*v*           ', &
                              ' km kh                         ', &
                            ( '                               ', i9 = 1, 45 ) /)

    INTEGER(iwp) ::  profile_columns = 2  !< number of coordinate systems on a profile plot per column
    INTEGER(iwp) ::  profile_rows = 3     !< number of coordinate systems on a profile plot per row

    INTEGER(iwp) ::  dopr_index(100) = 0                !< index number of respective profile quantity
    INTEGER(iwp) ::  dopr_initial_index(100) = 0        !< index number of initial profiles to be output

!    SAVE


 END MODULE profil_parameter

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of statistical quantities, e.g. global sums.
!------------------------------------------------------------------------------!
 MODULE statistics

    USE kinds

    INTEGER(iwp) ::  pr_palm = 41          !< maximum number of output profiles

    INTEGER(iwp) ::  u_max_ijk(3) = -1  !< index values (i,j,k) of location where u_max occurs
    INTEGER(iwp) ::  v_max_ijk(3) = -1  !< index values (i,j,k) of location where v_max occurs
    INTEGER(iwp) ::  w_max_ijk(3) = -1  !< index values (i,j,k) of location where w_max occurs

    LOGICAL ::  flow_statistics_called = .FALSE.  !< flag that tells other routines if flow statistics was executed
                                                  !< (after each timestep)

    REAL(wp) ::  u_max = 0.0_wp  !< maximum of absolute u-veloctiy in entire domain
    REAL(wp) ::  v_max = 0.0_wp  !< maximum of absolute v-veloctiy in entire domain
    REAL(wp) ::  w_max = 0.0_wp  !< maximum of absolute w-veloctiy in entire domain

    REAL(wp), DIMENSION(:), ALLOCATABLE ::  weight_substep             !< weighting factor for substeps in timestepping
    REAL(wp), DIMENSION(:), ALLOCATABLE ::  weight_pres                !< substep weighting factor for pressure solver

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sums             !< global sum array for the various output quantities
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sums_wsus_ws_l   !< subdomain sum of vertical momentum flux w'u' (5th-order advection scheme only)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sums_wsvs_ws_l   !< subdomain sum of vertical momentum flux w'v' (5th-order advection scheme only)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sums_vsus_ws_l   !< subdomain sum of horizontal momentum flux v'u' (5th-order advection scheme only)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sums_usvs_ws_l   !< subdomain sum of horizontal momentum flux u'v' (5th-order advection scheme only)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sums_us2_ws_l    !< subdomain sum of horizontal momentum flux u'u' (5th-order advection scheme only)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sums_vs2_ws_l    !< subdomain sum of horizontal momentum flux v'v' (5th-order advection scheme only)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sums_ws2_ws_l    !< subdomain sum of vertical momentum flux w'w' (5th-order advection scheme only)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sums_wspts_ws_l  !< subdomain sum of vertical sensible heat flux w'pt' (5th-order advection scheme only)
    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  sums_wssas_ws_l  !< subdomain sum of vertical salinity flux w'sa' (5th-order advection scheme only)

    REAL(wp), DIMENSION(:,:), ALLOCATABLE ::  hom_sum             !< sum array for horizontal mean
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  sums_l              !< subdomain sum (_l) gathered for various quantities

    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  hom  !< horizontal mean of various quantities (profiles/timeseries)

!    SAVE

 END MODULE statistics


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Definition of indices for transposed arrays.
!------------------------------------------------------------------------------!
 MODULE transpose_indices

    USE kinds

    INTEGER(iwp) ::  nxl_y   !< internal index bound for transpositions
    INTEGER(iwp) ::  nxl_yd  !< internal index bound for transpositions
    INTEGER(iwp) ::  nxl_z   !< internal index bound for transpositions
    INTEGER(iwp) ::  nxr_y   !< internal index bound for transpositions
    INTEGER(iwp) ::  nxr_yd  !< internal index bound for transpositions
    INTEGER(iwp) ::  nxr_z   !< internal index bound for transpositions
    INTEGER(iwp) ::  nyn_x   !< internal index bound for transpositions
    INTEGER(iwp) ::  nyn_z   !< internal index bound for transpositions
    INTEGER(iwp) ::  nys_x   !< internal index bound for transpositions
    INTEGER(iwp) ::  nys_z   !< internal index bound for transpositions
    INTEGER(iwp) ::  nzb_x   !< internal index bound for transpositions
    INTEGER(iwp) ::  nzb_y   !< internal index bound for transpositions
    INTEGER(iwp) ::  nzb_yd  !< internal index bound for transpositions
    INTEGER(iwp) ::  nzt_x   !< internal index bound for transpositions
    INTEGER(iwp) ::  nzt_y   !< internal index bound for transpositions
    INTEGER(iwp) ::  nzt_yd  !< internal index bound for transpositions

!    SAVE

 END MODULE transpose_indices
