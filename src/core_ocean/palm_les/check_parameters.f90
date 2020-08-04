!> @file check_parameters.f90
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
! -----------------
!
! 2018-11-01 cbegeman
! Add checks for bubble initial conditions
!
! 2018-10-31 cbegeman
! Add checks for profile output
!
! 2018-10-25 cbegeman
! Add checks for dirichlet bottom boundary conditions for salinity
!
! Former revisions:
! -----------------
! $Id: check_parameters.f90 2520 2017-10-05 13:50:26Z gronemeier &
! Add inital profile output for e (TG)
!
! 3065 2018-06-12 07:03:02Z Giersch
! dz was replaced by dz(1), error message revised
!
! 3049 2018-05-29 13:52:36Z Giersch
! add variable description
!
! 3046 2018-05-29 08:02:15Z Giersch
! Error messages revised
!
! 3045 2018-05-28 07:55:41Z Giersch
! Error messages revised
!
! 3035 2018-05-24 09:35:20Z schwenkel
! Add option to initialize warm air bubble close to surface
!
! 3034 2018-05-24 08:41:20Z raasch
! bugfix: check that initializing_actions has been set
!
! 2980 2018-04-17 15:19:27Z suehring
! Further improvement for spinup checks.
!
! 2974 2018-04-16 12:59:52Z gronemeier
! Bugfix: check if dt_data_output_av is zero in case of parallel NetCDF output
!
! 2970 2018-04-13 15:09:23Z suehring
! Bugfix in old large-scale forcing mode
!
! 2964 2018-04-12 16:04:03Z Giersch
! Calculation of fixed number of output time levels for parallel netcdf output
! has been revised (based on calculations in netcdf_interface_mod)
!
! 2938 2018-03-27 15:52:42Z suehring
! - Revise start and end indices for imposing random disturbances in case of
!   nesting or large-scale forcing.
! - Remove check for inifor initialization in case of nested runs.
! - Adapt call to check for synthetic turbulence geneartor settings.
!
! 2936 2018-03-27 14:49:27Z suehring
! Check spinup in case of nested runs, spinup time and timestep must be
! identical to assure synchronuous simulations.
!
! 2932 2018-03-26 09:39:22Z maronga
! renamed particles_par to particle_parameters
!
! 2918 2018-03-21 15:52:14Z gronemeier
! Add check for 1D model
!
! 2883 2018-03-14 08:29:10Z Giersch
! dt_dopr_listing is not set to the default value zero anymore
!
! 2851 2018-03-05 14:39:31Z maronga
! Bugfix: calculation of output time levels in case of restart runs (parallel
! NetCDF)
!
! 2836 2018-02-26 13:40:05Z Giersch
! dt_dopr_listing is set to the default value zero
!
! 2817 2018-02-19 16:32:21Z knoop
! Preliminary gust module interface implemented
!
! 2798 2018-02-09 17:16:39Z suehring
! Consider also default-type surfaces for surface temperature output.
!
! 2797 2018-02-08 13:24:35Z suehring
! Enable output of ground-heat flux also at urban surfaces.
!
! 2776 2018-01-31 10:44:42Z Giersch
! Variable synthetic_turbulence_generator has been abbreviated
!
! 2773 2018-01-30 14:12:54Z suehring
! Check for consistent initialization in nesting mode added.
!
! 2766 2018-01-22 17:17:47Z kanani
! Removed preprocessor directive __chem
!
! 2765 2018-01-22 11:34:58Z maronga
! Renamed simulation_time_since_reference to
! time_to_be_simulated_from_reference_point
!
! 2746 2018-01-15 12:06:04Z suehring
! Move flag plant canopy to modules
!
! 2743 2018-01-12 16:03:39Z suehring
! In case of natural- and urban-type surfaces output surfaces fluxes in W/m2.
!
! 2742 2018-01-12 14:59:47Z suehring
! Enable output of surface temperature
!
! 2735 2018-01-11 12:01:27Z suehring
! output of r_a moved from land-surface to consider also urban-type surfaces
!
! 2718 2018-01-02 08:49:38Z maronga
! Corrected "Former revisions" section
!
! 2696 2017-12-14 17:12:51Z kanani
! Change in file header (GPL part)
! Implementation of uv exposure model (FK)
! + new possible value for dissipation_1d
! Added checks for turbulence_closure_mod (TG)
! Implementation of chemistry module (FK)
!
! 2689 2017-12-12 17:46:55Z Giersch
! Bugfix in if query
!
! 2688 2017-12-12 17:27:04Z Giersch
! Check if humidity is set to TRUE in the _p3d file for coupled runs
!
! 2669 2017-12-06 16:03:27Z raasch
! mrun-string replaced by palmrun
!
! 2628 2017-11-20 12:40:38Z schwenkel
! Enabled particle advection with grid stretching -> Removed parameter check
!
! 2575 2017-10-24 09:57:58Z maronga
! Renamed phi --> latitude
!
! 2564 2017-10-19 15:56:56Z Giersch
! Variable wind_turbine was added to control_parameters.
!
! 2550 2017-10-16 17:12:01Z boeske
! Added checks for complex terrain simulations
!
! 2513 2017-10-04 09:24:39Z kanani
! Bugfix for some dopr(_initial)_index values and units connected to
! passive-scalar output
!
! 2508 2017-10-02 08:57:09Z suehring
! Bugfix, change default value of vertical_gradient level in order to consider
! also ocean runs
!
! 2422 2017-09-08 08:25:41Z raasch
! error message in case of missing "restart" file activation string
!
! 2375 2017-08-29 14:10:28Z schwenkel
! Added aerosol for bulk microphysics
!
! 2365 2017-08-21 14:59:59Z kanani
! Vertical grid nesting implemented: Check coupling mode. Generate file header
! (SadiqHuq)
!
! 2354 2017-08-17 10:49:36Z schwenkel
! Bugfix correlated to lsm_check_data_output_pr.
! If-statement for following checks is essential, otherwise units for lsm output
! are set to 'illegal' and palm will be aborted.
!
! 2348 2017-08-10 10:40:10Z kanani
! New: Check for simultaneous use of geostrophic wind and u_profile/v_profile
!
! 2345 2017-08-09 11:50:30Z Giersch
! Remove error message PA0156 and the conserve_volume_flow_mode option
! inflow_profile
!
! 2339 2017-08-07 13:55:26Z raasch
! corrected timestamp in header
!
! 2338 2017-08-07 12:15:38Z gronemeier
! Modularize 1D model
!
! 2329 2017-08-03 14:24:56Z knoop
! Bugfix: index corrected for rho_air and rho_air_zw output
!
! 2320 2017-07-21 12:47:43Z suehring
! Modularize large-scale forcing and nudging
!
! 2312 2017-07-14 20:26:51Z hoffmann
! PA0349 and PA0420 removed.
!
! 2300 2017-06-29 13:31:14Z raasch
! host-specific settings and checks removed
!
! 2292 2017-06-20 09:51:42Z schwenkel
! Implementation of new microphysic scheme: cloud_scheme = 'morrison'
! includes two more prognostic equations for cloud drop concentration (nc)
! and cloud water content (qc).
!
! 2274 2017-06-09 13:27:48Z Giersch
! Changed error messages
!
! 2271 2017-06-09 12:34:55Z sward
! roughness-length check altered
! Error messages fixed
!
! 2270 2017-06-09 12:18:47Z maronga
! Revised numbering (removed 2 timeseries)
!
! 2259 2017-06-08 09:09:11Z gronemeier
! Implemented synthetic turbulence generator
!
! 2251 2017-06-06 15:10:46Z Giersch
!
! 2250 2017-06-06 15:07:12Z Giersch
! Doxygen comment added
!
! 2248 2017-06-06 13:52:54Z sward
! Error message fixed
!
! 2209 2017-04-19 09:34:46Z kanani
! Check for plant canopy model output
!
! 2200 2017-04-11 11:37:51Z suehring
! monotonic_adjustment removed
!
! 2178 2017-03-17 11:07:39Z hellstea
! Index limits for perturbations are now set also in case of nested
! boundary conditions
!
! 2169 2017-03-06 18:16:35Z suehring
! Bugfix, move setting for topography grid convention to init_grid
!
! 2118 2017-01-17 16:38:49Z raasch
! OpenACC related parts of code removed
!
! 2088 2016-12-19 16:30:25Z suehring
! Bugfix in initial salinity profile
!
! 2084 2016-12-09 15:59:42Z knoop
! Removed anelastic + multigrid error message
!
! 2052 2016-11-08 15:14:59Z gronemeier
! Bugfix: remove setting of default value for recycling_width
!
! 2050 2016-11-08 15:00:55Z gronemeier
! Implement turbulent outflow condition
!
! 2044 2016-11-02 16:44:25Z knoop
! Added error code for anelastic approximation
!
! 2042 2016-11-02 13:47:31Z suehring
! Additional checks for wall_heatflux, wall_humidityflux and wall_scalarflux.
! Bugfix, check for constant_scalarflux.
!
! 2040 2016-10-26 16:58:09Z gronemeier
! Removed check for statistic_regions > 9.
!
! 2037 2016-10-26 11:15:40Z knoop
! Anelastic approximation implemented
!
! 2031 2016-10-21 15:11:58Z knoop
! renamed variable rho to rho_ocean
!
! 2026 2016-10-18 10:27:02Z suehring
! Bugfix, enable output of s*2
!
! 2024 2016-10-12 16:42:37Z kanani
! Added missing CASE, error message and unit for ssws*,
! increased number of possible output quantities to 500.
!
! 2011 2016-09-19 17:29:57Z kanani
! Flag urban_surface is now defined in module control_parameters,
! changed prefix for urban surface model output to "usm_",
! added flag lsf_exception (inipar-Namelist parameter) to allow explicit
! enabling of large scale forcing together with buildings on flat terrain,
! introduced control parameter varnamelength for LEN of var.
!
! 2007 2016-08-24 15:47:17Z kanani
! Added checks for the urban surface model,
! increased counter in DO WHILE loop over data_output (for urban surface output)
!
! 2000 2016-08-20 18:09:15Z knoop
! Forced header and separation lines into 80 columns
!
! 1994 2016-08-15 09:52:21Z suehring
! Add missing check for cloud_physics and cloud_droplets
!
! 1992 2016-08-12 15:14:59Z suehring
! New checks for top_scalarflux
! Bugfixes concerning data output of profiles in case of passive_scalar
!
! 1984 2016-08-01 14:48:05Z suehring
! Bugfix: checking for bottom and top boundary condition for humidity and scalars
!
! 1972 2016-07-26 07:52:02Z maronga
! Removed check of lai* (done in land surface module)
!
! 1970 2016-07-18 14:27:12Z suehring
! Bugfix, check for ibc_q_b and constant_waterflux in case of land-surface scheme
!
! 1962 2016-07-13 09:16:44Z suehring
! typo removed
!
! 1960 2016-07-12 16:34:24Z suehring
! Separate checks and calculations for humidity and passive scalar. Therefore,
! a subroutine to calculate vertical gradients is introduced, in order  to reduce
! redundance.
! Additional check - large-scale forcing is not implemented for passive scalars
! so far.
!
! 1931 2016-06-10 12:06:59Z suehring
! Rename multigrid into multigrid_noopt and multigrid_fast into multigrid
!
! 1929 2016-06-09 16:25:25Z suehring
! Bugfix in check for use_upstream_for_tke
!
! 1918 2016-05-27 14:35:57Z raasch
! setting of a fixed reference state ('single_value') for ocean runs removed
!
! 1914 2016-05-26 14:44:07Z witha
! Added checks for the wind turbine model
!
! 1841 2016-04-07 19:14:06Z raasch
! redundant location message removed
!
! 1833 2016-04-07 14:23:03Z raasch
! check of spectra quantities moved to spectra_mod
!
! 1829 2016-04-07 12:16:29Z maronga
! Bugfix: output of user defined data required reset of the string 'unit'
!
! 1826 2016-04-07 12:01:39Z maronga
! Moved checks for radiation model to the respective module.
! Moved checks for plant canopy model to the respective module.
! Bugfix for check of too large roughness_length
!
! 1824 2016-04-07 09:55:29Z gronemeier
! Check if roughness_length < dz/2. Moved location_message(finished) to the end.
!
! 1822 2016-04-07 07:49:42Z hoffmann
! PALM collision kernel deleted. Collision algorithms are checked for correct
! spelling.
!
! Tails removed.
! !
! Checks for use_sgs_for_particles adopted for the use of droplets with
! use_sgs_for_particles adopted.
!
! Unused variables removed.
!
! 1817 2016-04-06 15:44:20Z maronga
! Moved checks for land_surface model to the respective module
!
! 1806 2016-04-05 18:55:35Z gronemeier
! Check for recycling_yshift
!
! 1804 2016-04-05 16:30:18Z maronga
! Removed code for parameter file check (__check)
!
! 1795 2016-03-18 15:00:53Z raasch
! Bugfix: setting of initial scalar profile in ocean
!
! 1788 2016-03-10 11:01:04Z maronga
! Added check for use of most_method = 'lookup' in combination with water
! surface presribed in the land surface model. Added output of z0q.
! Syntax layout improved.
!
! 1786 2016-03-08 05:49:27Z raasch
! cpp-direktives for spectra removed, check of spectra level removed
!
! 1783 2016-03-06 18:36:17Z raasch
! netcdf variables and module name changed,
! check of netcdf precision removed (is done in the netcdf module)
!
! 1764 2016-02-28 12:45:19Z raasch
! output of nest id in run description header,
! bugfix: check of validity of lateral boundary conditions moved to parin
!
! 1762 2016-02-25 12:31:13Z hellstea
! Introduction of nested domain feature
!
! 1745 2016-02-05 13:06:51Z gronemeier
! Bugfix: check data output intervals to be /= 0.0 in case of parallel NetCDF4
!
! 1701 2015-11-02 07:43:04Z maronga
! Bugfix: definition of rad_net timeseries was missing
!
! 1691 2015-10-26 16:17:44Z maronga
! Added output of Obukhov length (ol) and radiative heating rates for RRTMG.
! Added checks for use of radiation / lsm with topography.
!
! 1682 2015-10-07 23:56:08Z knoop
! Code annotations made doxygen readable
!
! 1606 2015-06-29 10:43:37Z maronga
! Added check for use of RRTMG without netCDF.
!
! 1587 2015-05-04 14:19:01Z maronga
! Added check for set of albedos when using RRTMG
!
! 1585 2015-04-30 07:05:52Z maronga
! Added support for RRTMG
!
! 1575 2015-03-27 09:56:27Z raasch
! multigrid_fast added as allowed pressure solver
!
! 1573 2015-03-23 17:53:37Z suehring
! Bugfix: check for advection schemes in case of non-cyclic boundary conditions
!
! 1557 2015-03-05 16:43:04Z suehring
! Added checks for monotonic limiter
!
! 1555 2015-03-04 17:44:27Z maronga
! Added output of r_a and r_s. Renumbering of LSM PA-messages.
!
! 1553 2015-03-03 17:33:54Z maronga
! Removed check for missing soil temperature profile as default values were added.
!
! 1551 2015-03-03 14:18:16Z maronga
! Added various checks for land surface and radiation model. In the course of this
! action, the length of the variable var has to be increased
!
! 1504 2014-12-04 09:23:49Z maronga
! Bugfix: check for large_scale forcing got lost.
!
! 1500 2014-12-03 17:42:41Z maronga
! Boundary conditions changed to dirichlet for land surface model
!
! 1496 2014-12-02 17:25:50Z maronga
! Added checks for the land surface model
!
! 1484 2014-10-21 10:53:05Z kanani
! Changes due to new module structure of the plant canopy model:
!   module plant_canopy_model_mod added,
!   checks regarding usage of new method for leaf area density profile
!   construction added,
!   lad-profile construction moved to new subroutine init_plant_canopy within
!   the module plant_canopy_model_mod,
!   drag_coefficient renamed to canopy_drag_coeff.
! Missing KIND-attribute for REAL constant added
!
! 1455 2014-08-29 10:47:47Z heinze
! empty time records in volume, cross-section and masked data output prevented
! in case of non-parallel netcdf-output in restart runs
!
! 1429 2014-07-15 12:53:45Z knoop
! run_description_header exended to provide ensemble_member_nr if specified
!
! 1425 2014-07-05 10:57:53Z knoop
! bugfix: perturbation domain modified for parallel random number generator
!
! 1402 2014-05-09 14:25:13Z raasch
! location messages modified
!
! 1400 2014-05-09 14:03:54Z knoop
! Check random generator extended by option random-parallel
!
! 1384 2014-05-02 14:31:06Z raasch
! location messages added
!
! 1365 2014-04-22 15:03:56Z boeske
! Usage of large scale forcing for pt and q enabled:
! output for profiles of large scale advection (td_lsa_lpt, td_lsa_q),
! large scale subsidence (td_sub_lpt, td_sub_q)
! and nudging tendencies (td_nud_lpt, td_nud_q, td_nud_u and td_nud_v) added,
! check if use_subsidence_tendencies is used correctly
!
! 1361 2014-04-16 15:17:48Z hoffmann
! PA0363 removed
! PA0362 changed
!
! 1359 2014-04-11 17:15:14Z hoffmann
! Do not allow the execution of PALM with use_particle_tails, since particle
! tails are currently not supported by our new particle structure.
!
! PA0084 not necessary for new particle structure
!
! 1353 2014-04-08 15:21:23Z heinze
! REAL constants provided with KIND-attribute
!
! 1330 2014-03-24 17:29:32Z suehring
! In case of SGS-particle velocity advection of TKE is also allowed with
! dissipative 5th-order scheme.
!
! 1327 2014-03-21 11:00:16Z raasch
! "baroclinicity" renamed "baroclinity", "ocean version" replaced by "ocean mode"
! bugfix: duplicate error message 56 removed,
! check of data_output_format and do3d_compress removed
!
! 1322 2014-03-20 16:38:49Z raasch
! some REAL constants defined as wp-kind
!
! 1320 2014-03-20 08:40:49Z raasch
! Kind-parameters added to all INTEGER and REAL declaration statements,
! kinds are defined in new module kinds,
! revision history before 2012 removed,
! comment fields (!:) to be used for variable explanations added to
! all variable declaration statements
!
! 1308 2014-03-13 14:58:42Z fricke
! +netcdf_data_format_save
! Calculate fixed number of output time levels for parallel netcdf output.
! For masked data, parallel netcdf output is not tested so far, hence
! netcdf_data_format is switched back to non-paralell output.
!
! 1299 2014-03-06 13:15:21Z heinze
! enable usage of large_scale subsidence in combination with large_scale_forcing
! output for profile of large scale vertical velocity w_subs added
!
! 1276 2014-01-15 13:40:41Z heinze
! Use LSF_DATA also in case of Dirichlet bottom boundary condition for scalars
!
! 1241 2013-10-30 11:36:58Z heinze
! output for profiles of ug and vg added
! set surface_heatflux and surface_waterflux also in dependence on
! large_scale_forcing
! checks for nudging and large scale forcing from external file
!
! 1236 2013-09-27 07:21:13Z raasch
! check number of spectra levels
!
! 1216 2013-08-26 09:31:42Z raasch
! check for transpose_compute_overlap (temporary)
!
! 1214 2013-08-21 12:29:17Z kanani
! additional check for simultaneous use of vertical grid stretching
! and particle advection
!
! 1212 2013-08-15 08:46:27Z raasch
! checks for poisfft_hybrid removed
!
! 1210 2013-08-14 10:58:20Z raasch
! check for fftw
!
! 1179 2013-06-14 05:57:58Z raasch
! checks and settings of buoyancy parameters and switches revised,
! initial profile for rho_ocean added to hom (id=77)
!
! 1174 2013-05-31 10:28:08Z gryschka
! Bugfix in computing initial profiles for ug, vg, lad, q in case of Atmosphere
!
! 1159 2013-05-21 11:58:22Z fricke
! bc_lr/ns_dirneu/neudir removed
!
! 1115 2013-03-26 18:16:16Z hoffmann
! unused variables removed
! drizzle can be used without precipitation
!
! 1111 2013-03-08 23:54:10Z raasch
! ibc_p_b = 2 removed
!
! 1103 2013-02-20 02:15:53Z raasch
! Bugfix: turbulent inflow must not require cyclic fill in restart runs
!
! 1092 2013-02-02 11:24:22Z raasch
! unused variables removed
!
! 1069 2012-11-28 16:18:43Z maronga
! allow usage of topography in combination with cloud physics
!
! 1065 2012-11-22 17:42:36Z hoffmann
! Bugfix: It is not allowed to use cloud_scheme = seifert_beheng without
!         precipitation in order to save computational resources.
!
! 1060 2012-11-21 07:19:51Z raasch
! additional check for parameter turbulent_inflow
!
! 1053 2012-11-13 17:11:03Z hoffmann
! necessary changes for the new two-moment cloud physics scheme added:
! - check cloud physics scheme (Kessler or Seifert and Beheng)
! - plant_canopy is not allowed
! - currently, only cache loop_optimization is allowed
! - initial profiles of nr, qr
! - boundary condition of nr, qr
! - check output quantities (qr, nr, prr)
!
! 1036 2012-10-22 13:43:42Z raasch
! code put under GPL (PALM 3.9)
!
! 1031/1034 2012-10-22 11:32:49Z raasch
! check of netcdf4 parallel file support
!
! 1019 2012-09-28 06:46:45Z raasch
! non-optimized version of prognostic_equations not allowed any more
!
! 1015 2012-09-27 09:23:24Z raasch
! acc allowed for loop optimization,
! checks for adjustment of mixing length to the Prandtl mixing length removed
!
! 1003 2012-09-14 14:35:53Z raasch
! checks for cases with unequal subdomain sizes removed
!
! 1001 2012-09-13 14:08:46Z raasch
! all actions concerning leapfrog- and upstream-spline-scheme removed
!
! 996 2012-09-07 10:41:47Z raasch
! little reformatting

! 978 2012-08-09 08:28:32Z fricke
! setting of bc_lr/ns_dirneu/neudir
! outflow damping layer removed
! check for z0h*
! check for pt_damping_width
!
! 964 2012-07-26 09:14:24Z raasch
! check of old profil-parameters removed
!
! 940 2012-07-09 14:31:00Z raasch
! checks for parameter neutral
!
! 924 2012-06-06 07:44:41Z maronga
! Bugfix: preprocessor directives caused error during compilation
!
! 892 2012-05-02 13:51:44Z maronga
! Bugfix for parameter file check ( excluding __netcdf4 )
!
! 866 2012-03-28 06:44:41Z raasch
! use only 60% of the geostrophic wind as translation speed in case of Galilean
! transformation and use_ug_for_galilei_tr = .T. in order to mimimize the
! timestep
!
! 861 2012-03-26 14:18:34Z suehring
! Check for topography and ws-scheme removed.
! Check for loop_optimization = 'vector' and ws-scheme removed.
!
! 845 2012-03-07 10:23:05Z maronga
! Bugfix: exclude __netcdf4 directive part from namelist file check compilation
!
! 828 2012-02-21 12:00:36Z raasch
! check of collision_kernel extended
!
! 825 2012-02-19 03:03:44Z raasch
! check for collision_kernel and curvature_solution_effects
!
! 809 2012-01-30 13:32:58Z maronga
! Bugfix: replaced .AND. and .NOT. with && and ! in the preprocessor directives
!
! 807 2012-01-25 11:53:51Z maronga
! New cpp directive "__check" implemented which is used by check_namelist_files
!
! Revision 1.1  1997/08/26 06:29:23  raasch
! Initial revision
!
!
! Description:
! ------------
!> Check control parameters and deduce further quantities.
!------------------------------------------------------------------------------!
 SUBROUTINE check_parameters


    USE arrays_3d
    USE constants
    USE control_parameters
    USE grid_variables
    USE indices
    USE kinds

    USE netcdf_interface,                                                      &
        ONLY:  dopr_unit, do3d_unit, netcdf_data_format,            &
               netcdf_data_format_string

    USE pegrid

    USE profil_parameter
    USE statistics
    USE stokes_drift_mod,                                                      &
        ONLY: stokes_drift_check_parameters
    USE statistics
    USE transpose_indices
    USE turbulence_closure_mod,                                                &
        ONLY:  tcm_check_parameters

    IMPLICIT NONE

    CHARACTER (LEN=varnamelength)  ::  var           !< variable name
    CHARACTER (LEN=7)   ::  unit                     !< unit of variable
    CHARACTER (LEN=8)   ::  date                     !< current date string
    CHARACTER (LEN=10)  ::  time                     !< current time string
    CHARACTER (LEN=20)  ::  ensemble_string          !< string containing number of ensemble member
    CHARACTER (LEN=100) ::  action                   !< flag string

    INTEGER(iwp) ::  i                               !< loop index
    INTEGER(iwp) ::  ilen                            !< string length
    INTEGER(iwp) ::  j                               !< loop index
    INTEGER(iwp) ::  k                               !< loop index
    INTEGER(iwp) ::  kk                              !< loop index
    INTEGER(iwp) ::  netcdf_data_format_save         !< initial value of netcdf_data_format
    INTEGER(iwp) ::  position                        !< index position of string

    LOGICAL     ::  found                            !< flag, true if output variable is already marked for averaging

    REAL(wp)    ::  dum                              !< dummy variable
    REAL(wp)    ::  gradient                         !< local gradient
    REAL(wp)    ::  remote = 0.0_wp                  !< MPI id of remote processor
    REAL(wp)    ::  time_to_be_simulated_from_reference_point  !< time to be simulated from reference point


    CALL location_message( 'checking parameters', .FALSE. )

    !
!-- Check for overlap combinations, which are not realized yet
    IF ( transpose_compute_overlap )  THEN
       IF ( numprocs == 1 )  STOP '+++ transpose-compute-overlap not implemented for single PE runs'
    ENDIF
!
!-- User settings for restart times requires that "restart" has been given as
!-- file activation string. Otherwise, binary output would not be saved by
!-- palmrun.
    IF (  ( restart_time /= 9999999.9_wp  .OR.  dt_restart /= 9999999.9_wp )   &
         .AND.  .NOT. write_binary )  THEN
       WRITE( message_string, * ) 'manual restart settings requires file ',    &
                                  'activation string "restart"'
       CALL message( 'check_parameters', 'PA0001', 1, 2, 0, 6, 0 )
    ENDIF


!
!-- Generate the file header which is used as a header for most of PALM's
!-- output files
    CALL DATE_AND_TIME( date, time )
    run_date = date(7:8)//'-'//date(5:6)//'-'//date(3:4)
    run_time = time(1:2)//':'//time(3:4)//':'//time(5:6)
    ensemble_string = ''

    WRITE ( run_description_header,                                            &
            '(A,2X,A,2X,A,A,A,I2.2,2X,A,1X,A)' )                               &
          TRIM( version ), TRIM( revision ), 'run: ',                          &
          TRIM( run_identifier ), '.', runnr, run_date, run_time

!
!-- Check turbulence closure setup
    CALL tcm_check_parameters
!
!-- Timestep schemes: runge-kutta-3
    intermediate_timestep_count_max = 3

!
!-- Initializing actions must have been set by the user

    WRITE(message_string,*) 'initializing_actions = ',initializing_actions
    CALL location_message(adjustl(trim(message_string)),.TRUE.)

    IF ( TRIM( initializing_actions ) == '' )  THEN
       message_string = 'no value specified for initializing_actions'
       CALL message( 'check_parameters', 'PA0149', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( TRIM( initializing_actions ) /= 'read_restart_data'  .AND.            &
         TRIM( initializing_actions ) /= 'cyclic_fill' )  THEN
!
!--    No restart run: several initialising actions are possible
       action = initializing_actions
       DO  WHILE ( TRIM( action ) /= '' )
          position = INDEX( action, ' ' )
          SELECT CASE ( action(1:position-1) )

             CASE ( 'set_constant_profiles', 'set_1d-model_profiles',          &
                    'by_user', 'initialize_vortex', 			       &
                    'initialize_2D_bubble', 'initialize_3D_bubble', 'inifor' )
                action = action(position+1:)

             CASE DEFAULT
                message_string = 'initializing_action = "' //                  &
                                 TRIM( action ) // '" unknown or not allowed'
                CALL message( 'check_parameters', 'PA0030', 1, 2, 0, 6, 0 )

          END SELECT
       ENDDO
    ENDIF

    IF ( INDEX( initializing_actions, 'set_constant_profiles' ) /= 0  .AND.    &
         INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       message_string = 'initializing_actions = "set_constant_profiles"' //    &
                        ' and "set_1d-model_profiles" are not allowed ' //     &
                        'simultaneously'
       CALL message( 'check_parameters', 'PA0031', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( INDEX( initializing_actions, 'set_constant_profiles' ) /= 0  .AND.    &
         INDEX( initializing_actions, 'by_user' ) /= 0 )  THEN
       message_string = 'initializing_actions = "set_constant_profiles"' //    &
                        ' and "by_user" are not allowed simultaneously'
       CALL message( 'check_parameters', 'PA0032', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( INDEX( initializing_actions, 'by_user' ) /= 0  .AND.                  &
         INDEX( initializing_actions, 'set_1d-model_profiles' ) /= 0 )  THEN
       message_string = 'initializing_actions = "by_user" and ' //             &
                        '"set_1d-model_profiles" are not allowed simultaneously'
       CALL message( 'check_parameters', 'PA0033', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Check for Stokes drift settings if Stokes forcing in ocean mode is used
    IF ( stokes_force ) CALL stokes_drift_check_parameters
    !
!-- In case of no model continuation run, check initialising parameters and
!-- deduce further quantities
    IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN

      !
!--    Initial profiles for 1D and 3D model, respectively (u,v further below)
       pt_init = pt_surface
       sa_init = sa_surface
       !--

       !--    Let the initial wind profiles be the calculated ug/vg profiles or
!--    interpolate them from wind profile data (if given)
       IF ( u_profile(1) == 9999999.9_wp  .AND.  v_profile(1) == 9999999.9_wp )  THEN

          u_init = 0.0_wp
          v_init = 0.0_wp

       ELSEIF ( u_profile(1) == 0.0_wp  .AND.  v_profile(1) == 0.0_wp )  THEN

          IF ( uv_heights(1) /= 0.0_wp )  THEN
             message_string = 'uv_heights(1) must be 0.0'
             CALL message( 'check_parameters', 'PA0345', 1, 2, 0, 6, 0 )
          ENDIF

          IF ( omega /= 0.0_wp )  THEN
             message_string = 'Coriolis force must be switched off (by setting omega=0.0)' //  &
                              ' when prescribing the forcing by u_profile and v_profile'
             CALL message( 'check_parameters', 'PA0347', 1, 2, 0, 6, 0 )
          ENDIF

          use_prescribed_profile_data = .TRUE.

          kk = 1
          u_init(0) = 0.0_wp
          v_init(0) = 0.0_wp

          DO  k = 1, nz+1

             IF ( kk < 100 )  THEN
                DO  WHILE ( uv_heights(kk+1) <= zu(k) )
                   kk = kk + 1
                   IF ( kk == 100 )  EXIT
                ENDDO
             ENDIF

             IF ( kk < 100  .AND.  uv_heights(kk+1) /= 9999999.9_wp )  THEN
                u_init(k) = u_profile(kk) + ( zu(k) - uv_heights(kk) ) /       &
                                       ( uv_heights(kk+1) - uv_heights(kk) ) * &
                                       ( u_profile(kk+1) - u_profile(kk) )
                v_init(k) = v_profile(kk) + ( zu(k) - uv_heights(kk) ) /       &
                                       ( uv_heights(kk+1) - uv_heights(kk) ) * &
                                       ( v_profile(kk+1) - v_profile(kk) )
             ELSE
                u_init(k) = u_profile(kk)
                v_init(k) = v_profile(kk)
             ENDIF

          ENDDO

       ELSE

          message_string = 'u_profile(1) and v_profile(1) must be 0.0'
          CALL message( 'check_parameters', 'PA0346', 1, 2, 0, 6, 0 )

       ENDIF

!
!--    Compute initial temperature profile using the given temperature gradients
       CALL init_vertical_profiles( pt_vertical_gradient_level_ind,          &
                                    pt_vertical_gradient_level,              &
                                    pt_vertical_gradient, pt_init,           &
                                    pt_surface, bc_pt_t_val )

!
!--    If required, compute initial salinity profile using the given salinity
!--    gradients
       CALL init_vertical_profiles( sa_vertical_gradient_level_ind,          &
                                    sa_vertical_gradient_level,              &
                                    sa_vertical_gradient, sa_init,           &
                                    sa_surface, dum )

    ENDIF

    f  = 2.0_wp * omega * SIN( latitude / 180.0_wp * pi )
    fs = 2.0_wp * omega * COS( latitude / 180.0_wp * pi )

!
!-- Check time step and cfl_factor
    IF ( dt /= -1.0_wp )  THEN
       IF ( dt <= 0.0_wp )  THEN
          WRITE( message_string, * ) 'dt = ', dt , ' <= 0.0'
          CALL message( 'check_parameters', 'PA0044', 1, 2, 0, 6, 0 )
       ENDIF
       dt_3d = dt
       dt_fixed = .TRUE.
    ENDIF

    IF ( cfl_factor <= 0.0_wp  .OR.  cfl_factor > 1.0_wp )  THEN
       IF ( cfl_factor == -1.0_wp )  THEN
          cfl_factor = 0.9_wp
       ELSE
          WRITE( message_string, * ) 'cfl_factor = ', cfl_factor,              &
                 ' out of range &0.0 < cfl_factor <= 1.0 is required'
          CALL message( 'check_parameters', 'PA0045', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Store simulated time at begin
    simulated_time_at_begin = simulated_time

!
    time_since_reference_point = 0.0_wp

!
!-- Check boundary conditions and set internal variables:
!
!-- Bottom boundary condition for the turbulent Kinetic energy
    IF ( bc_e_b == 'neumann' )  THEN
       ibc_e_b = 1
    ELSEIF ( bc_e_b == '(u*)**2+neumann' )  THEN
       ibc_e_b = 2
    ELSE
       message_string = 'unknown boundary condition: bc_e_b = "' //            &
                        TRIM( bc_e_b ) // '"'
       CALL message( 'check_parameters', 'PA0058', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Boundary conditions for perturbation pressure
    IF ( bc_p_b == 'dirichlet' )  THEN
       ibc_p_b = 0
    ELSEIF ( bc_p_b == 'neumann' )  THEN
       ibc_p_b = 1
    ELSE
       message_string = 'unknown boundary condition: bc_p_b = "' //            &
                        TRIM( bc_p_b ) // '"'
       CALL message( 'check_parameters', 'PA0059', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( bc_p_t == 'dirichlet' )  THEN
       ibc_p_t = 0
!-- TO_DO: later set bc_p_t to neumann before, in case of nested domain
    ELSEIF ( bc_p_t == 'neumann' )  THEN
       ibc_p_t = 1
    ELSE
       message_string = 'unknown boundary condition: bc_p_t = "' //            &
                        TRIM( bc_p_t ) // '"'
       CALL message( 'check_parameters', 'PA0061', 1, 2, 0, 6, 0 )
    ENDIF

      IF ( bc_pt_b == 'dirichlet' )  THEN
          ibc_pt_b = 0
       ELSEIF ( bc_pt_b == 'neumann' )  THEN
          ibc_pt_b = 1
       ELSE
          message_string = 'unknown boundary condition: bc_pt_b = "' //        &
                           TRIM( bc_pt_b ) // '"'
          CALL message( 'check_parameters', 'PA0062', 1, 2, 0, 6, 0 )
       ENDIF

    IF ( bc_pt_t == 'dirichlet' )  THEN
       ibc_pt_t = 0
    ELSEIF ( bc_pt_t == 'neumann' )  THEN
       ibc_pt_t = 1
    ELSEIF ( bc_pt_t == 'initial_gradient' )  THEN
       ibc_pt_t = 2
    ELSE
       message_string = 'unknown boundary condition: bc_pt_t = "' //           &
                        TRIM( bc_pt_t ) // '"'
       CALL message( 'check_parameters', 'PA0063', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( top_heatflux == 9999999.9_wp )  constant_top_heatflux = .FALSE.

    IF ( top_momentumflux_u == 9999999.9_wp  .OR.                              &
           top_momentumflux_v == 9999999.9_wp )  THEN
       message_string = 'both, top_momentumflux_u AND top_momentumflux_v ' //  &
                        'must be set'
       CALL message( 'check_parameters', 'PA0064', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- A given temperature at the top implies Dirichlet boundary condition for
!-- temperature. In this case specification of a constant heat flux is
!-- forbidden.
    IF ( ibc_pt_t == 0  .AND.  constant_top_heatflux  .AND.                    &
         top_heatflux /= 0.0_wp )  THEN
       message_string = 'boundary_condition: bc_pt_t = "' // TRIM( bc_pt_t ) //&
                        '" is not allowed with constant_top_heatflux = .TRUE.'
       CALL message( 'check_parameters', 'PA0067', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Boundary conditions for salinity
       IF ( bc_sa_t == 'dirichlet' )  THEN
          ibc_sa_t = 0
       ELSEIF ( bc_sa_t == 'neumann' )  THEN
          ibc_sa_t = 1
       ELSE
          message_string = 'unknown boundary condition: bc_sa_t = "' //        &
                           TRIM( bc_sa_t ) // '"'
          CALL message( 'check_parameters', 'PA0068', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( bc_sa_b == 'dirichlet' )  THEN
          ibc_sa_b = 0
          CALL location_message('ib_sa_b assigned for dirichlet conditions',.TRUE.) !CB
       ELSEIF ( bc_sa_b == 'neumann' )  THEN
          ibc_sa_b = 1
          CALL location_message('ib_sa_b assigned for neumann conditions',.TRUE.) !CB
       ELSE
          message_string = 'unknown boundary condition: bc_sa_b = "' //        &
                           TRIM( bc_sa_b ) // '"'
          CALL message( 'check_parameters', 'PA0068', 1, 2, 0, 6, 0 )
       ENDIF

       IF ( ibc_sa_t == 1  .AND.  top_salinityflux == 9999999.9_wp )  THEN
          message_string = 'boundary condition: bc_sa_t = "' //                &
                           TRIM( bc_sa_t ) // '" requires to set ' //          &
                           'top_salinityflux'
          CALL message( 'check_parameters', 'PA0069', 1, 2, 0, 6, 0 )
       ENDIF
!
!--    A fixed salinity at the top implies Dirichlet boundary condition for
!--    salinity. In this case specification of a constant salinity flux is
!--    forbidden.
       IF ( top_salinityflux == 9999999.9_wp )  constant_top_salinityflux = .FALSE.
       IF ( ibc_sa_t == 0  .AND.  constant_top_salinityflux  .AND.             &
            top_salinityflux /= 0.0_wp )  THEN
          message_string = 'boundary condition: bc_sa_t = "' //                &
                           TRIM( bc_sa_t ) // '" is not allowed with ' //      &
                           'top_salinityflux /= 0.0'
          CALL message( 'check_parameters', 'PA0070', 1, 2, 0, 6, 0 )
       ENDIF
!
!-- Boundary conditions for horizontal components of wind speed
    IF ( bc_uv_b == 'dirichlet' )  THEN
       ibc_uv_b = 0
    ELSEIF ( bc_uv_b == 'neumann' )  THEN
       ibc_uv_b = 1
    ELSE
       message_string = 'unknown boundary condition: bc_uv_b = "' //           &
                        TRIM( bc_uv_b ) // '"'
       CALL message( 'check_parameters', 'PA0076', 1, 2, 0, 6, 0 )
    ENDIF
      IF ( bc_uv_t == 'dirichlet' .OR. bc_uv_t == 'dirichlet_0' )  THEN
          ibc_uv_t = 0
          IF ( bc_uv_t == 'dirichlet_0' )  THEN
!
!--          Velocities for the initial u,v-profiles are set zero at the top
!--          in case of dirichlet_0 conditions
             u_init(nzt+1)    = 0.0_wp
             v_init(nzt+1)    = 0.0_wp
          ENDIF
       ELSEIF ( bc_uv_t == 'neumann' )  THEN
          ibc_uv_t = 1
       ELSE
          message_string = 'unknown boundary condition: bc_uv_t = "' //        &
                           TRIM( bc_uv_t ) // '"'
          CALL message( 'check_parameters', 'PA0077', 1, 2, 0, 6, 0 )
       ENDIF

!
!-- Compute and check, respectively, the Rayleigh Damping parameter
    IF ( rayleigh_damping_factor == -1.0_wp )  THEN
       rayleigh_damping_factor = 0.0_wp
    ELSE
       IF ( rayleigh_damping_factor < 0.0_wp  .OR.                             &
            rayleigh_damping_factor > 1.0_wp )  THEN
          WRITE( message_string, * )  'rayleigh_damping_factor = ',            &
                              rayleigh_damping_factor, ' out of range [0.0,1.0]'
          CALL message( 'check_parameters', 'PA0078', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

    IF ( rayleigh_damping_height == -1.0_wp )  THEN
       rayleigh_damping_height = 0.66666666666_wp * zu(nzb)
    ELSE
       IF ( rayleigh_damping_height > 0.0_wp  .OR.                          &
            rayleigh_damping_height < zu(nzb) )  THEN
          WRITE( message_string, * )  'rayleigh_damping_height = ',         &
                rayleigh_damping_height, ' out of range [0.0,', zu(nzb), ']'
          CALL message( 'check_parameters', 'PA0079', 1, 2, 0, 6, 0 )
       ENDIF
    ENDIF

!
!-- Set the default intervals for data output, if necessary
!-- NOTE: dt_dosp has already been set in spectra_parin
    IF ( dt_data_output /= 9999999.9_wp )  THEN
       IF ( dt_dopr           == 9999999.9_wp )  dt_dopr           = dt_data_output
       IF ( dt_do3d           == 9999999.9_wp )  dt_do3d           = dt_data_output
       IF ( dt_data_output_av == 9999999.9_wp )  dt_data_output_av = dt_data_output
    ENDIF

!
!-- Set the default skip time intervals for data output, if necessary
    IF ( skip_time_dopr    == 9999999.9_wp )                                   &
                                       skip_time_dopr    = skip_time_data_output
    IF ( skip_time_do3d    == 9999999.9_wp )                                   &
                                       skip_time_do3d    = skip_time_data_output
    IF ( skip_time_data_output_av == 9999999.9_wp )                            &
                                skip_time_data_output_av = skip_time_data_output

!
!-- Check the average intervals (first for 3d-data, then for profiles)
    IF ( averaging_interval > dt_data_output_av )  THEN
       WRITE( message_string, * )  'averaging_interval = ',                    &
             averaging_interval, ' must be <= dt_data_output_av = ',           &
             dt_data_output_av
       CALL message( 'check_parameters', 'PA0085', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( averaging_interval_pr == 9999999.9_wp )  THEN
       averaging_interval_pr = averaging_interval
    ENDIF

    IF ( averaging_interval_pr > dt_dopr )  THEN
       WRITE( message_string, * )  'averaging_interval_pr = ',                 &
             averaging_interval_pr, ' must be <= dt_dopr = ', dt_dopr
       CALL message( 'check_parameters', 'PA0086', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Set the default interval for profiles entering the temporal average
    IF ( dt_averaging_input_pr == 9999999.9_wp )  THEN
       dt_averaging_input_pr = dt_averaging_input
    ENDIF

!
!-- Check the sample rate for averaging (first for 3d-data, then for profiles)
    IF ( dt_averaging_input > averaging_interval )  THEN
       WRITE( message_string, * )  'dt_averaging_input = ',                    &
                dt_averaging_input, ' must be <= averaging_interval = ',       &
                averaging_interval
       CALL message( 'check_parameters', 'PA0088', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( dt_averaging_input_pr > averaging_interval_pr )  THEN
       WRITE( message_string, * )  'dt_averaging_input_pr = ',                 &
                dt_averaging_input_pr, ' must be <= averaging_interval_pr = ', &
                averaging_interval_pr
       CALL message( 'check_parameters', 'PA0089', 1, 2, 0, 6, 0 )
    ENDIF
!
!-- Determine the number of output profiles and check whether they are
!-- permissible
    DO  WHILE ( data_output_pr(dopr_n+1) /= '          ' )

       dopr_n = dopr_n + 1
       i = dopr_n

!
!--    Determine internal profile number (for hom, homs)
!--    and store height levels
       SELECT CASE ( TRIM( data_output_pr(i) ) )

          CASE ( 'u' )
             dopr_index(i) = 1
             dopr_unit(i)  = 'm/s'
             hom(:,2,1)  = zu

          CASE ( 'v' )
             dopr_index(i) = 2
             dopr_unit(i)  = 'm/s'
             hom(:,2,2)  = zu

          CASE ( 'w' )
             dopr_index(i) = 3
             dopr_unit(i)  = 'm/s'
             hom(:,2,3)  = zw

          CASE ( 'pt' )
             dopr_index(i) = 4
             dopr_unit(i)  = 'K'
             hom(:,2,4)  = zu

          CASE ( 'sa' )
             dopr_index(i) = 5
             dopr_unit(i)  = 'psu'
             hom(:,2,5) = zu

          CASE ( 'p' )
             dopr_index(i) = 6
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,6) = zu

          CASE ( 'e' )
             dopr_index(i)  = 7
             dopr_unit(i)   = 'm2/s2'
             hom(:,2,7)   = zu

          CASE ( 'km' )
             dopr_index(i)  = 8
             dopr_unit(i)   = 'm2/s'
             hom(:,2,8)   = zu

          CASE ( 'kh' )
             dopr_index(i)   = 9
             dopr_unit(i)    = 'm2/s'
             hom(:,2,9)   = zu

          CASE ( 'w"u"' )
             dopr_index(i) = 10
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,10) = zw

          CASE ( 'w*u*' )
             dopr_index(i) = 11
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,11) = zw

          CASE ( 'w"v"' )
             dopr_index(i) = 12
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,12) = zw

          CASE ( 'w*v*' )
             dopr_index(i) = 13
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,13) = zw

          CASE ( 'w"pt"' )
             dopr_index(i) = 14
             dopr_unit(i)  = 'K m/s'
             hom(:,2,14) = zw

          CASE ( 'w*pt*' )
             dopr_index(i) = 15
             dopr_unit(i)  = 'K m/s'
             hom(:,2,15) = zw

          CASE ( 'w"sa"' )
             dopr_index(i) = 16
             dopr_unit(i)  = 'psu m/s'
             hom(:,2,16) = zw

          CASE ( 'w*sa*' )
             dopr_index(i) = 17
             dopr_unit(i)  = 'psu m/s'
             hom(:,2,17) = zw

          CASE ( 'u*2' )
             dopr_index(i) = 18
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,18) = zu

          CASE ( 'v*2' )
             dopr_index(i) = 19
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,19) = zu

          CASE ( 'w*2' )
             dopr_index(i) = 20
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,20) = zw

          CASE ( 'pt*2' )
             dopr_index(i) = 21
             dopr_unit(i)  = 'K2'
             hom(:,2,21) = zu

          CASE ( 'sa*2' )
             dopr_index(i) = 22
             dopr_unit(i)  = 'psu^2'
             hom(:,2,22) = zu

          CASE ( 'w*3' )
             dopr_index(i) = 23
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,23) = zw

          CASE ( 'w*p*' )
             dopr_index(i) = 24
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,24) = zu

          CASE ( 'w*e*' )
             dopr_index(i) = 25
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,25) = zu

          CASE ( 'rho_ocean' )
             dopr_index(i) = 26
             dopr_unit(i)  = 'kg/m3'
             hom(:,2,26) = zu

          CASE ( 'alpha_T' )
             dopr_index(i) = 27
             dopr_unit(i)  = '^oC^{-1}'
             hom(:,2,27) = zu

          CASE ( 'beta_S' )
             dopr_index(i) = 28
             dopr_unit(i)  = 'PSU^{-1}'
             hom(:,2,28) = zu

          CASE ( 'solar3d' )
             dopr_index(i) = 29
             dopr_unit(i)  = '^oC^{-1}/s'
             hom(:,2,29) = zu

          CASE ( 'prho' )
             dopr_index(i) = 30
             dopr_unit(i)  = 'kg/m3'
             hom(:,2,30) = zu

          CASE ( 'u_stk' )
             IF (  .NOT.  stokes_force ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for stokes_force = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 31
                dopr_unit(i)  = 'm/s'
                hom(:,2,31) = zu
             ENDIF

          CASE ( 'v_stk' )
             IF (  .NOT.  stokes_force ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for stokes_force = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 32
                dopr_unit(i)  = 'm/s'
                hom(:,2,32) = zu
             ENDIF

          CASE ( 'u_stk_zw' )
             IF (  .NOT.  stokes_force ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for stokes_force = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 33
                dopr_unit(i)  = 'm/s'
                hom(:,2,33) = zw
             ENDIF

          CASE ( 'v_stk_zw' )
             IF (  .NOT.  stokes_force ) THEN
                message_string = 'data_output_pr = ' //                        &
                                 TRIM( data_output_pr(i) ) // ' is not imp' // &
                                 'lemented for stokes_force = .FALSE.'
                CALL message( 'check_parameters', 'PA0091', 1, 2, 0, 6, 0 )
             ELSE
                dopr_index(i) = 34
                dopr_unit(i)  = 'm/s'
                hom(:,2,34) = zw
             ENDIF

          CASE ( 'u*v*' )
             dopr_index(i) = 35
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,35) = zu

          CASE ( 'v*u*' )
             dopr_index(i) = 36
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,36) = zu

          CASE ( 'eps' )
             dopr_index(i) = 37
             dopr_unit(i)  = 'm2/s3'
             hom(:,2,37) = zu

          CASE ( 'sgs_tran' )
             dopr_index(i) = 38
             dopr_unit(i)  = 'm2/s3'
             hom(:,2,38) = zu

          CASE ( 'w"u"u*' )
             dopr_index(i) = 39
             dopr_unit(i)  = 'm3/s3'
             hom(:,2,39) = zu

          CASE ( 'diss' )
             dopr_index(i) = 40
             dopr_unit(i)  = 'm2/s3'
             hom(:,2,40) = zw

          CASE ( 'e*' )
             dopr_index(i) = 41
             dopr_unit(i)  = 'm2/s2'
             hom(:,2,41) = zw

          CASE DEFAULT
             IF ( unit == 'illegal' )  THEN
                IF ( data_output_pr_user(1) /= ' ' )  THEN
                   message_string = 'illegal value for data_output_pr or ' //  &
                                    'data_output_pr_user = "' //               &
                                    TRIM( data_output_pr(i) ) // '"'
                   CALL message( 'check_parameters', 'PA0097', 1, 2, 0, 6, 0 )
                ELSE
                   message_string = 'illegal value for data_output_pr = "' //  &
                                    TRIM( data_output_pr(i) ) // '"'
                   CALL message( 'check_parameters', 'PA0098', 1, 2, 0, 6, 0 )
                ENDIF
             ENDIF

       END SELECT

    ENDDO
!
!-- Check and set steering parameters for 2d/3d data output and averaging
    i   = 1
    DO  WHILE ( data_output(i) /= ' '  .AND.  i <= 500 )
!
!--    Check for data averaging
       ilen = LEN_TRIM( data_output(i) )
       j = 0                                                 ! no data averaging
       IF ( ilen > 3 )  THEN
          IF ( data_output(i)(ilen-2:ilen) == '_av' )  THEN
             j = 1                                           ! data averaging
             data_output(i) = data_output(i)(1:ilen-3)
          ENDIF
       ENDIF
!
!--    Check for cross section or volume data
       ilen = LEN_TRIM( data_output(i) )
       k = 0                                                   ! 3d data
       var = data_output(i)(1:ilen)

!
!--    Check for allowed value and set units
       SELECT CASE ( TRIM( var ) )

          CASE ( 'solar3d' )
             unit = '^oC s^{-1}'

          CASE ( 'rho_ocean' )
             unit = 'kg/m3'

          CASE ( 'alpha_T' )
             unit = '^oC^{-1}'

          CASE ( 'beta_S' )
             unit = 'psu^{-1}'

          CASE ( 'sa' )
             unit = 'psu'

          CASE ( 'e', 'p', 'pt', 'u', 'v', 'w' )
             IF ( TRIM( var ) == 'e'  )  unit = 'm2/s2'
             IF ( TRIM( var ) == 'p'  )  unit = 'Pa'
             IF ( TRIM( var ) == 'pt' )  unit = 'K'
             IF ( TRIM( var ) == 'u'  )  unit = 'm/s'
             IF ( TRIM( var ) == 'v'  )  unit = 'm/s'
             IF ( TRIM( var ) == 'w'  )  unit = 'm/s'
             CONTINUE

          CASE ( 'shf_sol*', 'shf*')
             IF ( k == 0  )  THEN
                message_string = 'illegal value for data_output: "' //         &
                                 TRIM( var ) // '" & only 2d-horizontal ' //   &
                                 'cross sections are allowed for this value'
                CALL message( 'check_parameters', 'PA0111', 1, 2, 0, 6, 0 )
             ENDIF
             IF ( TRIM( var ) == 'shf*'   )  unit = 'K*m/s'
             IF ( TRIM( var ) == 'shf_sol*' ) unit = 'K*m/s'
          CASE DEFAULT

             CONTINUE

       END SELECT
!
!--    Set the internal steering parameters appropriately
       IF ( k == 0 )  THEN
          do3d_no(j)              = do3d_no(j) + 1
          do3d(j,do3d_no(j))      = data_output(i)
          do3d_unit(j,do3d_no(j)) = unit
       ENDIF

       IF ( j == 1 )  THEN
!
!--       Check, if variable is already subject to averaging
          found = .FALSE.
          DO  k = 1, doav_n
             IF ( TRIM( doav(k) ) == TRIM( var ) )  found = .TRUE.
          ENDDO

          IF ( .NOT. found )  THEN
             doav_n = doav_n + 1
             doav(doav_n) = var
          ENDIF
       ENDIF

       i = i + 1
    ENDDO

!
!-- Averaged 2d or 3d output requires that an averaging interval has been set
    IF ( doav_n > 0  .AND.  averaging_interval == 0.0_wp )  THEN
       WRITE( message_string, * )  'output of averaged quantity "',            &
                                   TRIM( doav(1) ), '_av" requires to set a ', &
                                   'non-zero averaging interval'
       CALL message( 'check_parameters', 'PA0323', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Upper plot limit for 3D arrays
    IF ( nz_do3d == -9999 )  nz_do3d = nzt + 1

!
!-- Set output format string (used in header)
    SELECT CASE ( netcdf_data_format )
       CASE ( 1 )
          netcdf_data_format_string = 'netCDF classic'
       CASE ( 2 )
          netcdf_data_format_string = 'netCDF 64bit offset'
       CASE ( 3 )
          netcdf_data_format_string = 'netCDF4/HDF5'
       CASE ( 4 )
          netcdf_data_format_string = 'netCDF4/HDF5 classic'
       CASE ( 5 )
          netcdf_data_format_string = 'parallel netCDF4/HDF5'
       CASE ( 6 )
          netcdf_data_format_string = 'parallel netCDF4/HDF5 classic'

    END SELECT

!
!-- Check the NetCDF data format
    IF ( netcdf_data_format > 2 )  THEN
#if defined( __netcdf4 )
       CONTINUE
#else
       message_string = 'netCDF: netCDF4 format requested but no ' //          &
                        'cpp-directive __netcdf4 given & switch '  //          &
                        'back to 64-bit offset format'
       CALL message( 'check_parameters', 'PA0171', 0, 1, 0, 6, 0 )
       netcdf_data_format = 2
#endif
    ENDIF
    IF ( netcdf_data_format > 4 )  THEN
#if defined( __netcdf4 ) && defined( __netcdf4_parallel )
       CONTINUE
#else
       message_string = 'netCDF: netCDF4 parallel output requested but no ' // &
                        'cpp-directive __netcdf4_parallel given & switch '   //&
                        'back to netCDF4 non-parallel output'
       CALL message( 'check_parameters', 'PA0099', 0, 1, 0, 6, 0 )
       netcdf_data_format = netcdf_data_format - 2
#endif
    ENDIF

!
!-- Calculate fixed number of output time levels for parallel netcdf output.
!-- The time dimension has to be defined as limited for parallel output,
!-- because otherwise the I/O performance drops significantly.
    IF ( netcdf_data_format > 4 )  THEN

!
!--    Check if any of the follwoing data output interval is 0.0s, which is
!--    not allowed for parallel output.
       CALL check_dt_do( dt_do3d,           'dt_do3d'           )
       CALL check_dt_do( dt_data_output_av, 'dt_data_output_av' )

!--    Set needed time levels (ntdim) to
!--    saved time levels + to be saved time levels.
       ntdim_3d(0) = do3d_time_count(0) + CEILING(                             &
                     ( end_time - MAX( skip_time_do3d,                         &
                                       simulated_time_at_begin )               &
                     ) / dt_do3d )
       IF ( do3d_at_begin ) ntdim_3d(0) = ntdim_3d(0) + 1

       ntdim_3d(1) = do3d_time_count(1) + CEILING(                             &
                     ( end_time - MAX( skip_time_data_output_av,               &
                                       simulated_time_at_begin )               &
                     ) / dt_data_output_av )
    ENDIF

!
!-- Check random generator
    IF ( (random_generator /= 'system-specific'      .AND.                     &
          random_generator /= 'random-parallel'   )  .AND.                     &
          random_generator /= 'numerical-recipes' )  THEN
       message_string = 'unknown random generator: random_generator = "' //    &
                        TRIM( random_generator ) // '"'
       CALL message( 'check_parameters', 'PA0135', 1, 2, 0, 6, 0 )
    ENDIF

!
!-- Determine upper and lower hight level indices for random perturbations
    IF ( disturbance_level_b == -9999999.9_wp )  THEN
       disturbance_level_b     = zu((nzt*2)/3)
       disturbance_level_ind_b = ( nzt * 2 ) / 3
    ELSEIF ( disturbance_level_b < zu(3) )  THEN
       WRITE( message_string, * )  'disturbance_level_b = ',                   &
                           disturbance_level_b, ' must be >= ', zu(3), '(zu(3))'
       CALL message( 'check_parameters', 'PA0126', 1, 2, 0, 6, 0 )
    ELSEIF ( disturbance_level_b > zu(nzt-2) )  THEN
       WRITE( message_string, * )  'disturbance_level_b = ',                   &
                   disturbance_level_b, ' must be <= ', zu(nzt-2), '(zu(nzt-2))'
       CALL message( 'check_parameters', 'PA0127', 1, 2, 0, 6, 0 )
    ELSE
       DO  k = 3, nzt-2
          IF ( disturbance_level_b <= zu(k) )  THEN
             disturbance_level_ind_b = k
             EXIT
          ENDIF
       ENDDO
    ENDIF

    IF ( disturbance_level_t == -9999999.9_wp )  THEN
       disturbance_level_t     = zu(nzt-3)
       disturbance_level_ind_t = nzt - 3
    ELSEIF ( disturbance_level_t > zu(nzt-2) )  THEN
       WRITE( message_string, * )  'disturbance_level_t = ',                   &
                   disturbance_level_t, ' must be <= ', zu(nzt-2), '(zu(nzt-2))'
       CALL message( 'check_parameters', 'PA0128', 1, 2, 0, 6, 0 )
    ELSEIF ( disturbance_level_t < disturbance_level_b )  THEN
       WRITE( message_string, * )  'disturbance_level_t = ',                   &
                   disturbance_level_t, ' must be >= disturbance_level_b = ',  &
                   disturbance_level_b
       CALL message( 'check_parameters', 'PA0129', 1, 2, 0, 6, 0 )
    ELSE
       DO  k = 3, nzt-2
          IF ( disturbance_level_t <= zu(k) )  THEN
             disturbance_level_ind_t = k
             EXIT
          ENDIF
       ENDDO
    ENDIF

!
!-- Check again whether the levels determined this way are ok.
!-- Error may occur at automatic determination and too few grid points in
!-- z-direction.
    IF ( disturbance_level_ind_t < disturbance_level_ind_b )  THEN
       WRITE( message_string, * )  'disturbance_level_ind_t = ',               &
                disturbance_level_ind_t, ' must be >= ',                       &
                'disturbance_level_ind_b = ', disturbance_level_ind_b
       CALL message( 'check_parameters', 'PA0130', 1, 2, 0, 6, 0 )
    ENDIF

    IF ( random_generator == 'random-parallel' )  THEN
       dist_nxl = nxl;  dist_nxr = nxr
       dist_nys = nys;  dist_nyn = nyn
    ELSE
       dist_nxl = 0;  dist_nxr = nx
       dist_nys = 0;  dist_nyn = ny
    ENDIF
!
!-- Set time for the next user defined restart (time_restart is the
!-- internal parameter for steering restart events)
    IF ( restart_time /= 9999999.9_wp )  THEN
       IF ( restart_time > time_since_reference_point )  THEN
          time_restart = restart_time
       ENDIF
    ELSE
!
!--    In case of a restart run, set internal parameter to default (no restart)
!--    if the NAMELIST-parameter restart_time is omitted
       time_restart = 9999999.9_wp
    ENDIF

!
!-- Check pressure gradient conditions
    IF ( dp_external )  THEN
       IF ( dp_level_b < zu(nzb)  .OR.  dp_level_b > zu(nzt) )  THEN
          WRITE( message_string, * )  'dp_level_b = ', dp_level_b, ' is out ', &
               ' of range [zu(nzb), zu(nzt)]'
          CALL message( 'check_parameters', 'PA0151', 1, 2, 0, 6, 0 )
       ENDIF
       IF ( .NOT. ANY( dpdxy /= 0.0_wp ) )  THEN
          WRITE( message_string, * )  'dp_external is .TRUE. but dpdxy is ze', &
               'ro, i.e. the external pressure gradient will not be applied'
          CALL message( 'check_parameters', 'PA0152', 0, 1, 0, 6, 0 )
       ENDIF
    ENDIF
    IF ( ANY( dpdxy /= 0.0_wp )  .AND.  .NOT.  dp_external )  THEN
       WRITE( message_string, * )  'dpdxy is nonzero but dp_external is ',     &
            '.FALSE., i.e. the external pressure gradient & will not be applied'
       CALL message( 'check_parameters', 'PA0153', 0, 1, 0, 6, 0 )
    ENDIF
!
!-- Prevent empty time records in volume, cross-section and masked data in case
!-- of non-parallel netcdf-output in restart runs
    IF ( netcdf_data_format < 5 )  THEN
       IF ( TRIM( initializing_actions ) == 'read_restart_data' )  THEN
          do3d_time_count    = 0
       ENDIF
    ENDIF

!
!-- Check for valid setting of most_method
    IF ( TRIM( most_method ) /= 'circular'  .AND.                              &
         TRIM( most_method ) /= 'newton'    .AND.                              &
         TRIM( most_method ) /= 'lookup' )  THEN
       message_string = 'most_method = "' // TRIM( most_method ) //            &
                        '" is unknown'
       CALL message( 'check_parameters', 'PA0416', 1, 2, 0, 6, 0 )
    ENDIF

    CALL location_message( 'finished', .TRUE. )
!

 CONTAINS

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Check the length of data output intervals. In case of parallel NetCDF output
!> the time levels of the output files need to be fixed. Therefore setting the
!> output interval to 0.0s (usually used to output each timestep) is not
!> possible as long as a non-fixed timestep is used.
!------------------------------------------------------------------------------!

    SUBROUTINE check_dt_do( dt_do, dt_do_name )

       IMPLICIT NONE

       CHARACTER (LEN=*), INTENT (IN) :: dt_do_name !< parin variable name

       REAL(wp), INTENT (INOUT)       :: dt_do      !< data output interval

       IF ( dt_do == 0.0_wp )  THEN
          IF ( dt_fixed )  THEN
             WRITE( message_string, '(A,F9.4,A)' )  'Output at every '  //     &
                    'timestep is wanted (' // dt_do_name // ' = 0.0).&'//      &
                    'The output interval is set to the fixed timestep dt '//   &
                    '= ', dt, 's.'
             CALL message( 'check_parameters', 'PA0060', 0, 0, 0, 6, 0 )
             dt_do = dt
          ELSE
             message_string = dt_do_name // ' = 0.0 while using a ' //         &
                              'variable timestep and parallel netCDF4 ' //     &
                              'is not allowed.'
             CALL message( 'check_parameters', 'PA0081', 1, 2, 0, 6, 0 )
          ENDIF
       ENDIF

    END SUBROUTINE check_dt_do




!------------------------------------------------------------------------------!
! Description:
! ------------
!> Inititalizes the vertical profiles of scalar quantities.
!------------------------------------------------------------------------------!

    SUBROUTINE init_vertical_profiles( vertical_gradient_level_ind,            &
                                       vertical_gradient_level,                &
                                       vertical_gradient,                      &
                                       pr_init, surface_value, bc_t_val )


       IMPLICIT NONE

       INTEGER(iwp) ::  i     !< counter
       INTEGER(iwp), DIMENSION(1:10) ::  vertical_gradient_level_ind !< vertical grid indices for gradient levels

       REAL(wp)     ::  bc_t_val      !< model top gradient
       REAL(wp)     ::  gradient      !< vertica gradient of the respective quantity
       REAL(wp)     ::  surface_value !< surface value of the respecitve quantity


       REAL(wp), DIMENSION(0:nz+1) ::  pr_init                 !< initialisation profile
       REAL(wp), DIMENSION(1:10)   ::  vertical_gradient       !< given vertical gradient
       REAL(wp), DIMENSION(1:10)   ::  vertical_gradient_level !< given vertical gradient level

       i = 1
       gradient = 0.0_wp

       vertical_gradient_level_ind(1) = nzt+1
       DO  k = nzt, 0, -1
          IF ( i < 11 )  THEN
             IF ( vertical_gradient_level(i) >= zu(k)  .AND.            &
                  vertical_gradient_level(i) <= 0.0_wp )  THEN
                gradient = vertical_gradient(i) / 100.0_wp
                vertical_gradient_level_ind(i) = k + 1
                i = i + 1
             ENDIF
          ENDIF
          IF ( gradient /= 0.0_wp )  THEN
             IF ( k /= nzt )  THEN
                pr_init(k) = pr_init(k+1) - dzu(k+1) * gradient
             ELSE
                pr_init(k)   = surface_value - 0.5_wp * dzu(k+1) * gradient
                pr_init(k+1) = surface_value + 0.5_wp * dzu(k+1) * gradient
             ENDIF
          ELSE
             pr_init(k) = pr_init(k+1)
          ENDIF
!
!--       Avoid negative humidities
          IF ( pr_init(k) < 0.0_wp )  THEN
             pr_init(k) = 0.0_wp
          ENDIF
       ENDDO


!
!--    In case of no given gradients, choose zero gradient conditions
       IF ( vertical_gradient_level(1) == -999999.9_wp )  THEN
          vertical_gradient_level(1) = 0.0_wp
       ENDIF
!
!--    Store gradient at the top boundary for possible Neumann boundary condition
       bc_t_val  = ( pr_init(nzt+1) - pr_init(nzt) ) / dzu(nzt+1)

    END SUBROUTINE init_vertical_profiles

 END SUBROUTINE check_parameters
