#include "GAMER.h"

#ifdef SUPPORT_LIBYT




//-------------------------------------------------------------------------------------------------------
// Function    :  YT_SetParameter
// Description :  Set YT-specific parameters for the inline analysis
//
// Note        :  1. This function must be called in advance **every time** we invoke the inline analysis
//                2. Invoked by YT_Inline()
//
// Parameter   :  NPatchAllLv : Total number of patches at all levels
//                NField      : Total number of fields
//                NPatchLocal : Number of local patches at all levels
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void YT_SetParameter( const int NPatchAllLv, const int NField, const int NPatchLocalLv)
{

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// 1. prepare the simulation information for libyt
   yt_param_yt param_yt;

   param_yt.frontend                = "gamer";           // simulation frontend
// param_yt.fig_basename            = "fig_basename";    // figure base name (default=Fig%09d)

   param_yt.length_unit             = UNIT_L;            // units are in cgs
   param_yt.mass_unit               = UNIT_M;
   param_yt.time_unit               = UNIT_T;

#ifdef MHD
   param_yt.magnetic_unit           = UNIT_B;
#endif

   param_yt.current_time            = Time[0];
   param_yt.dimensionality          = 3;
   param_yt.refine_by               = 2;
   param_yt.num_grids               = NPatchAllLv;

   param_yt.num_fields              = NField;
   param_yt.num_grids_local         = NPatchLocalLv;

#  ifdef FLOAT8
   param_yt.field_ftype             = YT_DOUBLE;
#  else
   param_yt.field_ftype             = YT_FLOAT;
#  endif

   for (int d=0; d<3; d++)
   {
      param_yt.domain_dimensions[d] = NX0_TOT[d];
      param_yt.domain_left_edge [d] = 0.0;
      param_yt.domain_right_edge[d] = amr->BoxSize[d];
      param_yt.periodicity      [d] = ( OPT__BC_FLU[0] == BC_FLU_PERIODIC ) ? 1 : 0;
   }

#  ifdef COMOVING
   param_yt.cosmological_simulation = 1;
   param_yt.current_redshift        = 1.0/Time[0] - 1.0;
   param_yt.omega_lambda            = 1.0 - OMEGA_M0;
   param_yt.omega_matter            = OMEGA_M0;
   param_yt.hubble_constant         = HUBBLE0;
#  else
   param_yt.cosmological_simulation = 0;
   param_yt.current_redshift        = 0.0;
   param_yt.omega_lambda            = 0.0;
   param_yt.omega_matter            = 0.0;
   param_yt.hubble_constant         = 0.0;
#  endif


// 2. transfer simulation information to libyt
   if ( yt_set_parameter( &param_yt ) != YT_SUCCESS )    Aux_Error( ERROR_INFO, "yt_set_parameter() failed !!\n" );

#  ifdef MHD
   const int mhd = 1;
   if (yt_add_user_parameter_int("mhd", 1, &mhd) != YT_SUCCESS)  Aux_Error( ERROR_INFO, "yt_add_user_parameter() add mhd failed !!\n" );
#  endif

#  if ( MODEL == HYDRO )
   const double gamma = (double) GAMMA;
   const double mu = (double) MOLECULAR_WEIGHT;
   if (yt_add_user_parameter_double("gamma", 1, &gamma) != YT_SUCCESS )  Aux_Error( ERROR_INFO, "yt_add_user_parameter() add GAMMA failed !!\n" );
   if (yt_add_user_parameter_double("mu", 1, &mu) != YT_SUCCESS )  Aux_Error( ERROR_INFO, "yt_add_user_parameter() add MOLECULAR_WEIGHT failed !!\n" );
#  endif

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );


} // FUNCTION : YT_SetParameter



#endif // #ifdef SUPPORT_LIBYT
