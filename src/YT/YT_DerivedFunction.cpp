#include "GAMER.h"

#ifdef MHD
#ifdef SUPPORT_LIBYT

//-------------------------------------------------------------------------------------------------------
// Function    :  MagX/Y/Z_DerivedFunc
// Description :  Derived function for MagX/Y/Z
//
// Note        :  1. This function's pointer will be passed into libyt.
//                2. The argument should be declared like this, in order to match the libyt API.
//
// Parameter   :  GID            : Grid GID
//                Converted_MagX : Store the converted field data here.
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void MagX_DerivedFunc(long GID, double *Converted_MagX){
    // Get the dimension of the grid, and the data array pointer of the grid
    int Dimensions[3];
    void *DataRaw;
    yt_getGridInfo_Dimensions( gid, &Dimensions );
    yt_getGridInfo_FieldData( gid, "MagX", &DataRaw );

    // Cast the DataRaw
    real *Data = (real *) DataRaw;

    // Compute converted field via data. [z, y, x] direction.
    for (int k=0; k<Dimensions[0]; k++){
        for (int j=0; j<Dimensions[1]; j++){
            for (int i=0; i<Dimensions[2]; i++){
                int idx_Bx = IDX321_BX(i, j, k, Dimensions[2], Dimensions[1]);
                Converted_MagX[inx_Bx] = (double) 0.5 * ( Data[ idx_Bx ] + Data[ idx_Bx + 1 ] );
            }
        }
    }
}

void MagY_DerivedFunc(long GID, double *Converted_MagY){
    // Get the dimension of the grid, and the data array pointer of the grid
    int Dimensions[3];
    void *DataRaw;
    yt_getGridInfo_Dimensions( gid, &Dimensions );
    yt_getGridInfo_FieldData( gid, "MagY", &DataRaw );

    // Cast the DataRaw
    real *Data = (real *) DataRaw;

    // Compute converted field via data. [z, y, x] direction.
    for (int k=0; k<Dimensions[0]; k++){
        for (int j=0; j<Dimensions[1]; j++){
            for (int i=0; i<Dimensions[2]; i++){
                int idx_By = IDX321_BY(i, j, k, Dimensions[2], Dimensions[1]);
                Converted_MagY[inx_By] = (double) 0.5 * ( Data[ idx_By ] + Data[ idx_By + Dimensions[2] ] );
            }
        }
    }
}

void MagZ_DerivedFunc(long GID, double *Converted_MagZ){
    // Get the dimension of the grid, and the data array pointer of the grid
    int Dimensions[3];
    void *DataRaw;
    yt_getGridInfo_Dimensions( gid, &Dimensions );
    yt_getGridInfo_FieldData( gid, "MagZ", &DataRaw );

    // Cast the DataRaw
    real *Data = (real *) DataRaw;

    // Compute converted field via data. [z, y, x] direction.
    for (int k=0; k<Dimensions[0]; k++){
        for (int j=0; j<Dimensions[1]; j++){
            for (int i=0; i<Dimensions[2]; i++){
                int idx_Bz = IDX321_BZ(i, j, k, Dimensions[2], Dimensions[1]);
                Converted_MagZ[inx_Bz] = (double) 0.5 * ( Data[ idx_Bz ] + Data[ idx_Bz + Dimensions[2] * Dimensions[1] ] );
            }
        }
    }
}

#endif
#endif