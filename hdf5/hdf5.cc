
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "hdf5.h"
#include "hdf5_hl.h"

///////////////////////////////////////////////////////////////
// Declaration of functions and structs we need
///////////////////////////////////////////////////////////////
typedef struct {
    double r;
    double i;
} complexDouble;

typedef struct {
    float r;
    float i;
} complexFloat;

typedef struct {
    int r;
    int i;
} complexInt;

typedef enum {
    CHAR,
    INT,
    FLOAT,
    DOUBLE,
    COMPLEXDOUBLE,
    COMPLEXFLOAT,
    COMPLEXINT
} TYPE;

int the_create(char *name);
char* fix_ext (char *name);
void readDataset(TYPE type, char *name, char* dataset, void* data);
void createDataset(TYPE type, char *name, char* dataset, int rank, int* dims, void* data);
hid_t complextype(TYPE type);

///////////////////////////////////////////////////////////////
//The external bindings
///////////////////////////////////////////////////////////////
extern "C" void createh5File(char *name)
{
    hid_t       file_id;   /* file identifier */
    herr_t      status;
    name = fix_ext(name);
        

    /* Create a new file using default properties. */
    file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Terminate access to the file. */
    status = H5Fclose(file_id); 
}

extern "C" void readDatasetDouble( char *name, char* dataset, double* data)
{
    readDataset(DOUBLE, name, dataset, (void*)data);
}

extern "C" void readDatasetComplex( char *name, char* dataset, complexDouble* data)
{
    readDataset(COMPLEXDOUBLE, name, dataset, (void*)data);
}

extern "C" void createDatasetDouble( char *name, char* dataset, int rank, int* dims, double* data)
{
    createDataset(DOUBLE, name, dataset, rank, dims, (void*)data);
}

extern "C" void createDatasetComplex(char *name, char* dataset, int rank, int* dims, complexDouble* data)
{
    createDataset(COMPLEXDOUBLE, name, dataset, rank, dims, data);
}

extern "C" int getRankDataset( char *name
                  , char* dataset
                  )
{
    hid_t      file_id;   /* file identifier */
    name = fix_ext(name);
    int rank[1];
    
     /* open file */
    file_id = H5Fopen (name, H5F_ACC_RDONLY, H5P_DEFAULT);
    /* get the dimensions of the dataset */
    H5LTget_dataset_ndims(file_id,dataset,rank);
    /* close file */
    H5Fclose (file_id);
    return rank[0];
}

extern "C" void getDimsDataset( char *name
                  , char* dataset
                  , int rank
                  , int* dims
                  )
{
    hid_t      file_id;   /* file identifier */
    hsize_t    h_dims[rank];
    name = fix_ext(name);
    
     /* open file */
    file_id = H5Fopen (name, H5F_ACC_RDONLY, H5P_DEFAULT);
    /* get the dimensions of the dataset */
    H5LTget_dataset_info(file_id,dataset,h_dims,NULL,NULL);
    /* get the dimension in ints */
    for (int i = 0; i < rank; i++)
        dims[i] = h_dims[i];
    /* close file */
    H5Fclose (file_id);
}

//Helpers 

//Make sure to close it after use with H5Tclose
hid_t complextype(TYPE type){
    herr_t status;
    hid_t  memtype;
    
    switch(type){
        case COMPLEXDOUBLE : memtype = H5Tcreate (H5T_COMPOUND, sizeof (complexDouble));
                            status = H5Tinsert (memtype, "r", HOFFSET (complexDouble, r), H5T_NATIVE_DOUBLE);
                            status = H5Tinsert (memtype, "i", HOFFSET (complexDouble, i), H5T_NATIVE_DOUBLE);
                            break;
        case COMPLEXFLOAT : memtype = H5Tcreate (H5T_COMPOUND, sizeof (complexFloat));
                            status = H5Tinsert (memtype, "r", HOFFSET (complexFloat, r), H5T_NATIVE_FLOAT);
                            status = H5Tinsert (memtype, "i", HOFFSET (complexFloat, i), H5T_NATIVE_FLOAT);
                            break;
        case COMPLEXINT   : memtype = H5Tcreate (H5T_COMPOUND, sizeof (complexInt));
                            status = H5Tinsert (memtype, "r", HOFFSET (complexInt, r), H5T_NATIVE_INT);
                            status = H5Tinsert (memtype, "i", HOFFSET (complexInt, i), H5T_NATIVE_INT);
                            break;
    }
    return memtype;
}

void createDataset( TYPE type
                  , char *name
                  , char* dataset
                  , int rank
                  , int* dims
                  , void* data
                  )
{
    hid_t       file_id, memtype;   /* file identifier */
    hsize_t    thedims[rank];
    switch(type){
        case CHAR: memtype = H5T_NATIVE_CHAR; break;
        case INT: memtype = H5T_NATIVE_INT ; break;
        case FLOAT: memtype = H5T_NATIVE_FLOAT ; break;
        case DOUBLE: memtype = H5T_NATIVE_DOUBLE ; break;
        case COMPLEXDOUBLE: 
        case COMPLEXFLOAT:
        case COMPLEXINT: memtype = complextype(type) ; break; 
    }
    
    for (int i = 0; i < rank; i++)
        thedims[i] = (hsize_t) dims[i];
    name = fix_ext(name);
     /* open file */
    file_id = H5Fopen (name, H5F_ACC_RDWR, H5P_DEFAULT);
    /* create and write a double type dataset */
    H5LTmake_dataset(file_id,dataset,rank,thedims,memtype,data);
    /* close file and memory type */
    switch(type){
        case COMPLEXDOUBLE: 
        case COMPLEXFLOAT:
        case COMPLEXINT: H5Tclose(memtype); break;
    }
    H5Fclose (file_id);
}

void readDataset( TYPE type
                , char *name
                , char* dataset
                , void* data
                )
{
    hid_t      file_id, memtype;   /* file identifier */
    name = fix_ext(name);
    switch(type){
        case CHAR: memtype = H5T_NATIVE_CHAR; break;
        case INT: memtype = H5T_NATIVE_INT ; break;
        case FLOAT: memtype = H5T_NATIVE_FLOAT ; break;
        case DOUBLE: memtype = H5T_NATIVE_DOUBLE ; break;
        case COMPLEXDOUBLE: 
        case COMPLEXFLOAT:
        case COMPLEXINT: memtype = complextype(type) ; break; 
    }
    
     /* open file */
    file_id = H5Fopen (name, H5F_ACC_RDONLY, H5P_DEFAULT);
    /* read dataset */
    H5LTread_dataset(file_id,dataset,memtype,data);
    /* close file and memory type */
    switch(type){
        case COMPLEXDOUBLE: 
        case COMPLEXFLOAT:
        case COMPLEXINT: H5Tclose(memtype); break;
    }
    H5Fclose (file_id);
    
}


// Helpers
char* fix_ext(char* name)
{
    char *dot = strrchr(name, '.');
    if (!(dot && !strcmp(dot, ".h5")))
        name = strcat(name, ".h5");
    return name;
}