
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <cstdlib>

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
    LLONG,
    FLOAT,
    DOUBLE,
    COMPLEXDOUBLE,
    COMPLEXFLOAT,
    COMPLEXINT
} TYPE;

typedef struct
{
    char  **names; //Names to be stored
    int   id;      //The id of where we are
} iterate;

int the_create(char *name);
char* fix_ext (char *name);
void readDataset(TYPE type, char *name, char* dataset, void* data);
void readDatasets(TYPE type, char *name, char** dataset, void* data);
void createDataset(TYPE type, char *name, char* dataset, int rank, int* dims, void* data);
hid_t complextype(TYPE type);
// Operator function to be called by H5Literate.
herr_t lister (hid_t loc_id, const char *name, const H5L_info_t *info, void *operator_data);
herr_t counter (hid_t loc_id, const char *name, const H5L_info_t *info,void *count);  

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

extern "C" void readDatasetLLong( char *name, char* dataset, long long* data)
{
    readDataset(LLONG, name, dataset, (void*)data);
}

extern "C" void readDatasetInt( char *name, char* dataset, int* data)
{
    readDataset(INT, name, dataset, (void*)data);
}

extern "C" void readDatasetDouble( char *name, char* dataset, double* data)
{
    readDataset(DOUBLE, name, dataset, (void*)data);
}

extern "C" void readDatasetComplex( char *name, char* dataset, complexDouble* data)
{
    readDataset(COMPLEXDOUBLE, name, dataset, (void*)data);
}

extern "C" void readDatasetsDouble( char *name, char** datasets, double* data)
{
    readDatasets(DOUBLE, name, datasets, (void*)data);
}

extern "C" void readDatasetsComplex( char *name, char** datasets, complexDouble* data)
{
    readDatasets(COMPLEXDOUBLE, name, datasets, (void*)data);
}


extern "C" void createDatasetInt( char *name, char* dataset, int rank, int* dims, int* data)
{
    createDataset(INT, name, dataset, rank, dims, (void*)data);
}

extern "C" void createDatasetLLong( char *name, char* dataset, int rank, int* dims, long long* data)
{
    createDataset(LLONG, name, dataset, rank, dims, (void*)data);
}

extern "C" void createDatasetDouble( char *name, char* dataset, int rank, int* dims, double* data)
{
    createDataset(DOUBLE, name, dataset, rank, dims, (void*)data);
}

extern "C" void createDatasetComplex(char *name, char* dataset, int rank, int* dims, complexDouble* data)
{
    createDataset(COMPLEXDOUBLE, name, dataset, rank, dims, data);
}

extern "C" int getRankDataset( char *name, char* dataset)
{
    hid_t      file_id;   /* file identifier */
    name = fix_ext(name);
    int rank;
    
     /* open file */
    file_id = H5Fopen (name, H5F_ACC_RDONLY, H5P_DEFAULT);
    /* get the dimensions of the dataset */
    H5LTget_dataset_ndims(file_id,dataset,&rank);
    /* close file */
    H5Fclose (file_id);
    return rank;
}

extern "C" void getDimsDataset( char *name, char* dataset, int rank, int* dims)
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

extern "C" char** listGroupMembers (char *name, char* groupname)
{
    hid_t           file, group;           /* Handle */
    herr_t          status;
    int count;
    count = 0;
    iterate dat;
    name = fix_ext(name);
    
    // Open file.
    file = H5Fopen (name, H5F_ACC_RDONLY, H5P_DEFAULT);
    group = H5Gopen(file, groupname, H5P_DEFAULT);
    
    // Begin iteration.
    status = H5Literate (group, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, counter, (void *) &count);
    
    char** operator_data;
    operator_data = (char **) malloc((count + 1) * sizeof(char *));
    dat.id = 0;
    dat.names = operator_data;
    
    status = H5Literate (group, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, lister, (void *) &dat);
    
    operator_data[count] = NULL;
    
    //Close and release resources.
    status = H5Gclose (group);
    status = H5Fclose (file);
    
    return operator_data;
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

void createDataset( TYPE type, char *name, char* dataset, int rank, int* dims, void* data)
{
    hid_t       file_id, memtype;   /* file identifier */
    hsize_t    thedims[rank];
    switch(type){
        case CHAR: memtype = H5T_NATIVE_CHAR; break;
        case INT: memtype = H5T_NATIVE_INT ; break;
        case LLONG : memtype = H5T_NATIVE_LLONG; break;
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

void readDataset( TYPE type, char *name, char* dataset, void* data)
{
    hid_t      file_id, memtype;   /* file identifier */
    name = fix_ext(name);
    switch(type){
        case CHAR: memtype = H5T_NATIVE_CHAR; break;
        case INT: memtype = H5T_NATIVE_INT ; break;
        case LLONG : memtype = H5T_NATIVE_LLONG; break;
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

void readDatasets( TYPE type, char *name, char** datasets, void* data)
{
    hid_t      file_id, memtype;   /* file identifier */
    name = fix_ext(name);
    int rank;
    int datasize;

    switch(type){
        case CHAR: memtype = H5T_NATIVE_CHAR;
                   datasize = sizeof(char); break;
        case INT: memtype = H5T_NATIVE_INT;
                  datasize = sizeof(int); break;
        case FLOAT: memtype = H5T_NATIVE_FLOAT ;
                    datasize = sizeof(float); break;
        case DOUBLE: memtype = H5T_NATIVE_DOUBLE ;
                     datasize = sizeof(double); break;
        case COMPLEXDOUBLE: memtype = complextype(type) ;
                            datasize = 2*sizeof(double); break;
        case COMPLEXFLOAT: memtype = complextype(type) ;
                           datasize = 2*sizeof(float); break;
        case COMPLEXINT: memtype = complextype(type) ;
                         datasize = 2*sizeof(int); break;
    }
    
     /* open file */
    file_id = H5Fopen (name, H5F_ACC_RDONLY, H5P_DEFAULT);
    /* get information about ranks and dims of the first dataset*/
    char* dataset0 = datasets[0];
    H5LTget_dataset_ndims(file_id,dataset0,&rank);
    hsize_t dims[rank];
    H5LTget_dataset_info(file_id,dataset0,dims,NULL,NULL);
    int sizedataset = 1;
    for (int i = 0; i < rank; i++)
        sizedataset *= (int) dims[i];
    //We cannot do pointer arithemtics on void, so cast it to char, which is always 1 byte big
    char *currentdata = (char *) data;
    int i = 0;
    /* read all the datasets */
    while(datasets[i] != NULL)
    {
        H5LTread_dataset(file_id,datasets[i++],memtype, (void *)currentdata);
        currentdata += datasize * sizedataset;
    }
    /* close file and memory type */
    switch(type){
        case COMPLEXDOUBLE: 
        case COMPLEXFLOAT:
        case COMPLEXINT: H5Tclose(memtype); break;
    }
    H5Fclose (file_id);
}

// Operator function for iterator
herr_t lister (hid_t loc_id, const char *name, const H5L_info_t *info,void *operator_data)
{
    iterate *dat = (iterate *) operator_data;
    char ** strlist = dat->names;
    
    strlist[dat->id] = (char *) malloc(strlen(name) + 1);
    strcpy(strlist[dat->id], name);
    dat->id++;
    return 0;
}

herr_t counter (hid_t loc_id, const char *name, const H5L_info_t *info,void *count)
{
    ((int *) count)[0]++;
    return 0;
}

char* fix_ext(char* name)
{
    char *dot = strrchr(name, '.');
    if (!(dot && !strcmp(dot, ".h5")))
        name = strcat(name, ".h5");
    return name;
}


/*
// For testing purposes
int main()
{
    char* in = (char*) malloc(strlen("test")+1); // +1 for the terminator
    strcpy(in,"test");
    char* group = (char*) malloc(strlen("/")+1); // +1 for the terminator
    strcpy(group,"/");
    
    char** data; // data[20];
    
    data = listGroupMembers(in, group);

    int rank = getRankDataset(in, *data);
    int dims[rank];
    getDimsDataset(in, *data, rank, dims);
    int sizedataset = 1;
    for (int i = 0; i < rank; i++)
        sizedataset *= (int) dims[i];
    double alldata[sizedataset * 2];
    readDatasets(DOUBLE, in, data, alldata);
    
    while(*data != NULL)
    {
        printf("%s\n", *data);
        data += 1;
    }
    for(int i =0; i < sizedataset * 2; i++)
    {
        printf("%f\n", alldata[i]);
    }
    return 0;
}*/