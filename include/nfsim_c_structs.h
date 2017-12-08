#ifndef NFSIM_C_STRUCTS_H 
#define NFSIM_C_STRUCTS_H 




#ifdef __cplusplus
extern "C" {
#endif

    //C wrapper classes for C++ std objects we use for nfsim api calls
    void* map_create(void); 
    const char* map_get(void* map, const char* key);
    void* mapvector_create(void);
    int mapvector_size(void*);
    void* mapvector_get(void* vector, int position);

    int mapvectormap_size(void*);
    void mapvector_delete(void*);
    void* mapvectormap_create(void);
    void mapvectormap_delete(void* container);
    void* mapvectormap_get(void* container, char* reactant);
    char** mapvectormap_getKeys(void* container);


#ifdef __cplusplus
}
#endif

#endif
