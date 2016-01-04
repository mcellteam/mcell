
#include "xmlrpc_client.h"


char const *const reset = "nfsim.reset";
char const *const query = "nfsim.query";
char const *const init = "nfsim.init";
char const *const fire = "nfsim.fire";
char const *const step = "nfsim.step";
char const *const serverUrl = "http://localhost:8080/RPC2";
xmlrpc_env env;


/*
    Initialize the error handling environment and the global client object (oneclient for the whole mcell process)
    check xml-rpc API if we wish to use more clients
*/
int init_client(){
    xmlrpc_env_init(&env);
    /* Create the global XML-RPC client object. */
    xmlrpc_client_init2(&env, XMLRPC_CLIENT_NO_FLAGS, NAME, VERSION, NULL, 0);
    dieIfFaultOccurred(&env);
    return 0;
}


/*
clean up the error handling environment and shutdown the client library
*/
int close_client(){
    xmlrpc_env_clean(&env);
    
    xmlrpc_client_cleanup();

    return 0;
}

int init_nfsim(char* initString){
    xmlrpc_value * resultP;

    resultP = xmlrpc_client_call(&env, serverUrl, init,
                                 "(i)", (xmlrpc_int32) 5);
    dieIfFaultOccurred(&env);
    xmlrpc_DECREF(resultP);

    return 0;
}


int reset_nfsim(){
    xmlrpc_value * resultP;

    resultP = xmlrpc_client_call(&env, serverUrl, reset,
                                 "()");
    dieIfFaultOccurred(&env);
    xmlrpc_DECREF(resultP);
    return 0;
}
