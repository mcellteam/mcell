/*  Character string handling functions */

#include <string.h>
#include <stdlib.h>
#include "strfunc.h"

char *my_strdup(char *s)
{ 
  char *temp;
  
  if ((temp=(char *)malloc(strlen(s)+1))!=NULL) {
    strcpy(temp,s);
  }
  return(temp);
} 
  

char *my_strcat(char *s1, char *s2)
{ 
  char *temp;
  size_t len1,len2;

  len1=strlen(s1);
  len2=strlen(s2);
  if ((temp=(char *)malloc(len1+len2+1))!=NULL) {
    strcpy(temp,s1);
    strcpy(temp+len1,s2);
  }
  return(temp);
}

char *strip_quotes(char *s)
{ 
  char *temp;
  
  if ((temp=(char *)malloc(strlen(s)-1))!=NULL) {
    strncpy(temp,s+1,strlen(s)-2);
    strncpy(temp+strlen(s)-2,"",1);
  }
  return(temp);
} 

#ifndef DEBUG
void no_printf(const char *format, ...)
{
}
#endif

