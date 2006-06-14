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

  len1 = (s1==NULL) ? 0 : strlen(s1);
  len2 = (s2==NULL) ? 0 : strlen(s2);
  if ((temp=(char *)malloc(len1+len2+1))!=NULL) {
    if (len1) strcpy(temp,s1);
    if (len2) strcpy(temp+len1,s2);
    temp[len1+len2] = '\0';
  }
  return(temp);
}

char *my_strclump(char **slist)
{
  int i,j,n,len;
  char **sp;
  char *s;
  char *temp;
  
  for (sp=slist,n=0 ; *sp!=NULL ; sp++,n++);
  
  for (i=0,len=0;i<n;i++) len += strlen(slist[i]);
  
  temp = (char*) malloc(len+1);
  if (temp==NULL) return NULL;
  
  j=0;
  for (sp=slist;*sp!=NULL;sp++)
  {
    for (s=*sp ; *s!=0 ; s++)
    {
      temp[j++] = *s;
      if (j==len) { temp[j]=0; return temp; }
    }
  }
  
  temp[j]=0;
  return temp;
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

