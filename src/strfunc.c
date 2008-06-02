/*  Character string handling functions */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "strfunc.h"
#include "mem_util.h"


/*************************************************************************
my_strcat:
  In: two strings
  Out: the two strings concatenated, in a newly malloced block of memory,
       or NULL if there isn't enough memory.
  Note: the calling function is responsible for freeing the memory.
*************************************************************************/
char *my_strcat(char const *s1, char const *s2)
{ 
  char *temp = NULL;
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

/*************************************************************************
my_strclump:
  In: a NULL-terminated array of strings
  Out: all of the strings concatenated, in a newly malloced block of
       memory, or NULL if there isn't enough memory.
  Note: the calling function is responsible for freeing the memory.
*************************************************************************/
char *my_strclump(char **slist)
{
  int i,j,n,len;
  char **sp = NULL;
  char *s = NULL;
  char *temp = NULL;
  
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

/*************************************************************************
strip_quotes:
  In: a string that must be at least two characters long
  Out: a copy of the string, newly malloced, missing the first and
       last characters (might even be quotes!); NULL is returned if
       malloc fails.
  Note: this function does NOT do any error checking!
*************************************************************************/
char *strip_quotes(char const *s)
{ 
  char *temp = NULL;
  int len = strlen(s);
  
  if ((temp=(char *)malloc(len-1))!=NULL) {
    strncpy(temp,s+1,len-2);
    temp[len-2]='\0';
  }

  return(temp);
} 

/*
 * Format a string into an allocated buffer.
 */
char *alloc_vsprintf(char const *fmt, va_list args)
{
  char stack_buffer[256];
  int len;
  char *retval = NULL;
  va_list saved_args;
  va_copy(saved_args, args);

  len = vsnprintf(stack_buffer, sizeof(stack_buffer), fmt, args);
  if (len >= (int) sizeof(stack_buffer))
  {
    retval = malloc(len + 1);
    if (retval != NULL)
      vsnprintf(retval, len + 1, fmt, saved_args);
  }
  else
    retval = strdup(stack_buffer);

  return retval;
}

/*
 * Format a string into an allocated buffer.
 */
char *alloc_sprintf(char const *fmt, ...)
{
  char *retval;
  va_list args;
  va_start(args, fmt);
  retval = alloc_vsprintf(fmt, args);
  va_end(args);
  return retval;
}

