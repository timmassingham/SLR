#ifndef _MYSTRING_H_
#define _MYSTRING_H_

struct __mystring_struct {
   char * string;
   int len, maxlen;
};
typedef struct __mystring_struct * Mystring;
#define Mystring_size sizeof(struct __mystring_struct)

Mystring new_mystring (const int len);
void free_mystring (Mystring string);
void append_char_to_mystring ( const char c, Mystring string);
char * cstring_of_mystring(const Mystring string);
Mystring mystring_of_cstring (const char * str);
#endif

