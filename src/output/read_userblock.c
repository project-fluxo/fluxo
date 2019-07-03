#include <string.h>
#include <stdio.h>
extern char userblock_start;
extern char userblock_end;
extern char userblock_size;

long get_userblock_size_(void) 
{
   return (unsigned long)(&userblock_size);
}

long get_inifile_size(char* filename) 
{
   FILE* fp = fopen(filename, "rb");
   fseek(fp,0,SEEK_END);
   long length=ftell(fp);
   fclose(fp);
   return length;
}

void insert_userblock(char* filename, char* inifilename) 
{
   FILE* fp = fopen(filename, "w");
   rewind(fp);
   fprintf(fp, "{[( START USERBLOCK )]}\n");

   // ini file
   fprintf(fp, "{[( INIFILE )]}\n");
   FILE* fini = fopen(inifilename, "rb");
   int c;
   do {
      c = fgetc (fini);
      if (c != EOF) fputc((char)c, fp);
   } while (c != EOF);
   fclose(fini);
   // compressed data ( as tar.xz )
   fprintf(fp, "{[( COMPRESSED )]}\n");
   fprintf(fp, "userblock.txt\n");      // filename
   fprintf(fp, "%ld\n", get_userblock_size_()); // filesize
   char* p = &userblock_start;
   while ( p != &userblock_end ) fputc(*p++, fp);
   fprintf(fp, "\n");

   fprintf(fp, "{[( END USERBLOCK )]}\n");
   fclose(fp);
}

void copy_userblock(char* outfilename, char* infilename)
{
   FILE* fout = fopen(outfilename,"rb+");
   FILE* fin  = fopen(infilename, "r");
   rewind(fout);
   int c;
   do {
      c = fgetc (fin);
      if (c != EOF) fputc((char) c, fout);
   } while (c != EOF);
   fclose(fin);
   fclose(fout);
}
