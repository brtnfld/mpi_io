#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdbool.h>

int
main (int argc, char *argv[] )
{
  FILE * fp;
  int *wdata;
  size_t dims[1];
  int GB=1;

  dims[0] = GB*1024*1024*1024*sizeof(int);

  wdata = malloc (dims[0]);

  /* open the file for writing*/
  fp = fopen ("scr","wb");
  fwrite(wdata, sizeof(*wdata),dims[0],fp);
  fclose(fp);
  
  return 0;
}
