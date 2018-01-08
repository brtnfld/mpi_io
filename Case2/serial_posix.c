#include<stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include<sys/types.h>
#include<time.h>

int main() {

  off_t expand_fs;
  FILE *fd;

  fd=fopen("datafile","w");
  fclose(fd);

  expand_fs = 33554432;

  clock_t tic = clock();
  truncate("datafile", expand_fs);
  clock_t toc = clock();

  printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

  return 0;
}
