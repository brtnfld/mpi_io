#include<stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include<sys/types.h>

int ctrunc(off_t expand_fs) {

  truncate("datafile", expand_fs);

/*   FILE *fd; */
/*   fd=fopen("datafile","w"); */
/*   ftruncate(fileno(fd), expand_fs); */
/*   fclose(fd); */

  return 0;
}

