## HDF5 VL Benchmarking Program

*Purpose*
 
 Determines writing and reading times for VL datatypes in HDF5. The program returns the time taken to:
        (1) Open, write and close a VL datatype dataset,
        (2) Open, read and close a VL datatype dataset.

*Function*
 
 The "-n" option specifies the size of the VL array. The "-v" option specifies the 
 maximum length of the variable length. The length of each array element is 
 randomly determined using the seed option, "-s". The "-p" option can set different 
 free space managers. The Paged FSM has additional parameters for the file space page 
 size (-p) and the buffer size of the page (-b). See the options below.

 Compile the program using the h5cc wrapper.

*OPTIONS*

```

   -r           read file [default - no]
   -w           write file [default - no]
   -n <int>     number of array elements [default = 1048576]
   -v <int>     maximum variable length [default = 1024]
   -s <seed>    seed for random number generator [default = 5]
   -f <int>     free space manager:
                      0 - FSM, Aggregators [default]
                      1 - Paged FSM
                      2 - Aggregators (no FSM)
                      3 - (no FSM)
   -p <int>      paged buffering (-f 1) option:
                      int - file space page size (in KiB)  [default = 4]
   -b <int>      paged buffering (-f 1) option:
                      int - buffer size of the page (in MiB) [default = 4]
   -h            help

``` 
