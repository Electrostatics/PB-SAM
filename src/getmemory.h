#ifndef _GETMEMORY_H_
#define _GETMEMORY_H_

#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>

//! process_mem_usage function
/*!  Function to calculate the process memory usage
 \param vm_usage a reference of a double of the virtual memory
 usage
 \param resident_set a reference of a double of the resident set size,
 the memory that is held in RAM */
void process_mem_usage(double& vm_usage, double& resident_set);

#endif
