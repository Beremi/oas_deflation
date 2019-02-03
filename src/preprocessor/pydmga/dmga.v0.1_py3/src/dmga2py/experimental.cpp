/*
 * dmga.cpp
 *
 * Python bindings C-file
 *
 *  Created on: 22-04-2013
 *      Author: robson
 */
#include <Python.h>
#include <stdint.h>

#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

int
main(int argc, char *argv[])
{
	std::cout << "Py_ssize_t = " << sizeof(Py_ssize_t) << "\n";
	std::cout << "intptr_t = " << sizeof(intptr_t) << "\n";
	std::cout << "int = " << sizeof(int) << "\n";
	std::cout << "long long = " << sizeof(long long) << "\n";
	std::cout << "long int = " << sizeof(long int) << "\n";

	return 0;
}



