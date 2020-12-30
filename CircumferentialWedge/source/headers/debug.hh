#ifndef DEBUGPRINT_H
#define DEBUGPRINT_H

#ifdef DEBUG

#include <iostream>
#define DEBUG_MSG(str) do { std::cout << str << std::endl; } while( false )

#else

#define DEBUG_MSG(str) do { } while ( false )

#endif // DEBUG
#endif // DEBUGPRINT_H


