#ifndef __wxtest__h__
#define __wxtest__h__

#ifdef _DO_USE_MPI_
#  include <mpi.h>
#endif

/**
 * Code in this file provides a simple test framework to rapidly test
 * various bits of code. Most functionality is provided via macros.
 */

// std includes
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

/**
 * Keeps track of how many tests have passed or failed
 */
class WxTestCounter
{
  public:
    WxTestCounter()
            : passed(0), failed(0) {
    }

    int passed;
    int failed;
    std::ofstream msgFile;

    std::vector<std::string> failedTests;

    void addFailedTest(int line, char *file, std::string s) {
        std::ostringstream lf;
        lf << s << " In file " << file ;
        failedTests.push_back(std::string(lf.str()));
    }

    void showFailedTests() {
        std::vector<std::string>::iterator i;
        for (i=failedTests.begin(); i!=failedTests.end(); ++i) {
            std::cout << "FAILED TEST: " << *i << std::endl;
            msgFile   << "FAILED TEST: " << *i << std::endl;
        }
    }
    
    void clearFailedTests() {
        failedTests.erase( failedTests.begin(), failedTests.end() );
    }
};
// `tc` is global to ease counting process
static WxTestCounter __tc;

/**
 * The WX_BEGIN_TESTS and WX_END_TEST macros can be used to sandwich a
 * set of WX_ASSERT and WX_RAISES macros defined below. These macros
 * count the number of passed and failed tests, printing them out
 * after the tests are done
*/
#ifdef _DO_USE_MPI_

#define WX_MPI_BEGIN_TESTS(file) \
 do {\
     __tc.passed = 0;\
     __tc.failed = 0;\
     int myrank; \
     std::ostringstream fname; \
     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);\
     fname << file << "-" << myrank << ".log";         \
     __tc.msgFile.open(fname.str().c_str(), std::fstream::out);\
 } while (0)

#define WX_MPI_END_TESTS \
 do {\
    MPI_Barrier(MPI_COMM_WORLD); \
    unsigned totalPassed, totalFailed; \
    MPI_Reduce(&__tc.passed, &totalPassed, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); \
    MPI_Reduce(&__tc.failed, &totalFailed, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); \
    int __myrank; \
    MPI_Comm_rank(MPI_COMM_WORLD, &__myrank); \
    if (__myrank == 0) { \
        std::cout << "PASSED = " << totalPassed <<  ". FAILED = " << totalFailed << std::endl; \
    } \
    __tc.showFailedTests(); \
    __tc.passed = 0;\
    __tc.failed = 0;\
    __tc.clearFailedTests(); \
 } while (0)
#else
# define WX_MPI_BEGIN_TESTS(file)
# define WX_MPI_END_TESTS
#endif // _DO_USE_MPI_

#define WX_BEGIN_TESTS(file)                    \
 do {\
     __tc.passed = 0;\
     __tc.failed = 0;\
     std::ostringstream fname;\
     fname << file << "-" << "0" << ".log";\
     __tc.msgFile.open(fname.str().c_str(), std::fstream::out); \
 } while (0)

#define WX_END_TESTS \
 do {\
     std::cout << "PASSED = " << __tc.passed <<  ". FAILED = " << __tc.failed << std::endl; \
     __tc.showFailedTests(); \
     __tc.passed = 0;\
     __tc.failed = 0;\
     __tc.clearFailedTests(); \
 } while (0)

/** 
 * The following macro can be used for running test cases. To test
 * some statement use, for example,
 *
 * WX_ASSERT("Running test 3", 2==2);
 *
 * The first macro parameter is the message while the second parameter
 * is the actual statement to test. This statement should be a
 * predicate, i.e. should evaluate to either true or false. If a test
 * fails the program does not halt.
 *
 */
#define WX_ASSERT(msg,expr) \
  do {\
   __tc.msgFile << msg << "\n";  \
   __tc.msgFile << "   Testing if " << #expr << "\n"; \
   if ((expr) == true)                               \
   {\
       __tc.msgFile << "     PASSED\n"; \
       __tc.passed++; \
   }\
   else \
   {\
       __tc.msgFile << "     FAILED\n"; \
       __tc.failed++; \
       __tc.addFailedTest(__LINE__, __FILE__, std::string(msg)+" [ "+std::string(#expr)+" ]"); \
   }\
  } while (0)

/**
 * The following macro can be used for running test cases. To test
 * some statement use, for example,
 *
 *  WX_RAISES("Running test 3", 2==2, const char *);
 *
 * The first macro parameter is the message while the second parameter
 * is the actual statement to test.  This statement should be a
 * predicate, i.e. should evaluate to either true or false. If a test
 * fails the program does not halt. The third parameter is a type
 * indicating the type of exception which needs to be caught.
 */
#define WX_RAISES(msg,expr,expc) \
do {\
    try\
    {\
        __tc.msgFile << msg << std::endl; \
        __tc.msgFile << "   Executing " << #expr << "...\n";                \
        __tc.msgFile << "   Expecting exception of type " << #expc << "...\n"; \
        expr; \
        __tc.msgFile << "     No exception thrown...\n";    \
        __tc.msgFile << "     FAILED\n";                    \
        __tc.failed++;\
       __tc.addFailedTest(__LINE__, __FILE__, std::string(msg)+" [ "+std::string(#expr)+" ]"); \
    }\
    catch(expc e)\
    {\
        __tc.msgFile << "     Caught exception of type " << #expc << "...\n";   \
        __tc.msgFile << "     PASSED\n"; \
        __tc.passed++;\
    }\
    catch(...)\
    {\
        __tc.msgFile << "     Caught exception of unexpected type\n"; \
        __tc.msgFile << "     FAILED\n";                              \
        __tc.failed++;\
       __tc.addFailedTest(__LINE__, __FILE__, std::string(msg)+" [ "+std::string(#expr)+" ]"); \
    }\
}while (0)

/**
 * Compares two arrays and returns true if they are equal.
 *
 * @param a first array to compare
 * @param b second array to compare
 * @param length size of arrays
 * @param true if arrays are the same
 */
template<typename T>
bool arraycmp(T a[], T b[], unsigned length)
{
    for (unsigned i = 0; i < length; ++i) 
    {
        if (a[i] != b[i]) 
            return false;
    }
    return true;
}

/**
 * Compares two arrays and returns true if they are equal.
 *
 * @param a first array to compare
 * @param b second array to compare
 * @return true if vectors are the same
 */
template<typename T>
bool arraycmp(const std::vector<T>& a, const std::vector<T>& b)
{
    if (a.size() != b.size())
        return false;
    for (unsigned i=0; i<a.size(); ++i)
    {
        if (a[i] != b[i]) 
            return false;
    }
    return true;
}

#endif // __wxtest__h__
