/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/TestMain.h>
#include <cxxtest/ErrorPrinter.h>

int main( int argc, char *argv[] ) {
 int status;
    CxxTest::ErrorPrinter tmp;
    status = CxxTest::Main<CxxTest::ErrorPrinter>( tmp, argc, argv );
    return status;
}
bool suite_DimensionalSplittingTest_init = false;
#include "DimensionalSplittingTest.t.h"

static DimensionalSplittingTest suite_DimensionalSplittingTest;

static CxxTest::List Tests_DimensionalSplittingTest = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_DimensionalSplittingTest( "DimensionalSplittingTest.t.h", 11, "DimensionalSplittingTest", suite_DimensionalSplittingTest, Tests_DimensionalSplittingTest );

static class TestDescription_suite_DimensionalSplittingTest_testCompareNetUpdates : public CxxTest::RealTestDescription {
public:
 TestDescription_suite_DimensionalSplittingTest_testCompareNetUpdates() : CxxTest::RealTestDescription( Tests_DimensionalSplittingTest, suiteDescription_DimensionalSplittingTest, 13, "testCompareNetUpdates" ) {}
 void runTest() { suite_DimensionalSplittingTest.testCompareNetUpdates(); }
} testDescription_suite_DimensionalSplittingTest_testCompareNetUpdates;

static class TestDescription_suite_DimensionalSplittingTest_testHorizontalSteadyState : public CxxTest::RealTestDescription {
public:
 TestDescription_suite_DimensionalSplittingTest_testHorizontalSteadyState() : CxxTest::RealTestDescription( Tests_DimensionalSplittingTest, suiteDescription_DimensionalSplittingTest, 63, "testHorizontalSteadyState" ) {}
 void runTest() { suite_DimensionalSplittingTest.testHorizontalSteadyState(); }
} testDescription_suite_DimensionalSplittingTest_testHorizontalSteadyState;

#include <cxxtest/Root.cpp>
const char* CxxTest::RealWorldDescription::_worldName = "cxxtest";
