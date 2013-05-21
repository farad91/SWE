ncgen test.cdl -o test.nc
cxxtestgen --error-printer SWE_TsunamiScenarioTest.hh -o runner.cpp
g++ -o runner -I$CXXTEST SWE_Scenario.hh SWE_TsunamiScenario.hh runner.cpp -lnetcdf
