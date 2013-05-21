ncgen test.cdl -o test.nc
cxxtestgen --error-printer SWE_TsunamiScenarioTest.hh -o runner.cpp
g++ -o runner runner.cpp -lnetcdf
