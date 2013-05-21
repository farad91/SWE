ncgen test.cdl -o test.nc
ncgen d_test.cdl -o d_test.nc
cxxtestgen --error-printer SWE_TsunamiScenarioTest.hh -o runner.cpp
g++ -o runner runner.cpp -lnetcdf
