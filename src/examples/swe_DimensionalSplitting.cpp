/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *         Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Univ.-Prof._Dr._Michael_Bader)
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * Basic setting of SWE, which uses a wave propagation solver and an artificial or ASAGI scenario on a single block.
 */

#include <cassert>
#include <cstdlib>
#include <string>
#include <iostream>

//just our block
#include "blocks/SWE_DimensionalSplitting.hh"

// input data
#ifdef TSUNAMINC
#include "scenarios/SWE_TsunamiScenario.hh"
#else
#include "scenarios/SWE_simple_scenarios.hh"
#endif

// output data just NetCDF
#include "writer/BoyeWriter.hh"
#include "writer/NetCdfWriter.hh"



#include "tools/help.hh"
#include "tools/Logger.hh"
#include "tools/ProgressBar.hh"

/**
 * Main program for the simulation on a single SWE_WavePropagationBlock.
 */
int main( int argc, char** argv ) {
  /**
   * Initialization.
   */

// check if the necessary command line input parameters are given 
#ifdef TSUNAMINC
  if(argc != 7 && argc != 8) {
    std::cout << "Aborting ... please provide proper input parameters." << std::endl
              << "Example: ./SWE_parallel 200 300 /work/openmp_out bat.nc dis.nc 20" << std::endl
              << "\tfor a single block of size 200 * 300 and simulation time 20sec" << std::endl
              << "Custom number of checkpoints (in this case 50):" << std::endl
              << "Example: ./SWE_parallel 200 300 /work/openmp_out bat.nc dis.nc 20 50" << std::endl;
    return 1;
  }
#else
  if(argc != 4) {
    std::cout << "Aborting ... please provide proper input parameters." << std::endl
              << "Example: ./SWE_parallel 200 300 /work/openmp_out" << std::endl
              << "\tfor a single block of size 200 * 300" << std::endl;
    return 1;
  }
#endif

  //! number of grid cells in x- and y-direction.
  int l_nX, l_nY;
  
  //! output file of a previous run is existing?
  bool checkpoint;

  //! time when the simulation will be stopped (in seconds)
  float l_endSimulation;
  
  //! number of checkpoints for visualization (at each checkpoint in time, an output file is written).
  int l_numberOfCheckPoints;
  
  //! l_baseName of the plots.
  std::string l_baseName;

  // read command line parameters
  l_nY = l_nX = atoi(argv[1]);
  l_nY = atoi(argv[2]);
  l_baseName = std::string(argv[3]);
  
#ifdef TSUNAMINC
  l_endSimulation  = atoi(argv[6]);
#else
  l_endSimulation = l_scenario->endSimulation();
#endif
  
  if(argc == 8)
    l_numberOfCheckPoints = atoi(argv[7]);
  else
    l_numberOfCheckPoints = 20;

  
#ifdef CHECKPOINT
  // check if an output/checkpoint file of a previous run is existing
  ifstream cp_file;
  string cp_file_name = l_fileName + ".nc";
  // try to open the checkpoint file
  cp_file.open(cp_file_name.c_str(), ifstream::out);
  // check for errors when opening the file
  checkpoint = cp_file.good();
  cp_file.close();
#else
  // assume that no output file is existing
  checkpoint = false;
#endif
  
#ifdef DEBUG
  cerr << "Creating scenario..." << endl;
#endif
 
#ifdef TSUNAMINC
  // load bathymetry (argv[4]) and displacement (argv[5]) from netCDF files
  SWE_NetCDFScenario *l_scenario = new SWE_TsunamiScenario();
  if(l_scenario->readNetCDF(argv[4],argv[5]) != 0) {
    cerr << "Please specify correct input files!" << endl
         << "Either " << argv[4] << " or " << argv[5] << " are unusable!" << endl;
    std::exit(1);
  }
#else
  // create a simple artificial scenario
  SWE_DamBreakScenario *l_scenario = new SWE_DamBreakScenario();
#endif

  
  //! size of a single cell in x- and y-direction
  float l_dX, l_dY;

  // compute the size of a single cell
  l_dX = (l_scenario->getBoundaryPos(BND_RIGHT) - l_scenario->getBoundaryPos(BND_LEFT)   )/l_nX;
  l_dY = (l_scenario->getBoundaryPos(BND_TOP) -   l_scenario->getBoundaryPos(BND_BOTTOM) )/l_nY;

  // create a single wave propagation block
  SWE_DimensionalSplitting l_wavePropgationBlock(l_nX, l_nY, l_dX, l_dY);

  //! origin of the simulation domain in x- and y-direction
  float l_originX, l_originY;

  // get the origin from the scenario
  l_originX = l_scenario->getBoundaryPos(BND_LEFT);
  l_originY = l_scenario->getBoundaryPos(BND_BOTTOM);

  
  // initialize the wave propagation block
  l_wavePropgationBlock.initScenario(l_originX, l_originY, *l_scenario);


  //! checkpoints when output files are written.
  float* l_checkPoints = new float[l_numberOfCheckPoints+1];

  // compute the checkpoints in time
  for(int cp = 0; cp <= l_numberOfCheckPoints; cp++) {
     l_checkPoints[cp] = cp*(l_endSimulation/l_numberOfCheckPoints);
  }

  // Init fancy progressbar
  tools::ProgressBar progressBar(l_endSimulation);

  // write the output at time zero
  tools::Logger::logger.printOutputTime((float) 0.);
  progressBar.update(0.);

  std::string l_fileName = generateBaseFileName(l_baseName,0,0);
  //boundary size of the ghost layers
  io::BoundarySize l_boundarySize = {{1, 1, 1, 1}};

  
  io::BoyeWriter l_boyeWriter( l_fileName,2);
  l_boyeWriter.initBoye(0,0,l_wavePropgationBlock,0);
  //l_boyeWriter.initBoye(10,0,l_wavePropgationBlock,1);
  //l_boyeWriter.initBoye(0,1,l_wavePropgationBlock,2);
  l_boyeWriter.writeBoye(0,l_wavePropgationBlock);
  
  #ifdef DYNAMIC
  l_wavePropgationBlock.updateBathymetry(*l_scenario, 0.f);
  // construct a NetCdfWriter
  io::NetCdfWriter l_writer( l_fileName,
                             l_wavePropgationBlock.getBathymetry(),
                             l_boundarySize,
                             l_nX, l_nY,
                             l_dX, l_dY,
                             l_endSimulation,
                             !checkpoint,true,
                             l_originX, l_originY,
                             0
                           );
  // Write zero time step
  l_writer.writeTimeStep( l_wavePropgationBlock.getWaterHeight(),
                          l_wavePropgationBlock.getDischarge_hu(),
                          l_wavePropgationBlock.getDischarge_hv(),
                          l_wavePropgationBlock.getBathymetry(),
                          (float) 0.);
  #else
    // construct a NetCdfWriter
  io::NetCdfWriter l_writer( l_fileName,
                             l_wavePropgationBlock.getBathymetry(),
                             l_boundarySize,
                             l_nX, l_nY,
                             l_dX, l_dY,
                             l_endSimulation,
                             !checkpoint,false,
                             l_originX, l_originY,
                             0
                           );
  l_writer.writeTimeStep( l_wavePropgationBlock.getWaterHeight(),
                          l_wavePropgationBlock.getDischarge_hu(),
                          l_wavePropgationBlock.getDischarge_hv(),
                          (float) 0.);
  #endif
  

  /**
   * Simulation.
   */
  // print the start message and reset the wall clock time
  progressBar.clear();
  tools::Logger::logger.printStartMessage();
  tools::Logger::logger.initWallClockTime(time(NULL));

  //! simulation time.
  float l_t = 0.f;
  int   c_h = 1;
#ifdef DYNAMIC
  int CP_Start = 0;
#endif
  progressBar.update(l_t);

  unsigned int l_iterations = 0;

  // loop over checkpoints
  for(int c=c_h; c<=l_numberOfCheckPoints; c++) {

    // do time steps until next checkpoint is reached
    while( l_t < l_checkPoints[c] ) {
      // set values in ghost cells:
      l_wavePropgationBlock.setGhostLayer();
      
      // reset the cpu clock
      tools::Logger::logger.resetCpuClockToCurrentTime();
      #ifdef DYNAMIC
      if(l_t<=l_scenario->getEruptionDuration()){
        l_wavePropgationBlock.runTimestep(l_scenario->getEruptionResolution());
      }
      else{
            l_wavePropgationBlock.runTimestep();
      }
      #else
      l_wavePropgationBlock.runTimestep();
      #endif
      float l_maxTimeStepWidth = l_wavePropgationBlock.getMaxTimestep();
      // update the cpu time in the logger
      tools::Logger::logger.updateCpuTime();

      // update simulation time with time step width.
      l_t += l_maxTimeStepWidth;
      l_iterations++;
      #ifdef DYNAMIC
      if(l_t<=l_scenario->getEruptionDuration()){
        l_wavePropgationBlock.updateBathymetry(*l_scenario, l_t);
      }
      #endif
      // print the current simulation time
      progressBar.clear();
      tools::Logger::logger.printSimulationTime(l_t);
      progressBar.update(l_t);
      l_boyeWriter.writeBoye(l_t, l_wavePropgationBlock);
      #ifdef DYNAMIC
      if(l_t<=l_scenario->getEruptionDuration() && CP_Start == 30){
        progressBar.clear();
        tools::Logger::logger.printOutputTime(l_t);
        progressBar.update(l_t);
        l_writer.writeTimeStep( l_wavePropgationBlock.getWaterHeight(),
                        l_wavePropgationBlock.getDischarge_hu(),
                        l_wavePropgationBlock.getDischarge_hv(),
                        l_wavePropgationBlock.getBathymetry(),
                        l_t);
        CP_Start = 0;
      }
        CP_Start++;
      #endif
    }
  #ifdef DYNAMIC
    if(l_t>l_scenario->getEruptionDuration()){
    // print current simulation time of the output
    progressBar.clear();
    tools::Logger::logger.printOutputTime(l_t);
    progressBar.update(l_t);
  
    // write output
    l_writer.writeTimeStep( l_wavePropgationBlock.getWaterHeight(),
                            l_wavePropgationBlock.getDischarge_hu(),
                            l_wavePropgationBlock.getDischarge_hv(),
                            l_wavePropgationBlock.getBathymetry(),
                            l_t);
    }                        
  #else
    progressBar.clear();
    tools::Logger::logger.printOutputTime(l_t);
    progressBar.update(l_t);
  
    l_writer.writeTimeStep( l_wavePropgationBlock.getWaterHeight(),
                            l_wavePropgationBlock.getDischarge_hu(),
                            l_wavePropgationBlock.getDischarge_hv(),
                            l_t);
  #endif
  }

  /**
   * Finalize.
   */
  // write the statistics message
  progressBar.clear();
  tools::Logger::logger.printStatisticsMessage();

  // print the cpu time
  tools::Logger::logger.printCpuTime();

  // print the wall clock time (includes plotting)
  tools::Logger::logger.printWallClockTime(time(NULL));

  // printer iteration counter
  tools::Logger::logger.printIterationsDone(l_iterations);
  
  // free dynamically allocated memory
  delete l_scenario;
  
  return 0;
}
