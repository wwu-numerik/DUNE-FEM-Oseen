#include "config.h"

#include<dune/common/mpihelper.hh>
#include<iostream>
int main(int argc, char** argv)
{
  typedef Dune::MPIHelper Helper;
  
  {
    Helper& mpi = Helper::instance(argc, argv);
    
    Helper::MPICommunicator comm = mpi.getCommunicator();
    comm= mpi.getCommunicator();
  }
  
  {
    Helper& mpi = Helper::instance(argc, argv);

    Helper::MPICommunicator comm= mpi.getCommunicator();
    comm= mpi.getCommunicator();
  }
  std::cout << "We are at the end!"<<std::endl;
  
}
