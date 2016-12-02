#include <mpi.h>
#include <iostream>
#include </working004a/an_sci2/Hao/eigen3/Eigen/Eigen>
#include <string>
#include <string>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>

//compile mpicxx -O3 -openmp crossproduct_mpi_mp.cpp -o crossproduct_mpi_mp.o
// run: ./cross.o 1000000 10000

using namespace std;


  void state(string statusfile,string message,clock_t ex, int rank);


int main ( int argc, char *argv[] )

{

    clock_t t1,t2;
    t1 = time(0);
    Eigen::initParallel();

    int n = atoi(argv[1]);  // number of observations
    int p = atoi(argv[2]);  // number of markers/covariates

    string OBS = argv[1];
    string MARKER = argv[2];

    string MPIS = argv[3];
    string MPS = argv[4];

    string statusfile("/working004a/an_sci2/Hao/OpenMP_MKL/output/" + OBS + "_" + MARKER + "_" + MPIS + "_" + MPS + ".txt");

    MPI::Init();
    int myRank=MPI::COMM_WORLD.Get_rank();
    int clusterSize=MPI::COMM_WORLD.Get_size();


//////// Hao's Stuff

    unsigned totalSize = n;
    unsigned batchSize = totalSize/clusterSize;
    unsigned iStart = myRank*batchSize;
    unsigned iEnd   = ( (myRank+1)==clusterSize ) ? totalSize :  iStart + batchSize;
    unsigned iSize  = iEnd - iStart;

    int recvcount[clusterSize];
    int displs[clusterSize];
    MPI::COMM_WORLD.Allgather(&iSize,1,MPI_INT,
                            &recvcount,1,MPI_INT);
    MPI::COMM_WORLD.Allgather(&iStart,1,MPI_INT,
                              &displs,1,MPI_INT);

///////////////

    t2 = time(0);
    state(statusfile,"Initializing: ",difftime(t2,t1),myRank);

// No need to globally assign ZtZFull - we can use MPI_IN_PLACE
    Eigen::MatrixXf ZtZ(p,p);

// Split Zfull in Z = Zfull(iSize, numberofmarkers)
    Eigen::MatrixXf Z = Eigen::MatrixXf::Random(iSize,p);

    t1 = time(0);
    state(statusfile,"Data Generation: ",difftime(t1,t2),myRank);

// This is OpenMP'd by Eigen
    ZtZ = Z.transpose()*Z;

    MPI::COMM_WORLD.Barrier();

    t2 = time(0);
    state(statusfile,"Crossproduct: ",difftime(t2,t1),myRank);


// Collapse subset ZtZ's to ZtZ of root-process (MPI_IN_PLACE)
    //if(myRank==0) {MPI::COMM_WORLD.Reduce(MPI_IN_PLACE, ZtZ.data(),ZtZ.size(), MPI_FLOAT,MPI_SUM,0); 
    //} else { MPI::COMM_WORLD.Reduce(ZtZ.data(), NULL ,ZtZ.size(), MPI_FLOAT,MPI_SUM,0); }
    
    if(myRank==0) {MPI::COMM_WORLD.Reduce(MPI_IN_PLACE, ZtZ.data(),1000, MPI_FLOAT,MPI_SUM,0); 
    } else { MPI::COMM_WORLD.Reduce(ZtZ.data(), NULL ,ZtZ.size(), MPI_FLOAT,MPI_SUM,0); }

    t1 = time(0);
    state(statusfile,"Collapsing: ",difftime(t1,t2),myRank);

    MPI::Finalize();



return(0);

}


// functions



void state(string statusfile,string message, clock_t ex, int rank)

{

  if(rank==0){
  std::ofstream status;
  status.open (statusfile.c_str(), std::ios_base::app);
  cout << message << ex <<  "\n\n";
  status << message <<  ex << "\n\n";
  status.close(); }

}














