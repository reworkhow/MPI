#include <iostream>
#include <Eigen>
#include <string>
#include <iostream>
#include <fstream>
#include <matvec/statdist.h>
#include <matvec/session.h>
#include <mpi.h>
#include <ctime>

typedef Eigen::Map<Eigen::VectorXf>     MVectorXf;

using namespace std;

int main(int argc, const char * argv[]){
    MPI::Init();

    int myRank=MPI::COMM_WORLD.Get_rank();
    int clusterSize=MPI::COMM_WORLD.Get_size();
    const unsigned numberMarkers = 3;
    const unsigned numObs        = 4;
    
    unsigned totalSize = numObs;
    unsigned batchSize = totalSize/clusterSize;
    unsigned iStart = myRank*batchSize;
    unsigned iEnd   = ( (myRank+1)==clusterSize ) ? totalSize :  iStart + batchSize;
    unsigned iSize  = iEnd - iStart;

    Eigen::MatrixXf Z(iSize,numberMarkers);
    Eigen::VectorXf means(numberMarkers);
    Eigen::VectorXf meansFull(numberMarkers);
 


    Z = Eigen::MatrixXf::Random(iSize,numberMarkers);
    means = Z.colwise().mean();
    MPI::COMM_WORLD.Reduce(&means, &meansFull, numberMarkers, MPI_FLOAT,MPI_SUM,0);     // MPI
    if (myRank==0) {
	cout<< "meansFull \n " << meansFull.head(1)<<endl;
	}

    MPI::Finalize();     
    return 0;
}

