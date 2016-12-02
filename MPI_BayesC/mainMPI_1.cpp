//  main.cpp
//  BayesCPi (MPI Parallel Computing)

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

float calcVariance(Eigen::VectorXf &var_y){ //computes the variance
    float sum_y = 0.0f, ss_y = 0.0f;
    sum_y = var_y.sum();
    ss_y  = var_y.squaredNorm();
    float mean_y = sum_y/var_y.size();
    float variance_y = ss_y/var_y.size() - mean_y*mean_y;
	if(variance_y < 0)variance_y = 0.0f;
    return variance_y;
}

float calcCovariance(Eigen::VectorXf &var_y, Eigen::VectorXf &var_x){ //computes the covariance
    float sum_y = 0.0f, sum_x=0.0f, scp = 0.0f;
    sum_y = var_y.sum();
    sum_x = var_x.sum();
    scp  = var_y.dot(var_x);
    float mean_y = sum_y/var_y.size();
    float mean_x = sum_x/var_y.size();
    float covariance_x_y = scp/var_y.size() - mean_y*mean_x;
    return covariance_x_y;
}

int main(int argc, const char * argv[]){

    matvec::SESSION.mtr.seed(31415);

    MPI::Init();
    int myRank=MPI::COMM_WORLD.Get_rank();
    int clusterSize=MPI::COMM_WORLD.Get_size();


    const unsigned numberMarkers = 5000;
    const unsigned numObs        = 1000000;
    const unsigned noGeno        = 100;

    unsigned totalSize = numObs;
    unsigned batchSize = totalSize/clusterSize;
    unsigned iStart = myRank*batchSize;
    unsigned iEnd   = ( (myRank+1)==clusterSize ) ? totalSize :  iStart + batchSize;
    unsigned iSize  = iEnd - iStart;

    Eigen::MatrixXf ZFull(1,numberMarkers);
    Eigen::VectorXf yFull(1);
    Eigen::VectorXf a; // for simulated breeding values
    Eigen::MatrixXf Z(iSize,numberMarkers);
    Eigen::VectorXf y(iSize);
    MVectorXf Zi(&Z(0,0),iSize);

    Eigen::VectorXf alpha;
    Eigen::VectorXf meanAlpha;
    Eigen::VectorXf iRes;
    Eigen::VectorXf meaniRes;
    Eigen::VectorXi delta;
    Eigen::VectorXf zpz(numberMarkers);


    unsigned chainLength = 100;
    unsigned outFreq     = 10;
    float  probFixed     = 0.00;
    float nuAlpha        = 4.0;
    float nuRes          = 4.0;
    float nuGen          = 100;
    float nuTilde1;
    float nuTilde2;
    float varGen;
    float varRes;

    float meanVarRes    = 0.0;
    float meanVarAlpha  = 0.0;
    float meanMu        = 0.0;
    float mu;
    float oldNewAlphai[2];
    float mean2pq;
    float varAlpha;
    float scaleAlpha;
    float scaleRes;
    float scaleGen;
    float logPi;
    float logPiComp;

    clock_t iter0Time;
    clock_t iteriTime;


    if (myRank==0) {
        ZFull = Eigen::MatrixXf::Random(numObs,numberMarkers);
//        ZFull.resize(numObs,numberMarkers);
//        string fileName = "/Users/rohan/latex/projects/AFRI2012/BigPed/medGen.bin";
//        std::fstream binaryFile;
//        Eigen::VectorXf genotypes(numberMarkers);
//        unsigned markerLength = numberMarkers*sizeof(genotypes(0));
//        //check for and open the binary genotype file
//        binaryFile.open(fileName.c_str(),  std::ios::in | std::ios::binary);
//        if(!binaryFile) {
//            std::cout << "ERROR: Couldn't open the binary file containing the genotypes " <<fileName << std::endl;
//            return(0);
//        }
//        unsigned count=0;
//        cout << "reading in Z" << endl;
//        //while (binaryFile.read((char*) &genotypes(0), markerLength)) {
//        while (count<numObs) {
//            binaryFile.read((char*) &genotypes(0), markerLength);
//            ZFull.row(count++) = genotypes;
//        }
//        binaryFile.close();
        Eigen::VectorXf Alpha = Eigen::ArrayXf::Random(numberMarkers);
        // simulate y
        a = ZFull*Alpha;
        varGen        = calcVariance(a);
        yFull = a + (Eigen::VectorXf::Random(numObs).array()*sqrt(varGen)).matrix();
        cout << yFull.head(10) << endl;

        varRes        = calcVariance(yFull) - varGen;

        cout << "heritability = " << varGen/(varGen+varRes) << endl;


        Eigen::VectorXf means = ZFull.colwise().mean();
        Eigen::VectorXf p = means/2.0;
        // center Z
        ZFull.rowwise()      -= p.transpose();
        // prior stuff
//        mean2pq    = (2*p.array()*(1-p.array())).mean();
        mean2pq    = 0.5;
        cout      << "mean2pq  = " << mean2pq   << endl;
        zpz = ZFull.colwise().squaredNorm();
        mu  = yFull.mean();
    }
    MPI::COMM_WORLD.Bcast(&mu,1,MPI_FLOAT,0);                                  // MPI broadcast mu
    MPI::COMM_WORLD.Bcast(&mean2pq,1,MPI_FLOAT,0);                             // MPI broadcast mmean2pq
    MPI::COMM_WORLD.Bcast(&varRes,1,MPI_FLOAT,0);                              // MPI broadcast varRes
    MPI::COMM_WORLD.Bcast(&varGen,1,MPI_FLOAT,0);                              // MPI broadcast varGen
    MPI::COMM_WORLD.Bcast(&zpz(0),numberMarkers,MPI_FLOAT,0);                  // MPI broadcast zpz
    varAlpha   = varGen/(numberMarkers*mean2pq);
    scaleAlpha = varAlpha*(nuAlpha-2.0)/nuAlpha;
    scaleRes   = varRes*(nuRes-2.0)/nuRes;
    nuTilde1    = nuRes+float(numObs);
    scaleGen   = varGen*(nuGen-2.0)/nuGen;
    logPi      = log(probFixed);
    logPiComp  = log(1.0-probFixed);

    meanVarRes    = 0.0;
    meanVarAlpha  = 0.0;
    meanMu        = 0.0;

    meanAlpha = Eigen::VectorXf::Zero(numberMarkers);
    meaniRes  = Eigen::VectorXf::Zero(noGeno);
    alpha     = Eigen::VectorXf::Zero(numberMarkers);
    iRes.resize(noGeno);
    delta     = Eigen::VectorXi::Zero(numberMarkers);
    for (unsigned i=0;i<numberMarkers; i++){
        MPI::COMM_WORLD.Scatter(&ZFull(0,i),iSize,MPI_FLOAT,
                                &Z(0,i)    ,iSize,MPI_FLOAT,
                                0);                                       // MPI scatter column i of Z
    }
    //MPI::COMM_WORLD.Barrier();                                            // MPI
    ZFull.resize(0,0);
    MPI::COMM_WORLD.Scatter(&yFull(0),iSize,MPI_FLOAT,
                            &y(0)    ,iSize,MPI_FLOAT,
                            0);                                           // MPI scatter y
    Eigen::VectorXf ycorr = y - Eigen::VectorXf::Constant(iSize,mu);
    // loop for drawing mcmc samples
    iter0Time = time(0);
    for (unsigned iter=1; iter <= chainLength; iter++){
        // sample residual variance
        float ssqi = ycorr.squaredNorm();
        float ssq;
        MPI::COMM_WORLD.Allreduce(&ssqi, &ssq, 1, MPI_FLOAT,MPI_SUM);     // MPI
        // at this point ssq is the sum of squares of all the residuals
        float sampleVal = 2.0f*matvec::sgamma(nuTilde1/2.0f);
        varRes = (ssq + nuRes*scaleRes)/sampleVal;
        meanVarRes += varRes;
        // sample intercept
        float oldMu = mu;
        float rhsi = ycorr.sum();
        float rhs;
        MPI::COMM_WORLD.Allreduce(&rhsi, &rhs, 1, MPI_FLOAT,MPI_SUM);      // MPI
        // rhs now has 1'ycorr from all observations
        rhs += numObs*mu; // unadjusting the rhs for mu
        float invLhs = 1.0/numObs;
        float mean   = rhs*invLhs;
        mu = mean + matvec::snorm()*sqrt(invLhs*varRes);
        meanMu += mu;
        ycorr.array() += (oldMu - mu);
        // sample delta and then the effect for each locus
        unsigned countLoci = 0;
        for (unsigned i=0; i<numberMarkers; i++) {
            new (&Zi) MVectorXf(&Z(0,i),iSize);
            rhsi = Zi.dot(ycorr);
            MPI::COMM_WORLD.Allreduce(&rhsi, &rhs, 1, MPI_FLOAT, MPI_SUM);   // MPI
            // rhs now has Zi'ycorr from all observations
            float oldAlpha = alpha(i);
            float zpzi = zpz(i);
            rhs += zpzi*alpha(i); // unadjusting the rhs for alpha[i]
            float v0 = zpzi*varRes;
            float v1 = zpzi*zpzi*varAlpha + v0;
            float logDelta0 = -0.5*(log(v0) + rhs*rhs/v0) + logPi;
            float logDelta1 = -0.5*(log(v1) + rhs*rhs/v1) + logPiComp;
            float probDelta1 = 1.0/(1.0 + exp(logDelta0-logDelta1));
            if(matvec::ranf()<probDelta1) {
                countLoci++;
                float mmeLhs = zpzi + varRes/varAlpha;
                float invLhs =  1.0/mmeLhs;
                float mean = invLhs*rhs;
                alpha(i) = mean + matvec::snorm()*sqrt(invLhs*varRes);
                meanAlpha(i) += alpha(i);
                delta(i) = 1;
                ycorr   += Zi*(oldAlpha - alpha(i));
            }
            else {
                alpha(i) = 0.0;
                delta(i) = 0;
            }

        }
        cout << "my rank, number of loci = "<< myRank <<", "<<countLoci << endl;
        // sample varAlpha
        nuTilde2  = nuAlpha + float(countLoci);
        sampleVal = 2.0f*matvec::sgamma(nuTilde2/2.0f);
        varAlpha        = (alpha.squaredNorm() + nuAlpha*scaleAlpha)/sampleVal;
        meanVarAlpha   += varAlpha;
        if (!(iter%10)) {
            Eigen::VectorXf aHat = (Z*meanAlpha).array()/float(iter);
            Eigen::VectorXf aHatFull(numObs);
            MPI::COMM_WORLD.Gather(&aHat(0),iSize,MPI_FLOAT,
                                   &aHatFull(0),iSize,MPI_FLOAT,0);
            if (myRank==0) {
                iteriTime = time(0);
                clock_t timePerIter = (iteriTime - iter0Time);
                cout << "iteration: " << iter;
                float var1 = calcVariance(a);
                float var2 = calcVariance(aHatFull);
                float cov  = calcCovariance(a, aHatFull);
                float cor = cov/sqrt(var1*var2);
                cout << " r = " << cor << " time in seconds " << timePerIter << endl;
            }
        }
    }

    MPI::Finalize();
    return 0;
}
