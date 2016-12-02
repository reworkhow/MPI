//  ZPZ full-stored, single-site BayesC
#include <iostream>
#include <Eigen>
#include <string>
#include <iostream>
#include <fstream>
#include <matvec/statdist.h>
//#include <matvec/session.h>

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

    //matvec::SESSION.mtr.seed(12345);
    srand (31415);

    clock_t iter0Time;
    unsigned numObs          = 1000;
    unsigned numberMarkers   = 5000;


    Eigen::MatrixXf Z = Eigen::MatrixXf::Random(numObs,numberMarkers);
    Eigen::VectorXf Alpha = Eigen::VectorXf::Random(numberMarkers);
    unsigned chainLength = 1000;
    float  probFixed     = 0.0;
    float  nuAlpha       = 4.1;
    float  nuRes         = 100;
    float  nuGen         = 100;
    unsigned outFreq     = 10;


    Eigen::VectorXf a = Z*Alpha;
    Eigen::VectorXf y = a + Eigen::VectorXf::Random(numObs)*2;
    cout << y.head(10) << endl;
    float varGen        = calcVariance(a);
    float varRes        = calcVariance(y) - varGen;
    cout << "probFixed = " << probFixed << endl;
    cout << "varRes = " << varRes << endl;
    float mean2pq       = 0.5;
    cout << "mean2pq  = " << mean2pq   << endl;
    float varAlpha    =  calcVariance(Alpha);
    cout << "varAlpha = " << varAlpha << endl;
    cout << "heritability = " << varAlpha/(varAlpha + varRes) << endl;
    float scaleAlpha = varAlpha*(nuAlpha-2.0)/nuAlpha;
    float scaleRes   = varRes*(nuRes-2.0)/nuRes;
    float scaleGen   = varGen*(nuGen-2.0)/nuGen;
    float logPi      = log(probFixed);
    float logPiComp  = log(1.0-probFixed);

    float meanVarRes    = 0.0;
    float meanVarAlpha  = 0.0;
    float meanMu        = 0.0;

    Eigen::VectorXf alpha       = Eigen::VectorXf::Zero(numberMarkers);
    Eigen::VectorXf meanAlpha   = Eigen::VectorXf::Zero(numberMarkers);
    Eigen::VectorXf iRes(numberMarkers);
    Eigen::VectorXf meaniRes(numberMarkers);
    Eigen::VectorXi delta       = Eigen::VectorXi::Zero(numberMarkers);

    cout << "computing ZPZ " << endl;
    Eigen::MatrixXf C = Z.transpose() * Z;


    Eigen::VectorXf rhs = Z.transpose() * y;
    MVectorXf Ci(&C(0,0),numberMarkers);
    float ypy = y.squaredNorm();

    cout << "start sampling " << endl;
    iter0Time = time(0);
     // loop for drawing mcmc samples
    for (unsigned iter=1; iter<=chainLength; iter++){
        float invIter   = 1.0f/iter;
        float invVarRes = 1.0f/varRes;
        float lambda   = varRes/varAlpha;
        float logVarAlpha = log(varAlpha);
        float logDeltaZero = logPi;
        for (unsigned i=0; i<numberMarkers; i++) {
            new (&Ci) MVectorXf(&C(0,i),numberMarkers);
            float lhs = (C(i,i) + lambda);//problem before
            float invLhs = 1.0f/lhs;
            float mean = (rhs(i) - Ci.dot(alpha) + alpha(i)*C(i,i))*invLhs;
            float logDeltaOne = -0.5f*(log(lhs*invVarRes) + logVarAlpha - mean*rhs(i)*invVarRes) + logPiComp;//problem
            float probDeltaOne = 1.0f/(1.0f + exp(logDeltaZero-logDeltaOne));
            if(matvec::ranf()<probDeltaOne) {
                delta[i] = 1;
                alpha[i] = mean + matvec::snorm()*sqrt(invLhs*varRes);//
                meanAlpha[i] += (alpha[i] - meanAlpha[i])*invIter;
            }
            else {
                delta[i] = 0;
                alpha[i] = 0.0f;
            }
        }
        // sample varAlpha
        float nuTilde   = nuAlpha + delta.sum();
        float sampleVal = 2.0f*matvec::sgamma(nuTilde/2.0f);
        varAlpha        = (alpha.squaredNorm() + nuAlpha*scaleAlpha)/sampleVal;
        meanVarAlpha   += (varAlpha - meanVarAlpha)*invIter;

        // sample residual variance
        float ssq = ypy - rhs.dot(alpha)*2 + alpha.dot(C*alpha);
        nuTilde   = nuRes+ float(numObs);
        sampleVal = 2.0f*matvec::sgamma(nuTilde/2.0f);
        varRes = (ssq + nuRes*scaleRes)/sampleVal;
        meanVarRes += (varRes - meanVarRes)*invIter;
        if (!(iter%10)) {
            Eigen::VectorXf aHat = (Z*meanAlpha).array()/float(iter);
            float var1 = calcVariance(a);
            float var2 = calcVariance(aHat);
            float cov  = calcCovariance(a, aHat);
            cout <<var1 << " " << var2 << " " << cov << endl;
            float cor = cov/sqrt(var1*var2);
            clock_t timeTotal = (time(0) - iter0Time);
            cout << "iteration: " << iter << " r = " << cor << " time: " << timeTotal << endl;
            cout << "delta.sum() " << delta.sum() << " varRes " << " " << meanVarRes << " varAlpha " << " " << meanVarAlpha << endl;
        }
    }
    return 0;
}
