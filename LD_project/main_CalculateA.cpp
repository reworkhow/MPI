//install library? some need, some doesn't.
//mpi
//using boost::format;
//vector map

#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <vector>
#include <mpi.h>
#include <Eigen/Eigen>    
#include <boost/format.hpp> 

using namespace std;
using boost::format;
using namespace Eigen;

class HapVal{
public:
 static Vector4d p;
 static Vector4d x;
 static int Ne;
 int h00, h01, h10;
 double getProb();
};

class MultinomialDist {
public:
	static VectorXd gammalnVec;
	
	static void makeGammalnVec(void ) {
	    int n = HapVal::Ne;
		gammalnVec.resize(n+1);
		gammalnVec[0] = 0.0;
		for (int i=1;i<=n;i++){
			gammalnVec[i] = gammalnVec[i-1] + log(double(i));
		}
	}
	static double getProb(Vector4d x, Vector4d p){
        double logP;
	    unsigned n = x.size();
		double pr = gammalnVec[HapVal::Ne];
		for (unsigned i=0;i<n;i++){
			if(p[i]==0){
				logP = (x[i]==0) ? 0.0 : -10000.0;
			}
			else {
				logP = log(p[i]);
			}
			pr += logP*x[i] - gammalnVec[x[i]];
		}
		if (pr > 0.0) return 1.0;
		if (pr < -600.00) return 0.0;
		return std::exp(pr);
	}
};
double HapVal::getProb() {
	x[0] = h00;
	x[1] = h01;
	x[2] = h10;
	x[3] = Ne - h00 - h01 - h10;
	if (x[3]<0){
		cout << "h11 = " << x[3] << endl;
		exit (1);
	}
	return (MultinomialDist::getProb(x,p));
 }


class RDist{
public:
    int myRank, clusterSize;
    VectorXd dist, newDist;
	MatrixXd A;
	vector<HapVal> multDist;//a vector of several HapVal objects
	double r, mut, coeffSel, maf1, maf2; //MAF -- not in R
	unsigned genStart,gen;

	void initialize(){
		HapVal::p.resize(4);
		HapVal::p[0] = 0.25;
		HapVal::p[1] = 0.25;
		HapVal::p[2] = 0.25;
		HapVal::p[3] = 0.25;
		HapVal::x.resize(4);
		HapVal hapVal;
		for (unsigned i=0;i<=HapVal::Ne;i++){
			for (unsigned j=0;j<=HapVal::Ne-i;j++){
				for (unsigned k=0;k<=HapVal::Ne-i-j;k++){
					hapVal.h00 = i;
					hapVal.h01 = j;
					hapVal.h10 = k;
					multDist.push_back(hapVal); //std::vector::push_back
                                                //Adds a new element at the end of the vector, after its current last element.
				}
			}
		}
		dist.resize(multDist.size());
		for (unsigned i=0;i<multDist.size();i++){
			dist[i] = multDist[i].getProb();
		}									 
	}
	
	void openMPI(void){
		MPI::Init();
		myRank=MPI::COMM_WORLD.Get_rank();
		clusterSize=MPI::COMM_WORLD.Get_size();
	}
	
	void display (){
	    if (myRank) return;
		double mur2 = 0.0, sum = 0.0, sumClass1 = 0.0, sumClass2 = 0.0, sumClass3;
		double mud2 = 0.0, muDenom = 0.0;
		map<int,double> pDist;
		map<int,double>::iterator mapit;
		map<double,double> r2Dist;//rsquare
		map<double,double>::iterator r2Mapit;
		int h11;
		double c1 = maf1;
		double c2 = 1.0 - maf1;
		double c3 = maf2;
		double c4 = 1.0 - maf2;
		double mup1=0.0,mup2=0.0;
        
		for (unsigned i=0;i<dist.size();i++){
			h11 = HapVal::Ne - multDist[i].h00 - multDist[i].h01 - multDist[i].h10;
			double p1 = multDist[i].h10/double(HapVal::Ne) + h11/double(HapVal::Ne);
			double p2 = multDist[i].h01/double(HapVal::Ne) + h11/double(HapVal::Ne);
//			if (p1==0||p1==1||p2==0||p2==1) {
//				continue;
//			}
//			if(p1<c1||p1>c2 && p2<c3||p2>c4){
			if (p1>c1&&p1<c2 && p2>c3&&p2<c4) { 
				double d2 = (h11/double(HapVal::Ne) - p1*p2)*(h11/double(HapVal::Ne) - p1*p2);
				double r2 = d2/(p1*(1-p1)*p2*(1-p2));
				if (r2 < 0.999999) {
					sumClass1 += dist[i]; 
				}
				else if(r2 >= 0.999999){
					sumClass2 += dist[i]; 
				}
				mur2 += r2*dist[i];
				mud2 += d2*dist[i];
				mup1 += p1*dist[i];
				mup2 += p2*dist[i];
				muDenom += p1*(1-p1)*p2*(1-p2)*dist[i];
				sum += dist[i];
			}	
		}
		mur2 /= sum;
		mud2 /= sum;
		muDenom /= sum;
		double s2 = mud2/muDenom;
		mup1 /= sum;
		mup2 /= sum;
		sumClass3 = 1.0 - sumClass1 - sumClass2;
		double sumClass1PlusClass2 = sumClass1 + sumClass2;
		double relClass1 = sumClass1/sumClass1PlusClass2;
		double relClass2 = sumClass2/sumClass1PlusClass2;
		cout << "Gen = "<< setw(10) << setiosflags (ios::right | ios::fixed) << gen                     << "  ";
		cout << setw(15) << setprecision (10)  << setiosflags (ios::right | ios::fixed) << relClass1    << "  "; 
		cout << setw(15) << setprecision (10)  << setiosflags (ios::right | ios::fixed) << relClass2    << "  "; 
		cout << setw(15) << setprecision (10)  << setiosflags (ios::right | ios::fixed) << sumClass3    << "  "; 
		cout << setw(15) << setprecision (10)  << setiosflags (ios::right | ios::fixed) << mur2         << "  "; 
		cout << setw(15) << setprecision (10)  << setiosflags (ios::right | ios::fixed) << s2           << "  "; 
		cout << setw(15) << setprecision (10)  << setiosflags (ios::right | ios::fixed) << mup1         << "  "; 
		cout << setw(15) << setprecision (10)  << setiosflags (ios::right | ios::fixed) << mup2         << "\n";
        
        cout << "E(r2) = " << mur2 << endl;
		
        //		for (mapit=pDist.begin();mapit!=pDist.end();mapit++){
		//			cout << mapit->first << "\t" << mapit->second << endl;
		//		}
		//		for (r2Mapit=r2Dist.begin();r2Mapit!=r2Dist.end();r2Mapit++){
		//			cout << setw(16) << setprecision (14) << r2Mapit->first << "\t" << r2Mapit->second/sum << endl;
		//		}
		
	}
	
	void calcA(){//OKAY!!! get A matrix, the transform matrix
		Vector4d p;			
		
		double a = (1-mut)*(1-mut);
		double b = (1-mut)*mut;
		double c = mut*mut;
		
		Matrix4d Mu;
		Mu <<	a, b, b, c,
		b, a, c, b,
		b, c, a, b,
		c, b, b, a
		;
		
		unsigned totalSize = multDist.size();
		unsigned batchSize = totalSize/clusterSize;
		unsigned iStart = myRank*batchSize;
		unsigned iEnd   = ( (myRank+1)==clusterSize ) ? totalSize :  iStart + batchSize;
		unsigned rows   = multDist.size();
		unsigned cols   = iEnd - iStart;
		A.setConstant(rows, cols,0.0);		
		//cout << "rows = " << rows << " cols " << cols << endl;
		for (unsigned i=iStart;i<iEnd;i++){
			int h11 = HapVal::Ne - multDist[i].h00 - multDist[i].h01 - multDist[i].h10;
			p[0] = multDist[i].h00/double(HapVal::Ne);
			p[1] = multDist[i].h01/double(HapVal::Ne);
			p[2] = multDist[i].h10/double(HapVal::Ne);
			p[3] =         h11/double(HapVal::Ne);
			
			double p00 = p[0]*(1-r) + r/(HapVal::Ne - 1) * 
			(p[0]*(multDist[i].h00 - 1 + multDist[i].h10)+
			 p[1]*(multDist[i].h00     + multDist[i].h10));
			
			double p01 = p[1]*(1-r) + r/(HapVal::Ne - 1) * 
			(p[0]*(multDist[i].h01     + h11)+
			 p[1]*(multDist[i].h01 - 1 + h11));
			
			double p10 = p[2]*(1-r) + r/(HapVal::Ne - 1) * 
			(p[2]*(multDist[i].h10 - 1 + multDist[i].h00)+
			 p[3]*(multDist[i].h10     + multDist[i].h00));
			
			double p11 = p[3]*(1-r) + r/(HapVal::Ne - 1) * 
			(p[2]*(h11     + multDist[i].h01)+
			 p[3]*(h11 - 1 + multDist[i].h01));
			
			p[0] = p00;
			p[1] = p01;
			p[2] = p10;
			p[3] = p11;
			
			p = Mu*p;
			
			p00 = p[0];
			p01 = p[1];
			p10 = p[2];
			p11 = p[3];
			
			double wbar = (1-coeffSel)*(p00 + p01) + p10 + p11;//?what's this?
			p[0] = (1-coeffSel)*p00/wbar;                      //selection?
			p[1] = (1-coeffSel)*p01/wbar;
			p[2] = p10/wbar;
			p[3] = p11/wbar;
			
			
			HapVal::p = p;
			
			for (unsigned j=0;j<multDist.size();j++){
				A(j,i-iStart)  += multDist[j].getProb();
			}
		}
	}
	
	void calcDistNextGen() {
		
		gen++;
		VectorXd sumDist;
		sumDist.setConstant(multDist.size(),0.0);
		
		unsigned totalSize = multDist.size();
		unsigned batchSize = totalSize/clusterSize;
		unsigned iStart = myRank*batchSize;
		unsigned iEnd   = ( (myRank+1)==clusterSize ) ? totalSize :  iStart + batchSize;
		unsigned cols   = iEnd - iStart;
		
		newDist = A*dist.segment(iStart, cols);
		
		if (clusterSize > 1) {
			double *sendbuf = &newDist[0];
			double *recvbuf = &sumDist[0];
			MPI::COMM_WORLD.Allreduce(sendbuf, recvbuf, totalSize, MPI_DOUBLE, MPI_SUM);//??MPI_DOUBLE MPI_SUM??
		}
		else {
			sumDist = newDist;
		}
		double sum = 0.0;
		for (unsigned j=0;j<multDist.size();j++){
			sum += sumDist[j];
		}
		for (unsigned j=0;j<multDist.size();j++){
			dist[j] = sumDist[j]/sum;
		}
	}
	void input(string fileName){
		ifstream inFile;
		inFile.open(fileName.c_str());
		inFile >> HapVal::Ne >> genStart;
		gen = genStart;
		inFile >> r    >> mut; 
		double prob;
		unsigned i = 0;
		while (inFile >> prob){
			dist[i++] = prob;
		}
	}
	void save(string fileName, unsigned genEnd){
		if (myRank) return;
		ofstream outFile;
		outFile.open(fileName.c_str());
		outFile << HapVal::Ne << " " << genEnd-1 << endl;
		outFile << r << " " << mut << endl;
		for (unsigned i=0;i<dist.size();i++){
			outFile << setw(18) << setprecision (12)  << setiosflags (ios::right | ios::fixed) << dist[i] << "\n";
		}
	}
	void saveDist(string fileName){

		if (myRank) return;////????
//		MatrixXd D;
//		D.setConstant(HapVal::Ne+1, HapVal::Ne+1,0.0);
//		ofstream outFile;
//		outFile.open(fileName.c_str());
//		for (unsigned i=0;i<dist.size();i++){
//			int h11 = HapVal::Ne - multDist[i].h00 - multDist[i].h01 - multDist[i].h10;
//			int p1 = multDist[i].h10 + h11;
//			int p2 = multDist[i].h01 + h11;
//			D(p1,p2) += dist[i];
//			
//		}
//		double sum=0.0;
//		for (unsigned i=1; i<HapVal::Ne;i++){
//			for (unsigned j=1; j<HapVal::Ne;j++){
//				sum +=  D(i,j);
//			}
//		}
//		for (unsigned i=1; i<HapVal::Ne;i++){
//			for (unsigned j=1; j<HapVal::Ne;j++){
//					outFile << setw(18) << setprecision (12)  << setiosflags (ios::right | ios::fixed)<< i/double(HapVal::Ne)  << "  "
//					        << setw(18) << setprecision (12)  << setiosflags (ios::right | ios::fixed)<< j/double(HapVal::Ne)  << "  ";
//				    outFile << setw(18) << setprecision (12)  << setiosflags (ios::right | ios::fixed)<< D(i,j)/sum  << "\n";
//			}
//		}
//	}
	
		ofstream outFile;
		outFile.open(fileName.c_str());
		for (unsigned i=0;i<dist.size();i++){
			int h11 = HapVal::Ne - multDist[i].h00 - multDist[i].h01 - multDist[i].h10;
			double p1 = multDist[i].h10/double(HapVal::Ne) + h11/double(HapVal::Ne);
			double p2 = multDist[i].h01/double(HapVal::Ne) + h11/double(HapVal::Ne);
			double maf1 = (p1<0.5) ? p1 : 1.0-p1;
			double maf2 = (p2<0.5) ? p2 : 1.0-p2;
			outFile << setw(18) << setprecision (12)  << setiosflags (ios::right | ios::fixed) << multDist[i].h00  << "  "
			        << setw(18) << setprecision (12)  << setiosflags (ios::right | ios::fixed) << multDist[i].h01  << "  "
			        << setw(18) << setprecision (12)  << setiosflags (ios::right | ios::fixed) << multDist[i].h10  << "  "
			        << setw(18) << setprecision (12)  << setiosflags (ios::right | ios::fixed) << h11  << "  "
			        << setw(18) << setprecision (12)  << setiosflags (ios::right | ios::fixed) << maf1 << "  "
			        << setw(18) << setprecision (12)  << setiosflags (ios::right | ios::fixed) << maf2 << "  ";
			outFile << format("%17.14e \n") % dist[i];
		}
	}
	
};

Vector4d HapVal::p, HapVal::x;
VectorXd MultinomialDist::gammalnVec;

//  Effective pop size is Ne/2;  //
int HapVal::Ne = 5; 
int main (int argc, char *argv[]) {
	
	
	RDist rDist;
	
	//default vallues 
	rDist.r        = 0.002;
	rDist.mut      = 1e-9;
	rDist.coeffSel = 0.00;
	rDist.maf1 = rDist.maf2    = 0.0;
	unsigned ngen  = 2000;
	
	if(argc==7){//???
		rDist.r   = atof(argv[1]); //convert string to double
		rDist.mut = atof(argv[2]);
		rDist.coeffSel = atof(argv[3]);
		rDist.maf1 = atof(argv[4]);
		rDist.maf2 = atof(argv[5]);
		ngen       = atoi(argv[6]);
	}
	
	rDist.clusterSize = 1;
	rDist.myRank = 0;
	rDist.openMPI();
	
	if (rDist.myRank==0) cout << "Ne= " << HapVal::Ne/2.0 
		                      << "\t recomb= " << rDist.r 
		                      << "\t Mu= "     << rDist.mut 
		                      << "\t S= "      << rDist.coeffSel
		                      << "\t maf1= "    << rDist.maf1 
		                      << "\t maf1= "    << rDist.maf2 << endl;

	
	MultinomialDist::makeGammalnVec();
 
	string inFileName = "";//????
	//string inFileName = "RDist.save";
    if (inFileName!=""){
	    rDist.initialize();
		rDist.input(inFileName);
	}
	else {
		rDist.initialize();
		rDist.gen = 1;
		rDist.display();
		rDist.genStart = 1;
	}
	unsigned i;
	rDist.calcA();
	for (i=rDist.genStart+1;i<=rDist.genStart+ngen;i++){
		rDist.calcDistNextGen();
		//if (rDist.myRank==0) cout << "Gen = "<< setw(10) << setiosflags (ios::right | ios::fixed) << i   << " ";
		rDist.display();
		//if ((i%10)==0) rDist.save("Rdist.save",i+1);
	}
	rDist.saveDist("RDist.save");	
	MPI::Finalize();
    return 0;
}
