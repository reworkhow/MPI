#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <stat.h>
#include <safe_vectors.h>
#include <doublematrix.h>
#include <map>
#include <mim.h>
#include <mpi.h>
#include <vector>
#include <mpi.h>
#include <Eigen/Eigen>
#include <boost/format.hpp>


class HapVal{
public:
 static Vector4d p;
 static Vector4d x;
 static int Ne;
 int h00, h01, h10;
 double pr;
 void getProb();
};

class MultinomialDist {
public:
	static double getProb(Vector4d x, Vector4d p){
		double logP;
	    unsigned n = x.size();
		double pr = matvec::gammln(double(HapVal::Ne + 1));
		for (unsigned i=0;i<n;i++){
			if(p[i]==0){
				logP = (x[i]==0) ? 0.0 : -10000.0;
			}
			else {
				logP = log(p[i]);
			}
			pr += logP*x[i] - matvec::gammln(double(x[i] + 1));
		}
		if (pr > 0.0) return 1.0;
		if (pr < -600.00) return 0.0;
		return std::exp(pr);
	}
};
void HapVal::getProb() {
	x[0] = h00;
	x[1] = h01;
	x[2] = h10;
	x[3] = Ne - h00 - h01 - h10;
	if (x[3]<0){
		cout << "h11 = " << x[3] << endl;
		exit (1);
	}
	pr = MultinomialDist::getProb(x,p);
 }


class RDist{
public:
    int myRank, clusterSize;
    SafeSTLVector<double> dist, newDist;
	SafeSTLVector<HapVal> multDist;
	double r, mut,coeffSel, maf1, maf2;;
	unsigned genStart;
    
	void initialize(){
		HapVal::p.resize(4);
		HapVal::p[0] = 0.25;
		HapVal::p[1] = 0.25;
		HapVal::p[2] = 0.25;
		HapVal::p[3] = 0.25;
		HapVal::x.resize(4);
		HapVal hapVal;
		multDist.name = "multDist";
		dist.name = "dist";
		for (unsigned i=0;i<=HapVal::Ne;i++){
			for (unsigned j=0;j<=HapVal::Ne-i;j++){
				for (unsigned k=0;k<=HapVal::Ne-i-j;k++){
					hapVal.h00 = i;
					hapVal.h01 = j;
					hapVal.h10 = k;
					hapVal.getProb();
					multDist.push_back(hapVal);     //a vector composed of objects for hapVal
					dist.push_back(hapVal.pr);      //construct a vector "dist" composed of pr for different haplotype combinations
				}
			}
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
//			sum11 += h11/double(HapVal::Ne)*dist[i];
//			cout << multDist[i].h00 << "\t"
//			     << multDist[i].h01 << "\t"
//				 << multDist[i].h10 << "\t"
//				 <<             h11 << "\t"
//				 << dist[i] << endl;
				 double p1 = multDist[i].h10/double(HapVal::Ne) + h11/double(HapVal::Ne);
				 double p2 = multDist[i].h01/double(HapVal::Ne) + h11/double(HapVal::Ne);
//				 mapit = pDist.find(c1);
//				 if (mapit!=pDist.end()) {
//					pDist[c1] += dist[i];
//				 }
//				 else {
//				    pDist[c1] = dist[i];
//				 }	
			if(p1>c1&&p1<c2 && p2>c3&&p2<c4){                
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
            
//				int r21000 = r2*1000;
//				double r2Rounded = r21000/1000.0;
//				r2Mapit = r2Dist.find(r2Rounded);
//				if (r2Mapit!=r2Dist.end()) {
//					r2Dist[r2Rounded] += dist[i];
//				}
//				else {
//					r2Dist[r2Rounded] = dist[i];
//				}	
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
//		cout << "Gen = "<< setw(10) << setiosflags (ios::right | ios::fixed) << gen                     << "  ";
//		cout << setw(15) << setprecision (10)  << setiosflags (ios::right | ios::fixed) << relClass1    << "  ";
//		cout << setw(15) << setprecision (10)  << setiosflags (ios::right | ios::fixed) << relClass2    << "  ";
//		cout << setw(15) << setprecision (10)  << setiosflags (ios::right | ios::fixed) << sumClass3    << "  ";
//		cout << setw(15) << setprecision (10)  << setiosflags (ios::right | ios::fixed) << mur2         << "  ";
//		cout << setw(15) << setprecision (10)  << setiosflags (ios::right | ios::fixed) << s2           << "  ";
//		cout << setw(15) << setprecision (10)  << setiosflags (ios::right | ios::fixed) << mup1         << "  ";
//		cout << setw(15) << setprecision (10)  << setiosflags (ios::right | ios::fixed) << mup2         << "\n";
        
        cout << "E(r2) = " << mur2 << endl;
//		for (mapit=pDist.begin();mapit!=pDist.end();mapit++){
//			cout << mapit->first << "\t" << mapit->second << endl;
//		}
//		for (r2Mapit=r2Dist.begin();r2Mapit!=r2Dist.end();r2Mapit++){
//			cout << setw(16) << setprecision (14) << r2Mapit->first << "\t" << r2Mapit->second/sum << endl;
//		}

	}
    
    
	void calcMultDist(){
		for (unsigned i=0;i<multDist.size();i++){
			multDist[i].getProb();
		}
	}
	void calcDistNextGen(){
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
		
		newDist.resize(multDist.size(),0.0); 
		for (unsigned i=0;i<multDist.size();i++) newDist[i] = 0.0;
		double sum = 0.0;
		SafeSTLVector<double> sumDist;
		sumDist.resize(multDist.size(),0.0);
		
		unsigned totalSize = sumDist.size();
		unsigned batchSize = totalSize/clusterSize;
		unsigned iStart = myRank*batchSize;
		unsigned iEnd   = ( (myRank+1)==clusterSize ) ? totalSize :  iStart + batchSize; //for last cluster,from istart to the totalsize
		
		for (unsigned i=iStart;i<iEnd;i++){
			int h11 = HapVal::Ne - multDist[i].h00 - multDist[i].h01 - multDist[i].h10;
			p[0] = multDist[i].h00/double(HapVal::Ne);
			p[1] = multDist[i].h01/double(HapVal::Ne);
			p[2] = multDist[i].h10/double(HapVal::Ne);
			p[3] =             h11/double(HapVal::Ne);
			
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
			
			double wbar = (1-coeffSel)*(p00 + p01) + p10 + p11;//selection?
			p[0] = (1-coeffSel)*p00/wbar;                      
			p[1] = (1-coeffSel)*p01/wbar;
			p[2] = p10/wbar;
			p[3] = p11/wbar;
			
			
			HapVal::p = p;
			
			
			calcMultDist();
			for (unsigned j=0;j<multDist.size();j++){
				newDist[j] += multDist[j].pr*dist[i];
			}
		}
		if (clusterSize > 1) {
			double *sendbuf = &newDist[0];
			double *recvbuf = &sumDist[0];
			MPI::COMM_WORLD.Allreduce(sendbuf, recvbuf, totalSize, MPI_DOUBLE, MPI_SUM);
		}
		else {
			sumDist = newDist;
		}
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
			outFile << dist[i] << "\n";
		}
	}
};

int HapVal::Ne = 200;                                                 // effective pop size is Ne/2;
Vector4d HapVal::p, HapVal::x;

int main (int argc, char * const argv[]) {
	
	RDist rDist;
	rDist.clusterSize = 1;
	rDist.myRank = 0;
	rDist.openMPI();
    
    unsigned ngen   = 2000;

 
	string inFileName = "";
    if (inFileName!=""){
	    rDist.initialize();
		rDist.input(inFileName);

	}
	else {
		rDist.r     = atof(argv[1]);
		rDist.mut   = atof(argv[2]);
        rDist.coeffSel  = 0.00;
        rDist.maf1 = rDist.maf2  = 0.0;
        unsigned ngen   = atoi(argv[3]);

		rDist.initialize();
		rDist.display();
		rDist.genStart = 1;
	}
	cout << "Ne = "         << HapVal::Ne/2.0
         << "\t recomb = "  << rDist.r
         << "\t Mu = "      << rDist.mut
         << "\t S= "        << rDist.coeffSel
         << "\t maf1= "     << rDist.maf1
         << "\t maf1= "     << rDist.maf2   << endl;
	
    unsigned i;
	for (i=rDist.genStart+1;i<=rDist.genStart+ngen;i++){
		rDist.calcDistNextGen();
		//if (rDist.myRank==0) cout << "Gen = " << i << " ";
		rDist.display();
        if ((i%100)==0) rDist.save("Rdist.save",i+1);
	}
	rDist.save("RDist.save",i);
	MPI::Finalize();
    return 0;
}
