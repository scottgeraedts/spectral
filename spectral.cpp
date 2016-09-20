#include "spectral.h"
//load eigenvectors from nicolas file
void Spectral::load(double h){
	stringstream filename;
	char buffer[300];
	sprintf(buffer,"spin_1_2_closedchain_n_16_jx_1.000000_jy_1.000000_jz_1.000000_randomhx_%8.6f__randomhy_0.000000_randomhz_%8.6f_runid_0.eigenvec.mat",h,h);
	filename<<buffer;
	cout<<filename.str()<<endl;
	return;
	
	ifstream File;
	File.open(filename.str().c_str(), ios::binary | ios::in);
	int TmpType;
	int TmpNbrRow;
	int TmpNbrColumn;
	File.read ((char*) &(TmpType), sizeof(int));
	File.read ((char*) &(TmpNbrRow), sizeof(int));
	File.read ((char*) &(TmpNbrColumn), sizeof(int));
	if ((TmpType & 1) != 0)
	{
	double Tmp;
	for (int i = 0; i < TmpNbrRow; ++i)
		for (int j = 0; j < TmpNbrColumn; ++j)
		{
			File.read ((char*) &(Tmp), sizeof(double));
			eigvecs(i,j)=Tmp;
		}
	}
	else
	 {
	 	cout<<"trying to read complex eigenvectors"<<endl;
	 	exit(0);
	 }
//	   double TmpRe;
//	   double TmpIm;
//	   for (int i = 0; i < TmpNbrRow; ++i)
//		 for (int j = 0; j < TmpNbrColumn; ++j)
//		   {
//		     File.read ((char*) &(TmpRe), sizeof(double));
//		     File.read ((char*) &(TmpIm), sizeof(double));
//		     cout << i << " " << j << " " << TmpRe << " " << TmpIm <<endl;
//		   }
//	 }
}


//  Matrix-vector multiplication w <- M*v.
void Spectral::HeisenbergXZ(double* v, double* w, double J){
	int countJ=0;
	int sign;
	double countH;
	for(int in=0; in<rows; in++) w[in]=0.;
	for(int in=0; in<rows; in++){
		//off-diagonal elements
		countH=0.;
		for(int i=0; i<N; i++){
			if(bittest(states[in],i)) sign=1;
			else sign=-1;

			w[lookup_flipped(in, states, 1, i)]+=0.5*alphax[i];
			//hz
			countH+=0.5*sign*alphaz[i];
		}
		
		//diagonal elements
		countJ=0;
		for(int i=0; i<N;i++){
			if(bittest(states[in],i) != bittest(states[in],next(i)) ){
				w[ lookup_flipped(in , states ,2, i,next(i))  ]+=2*v[in]*0.25;
				countJ--;
			}else
				countJ++;	
		}
		w[in]+=(countJ*J*0.25+countH)*v[in];	
	}
} //  MultMv.

//  Matrix-vector multiplication w <- M*v.
void Spectral::HeisenbergZ(double* v, double* w, double J){
	int countJ=0;
	int sign;
	double countH;
	for(int in=0; in<rows; in++) w[in]=0.;
	for(int in=0; in<rows; in++){
		//off-diagonal elements
		countH=0.;
		for(int i=0; i<N; i++){
			if(bittest(states[in],i)) sign=1;
			else sign=-1;

			//hz
			countH+=0.5*sign*alphaz[i];
		}
		
		//diagonal elements
		countJ=0;
		for(int i=0; i<N;i++){
			if(bittest(states[in],i) != bittest(states[in],next(i)) ){
				w[ lookup_flipped(in , states ,2, i,next(i))  ]+=2*v[in]*0.25;
				countJ--;
			}else
				countJ++;	
		}
		w[in]+=(countJ*J*0.25+countH)*v[in];	
	}
} //  MultMv.

void Spectral::Sz(double* v, double* w, int site){
	int sign;
	for(int in=0; in<rows; in++) w[in]=0.;
	for(int in=0; in<rows; in++){
		if(!bittest(states[in],site)) sign=1;
		else sign=-1;
		w[in]+=sign*v[in];	
	}
}
void Spectral::make_states(int charge=-1){
	states.clear();
	for(int i=0;i<1<<N;i++)
		if(charge==-1 || count_bits(i)==charge) states.push_back(i);
	rows=states.size();
}
void Spectral::makeDense( function<void(double  *v, double *w)> matvec, Eigen::MatrixXd &EigenDense){
    EigenDense=Eigen::Matrix<double,-1,-1>::Zero(rows,rows);
//	dense=new ART[n*n];
	double *v=new double[rows];
	double *w=new double[rows];
	for(int i=0;i<rows;i++){
		for(int j=0; j<rows; j++){
			if(i==j) v[j]=1;
			else v[j]=0;
			w[j]=0;
		}
		matvec(v,w);
		for(int j=0; j<rows; j++){
//			dense[i+j*n]=w[j];
			EigenDense(j,i)=w[j];
		}
	}
	delete [] v;
	delete [] w;
}
void Spectral::makeSparse( function<void(double  *v, double *w)> matvec, Eigen::SparseMatrix<double> &EigenSparse){
    EigenSparse=Eigen::SparseMatrix<double>(rows,rows);
    vector<Eigen::Triplet<double> > triplets;
//	dense=new ART[n*n];
	double *v=new double[rows];
	double *w=new double[rows];
	for(int i=0;i<rows;i++){
		for(int j=0; j<rows; j++){
			if(i==j) v[j]=1;
			else v[j]=0;
			w[j]=0;
		}
		matvec(v,w);
		for(int j=0; j<rows; j++){
			triplets.push_back(Eigen::Triplet<double>(j,i,w[j]));
		}
	}
	delete [] v;
	delete [] w;
	EigenSparse.setFromTriplets(triplets.begin(),triplets.end());
}
int Spectral::next(int i){
	if(i%N==N-1) return 0;
	else return i+1;
}

void Spectral::check(Eigen::VectorXd eigvals){
	for(int n=0;n<rows;n++){
		if( abs(eigvecs.col(n).adjoint()*Hnn*eigvecs.col(n) - eigvals(n)) > 1e-10){
			cout<<"check failed!"<<n<<" "<<eigvecs.col(n).adjoint()*Hnn*eigvecs.col(n)<<" "<<eigvals(n)<<endl;
		}
	}
}
Spectral::Spectral(){
	ifstream cin("params");
	int seed;
	cin>>N;
	cin>>h;
	cin>>loading;
	cin>>seed;

	MTRand ran(seed);
	alphax=vector<double>(N);
	alphaz=vector<double>(N);
	
	make_states(N/2);

	//make spin operators
	vector<Eigen::SparseMatrix<double> > szs(N);
//	vector<Eigen::MatrixXd > szs(N);
	for(int i=0;i<N;i++) makeSparse(bind(&Spectral::Sz,this,placeholders::_1,placeholders::_2,i),szs[i]);

	//if dense, make Hamiltonian
	//if sparse, do nothing
	if(loading){
		
	}
	else{
		for(int i=0;i<N;i++){
			alphax[i]=(2*ran.rand()-1)*h;
			alphaz[i]=(2*ran.rand()-1)*h;
		}
		makeDense(bind(&Spectral::HeisenbergZ,this,placeholders::_1,placeholders::_2,1.),Hnn); 
	}
	
	//put energies into eigvals, eigenvectors into eigvals
	if(loading){
		load(h);
	}
	else{
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Hnn);
		eigvecs=es.eigenvectors();				
		
		check(es.eigenvalues());
		eigvals=es.eigenvalues();
	}

	//make the bins on Ebar, omega and expectation values
	//I'm trying to use an adaptive grid for these so that I don't wind up with a lot of empty bins
	
	//Ebar
	//the adaptive grid strategy is: use 25 bins for the top and bottom quarters, 100 bins for the middle half
	double Ebarmin=eigvals(0);
	double Ediff=eigvals(eigvals.size()-2)-Ebarmin;
	vector<double> Ebar_bins;
	for(double e=Ebarmin;e<Ebarmin+Ediff/4.;e+=Ediff/4./25.) Ebar_bins.push_back(e);
	for(double e=Ebarmin+Ediff/4.;e<Ebarmin+3*Ediff/4.;e+=Ediff/2./100.) Ebar_bins.push_back(e);
	for(double e=Ebarmin+3*Ediff/4.;e<Ebarmin+Ediff;e+=Ediff/4./25.) Ebar_bins.push_back(e);
	Ebar_bins.push_back(Ebarmin+Ediff);
	
	//omega
	//the adaptive grid strategy is: use 100 bins for the first quarter, 100 bins for the last 3/4
	vector<double> omega_bins;
	for(double e=0;e<Ediff/4.;e+=Ediff/4./100.) omega_bins.push_back(e);
	for(double e=Ediff/4.;e<Ediff;e+=3*Ediff/4./100.) omega_bins.push_back(e);
	omega_bins.push_back(Ediff);
	
	//M
	//here there's no variable widths since its too much of a headach
	vector<double> M_bins;
	int nbins=1000;
	for(int e=0;e<nbins;e++) M_bins.push_back(2*e/(1.*nbins)-1.);
	M_bins.push_back(1);
//	for(double e=-1;e<-0.5;e+=0.05) M_bins.push_back(e);
//	for(double e=-0.5;e<-0.1;e+=0.008) M_bins.push_back(e);
//	for(double e=-0.1;e<0;e+=0.002) M_bins.push_back(e);
//	for(double e=0;e<0.1;e+=0.002) M_bins.push_back(e);
//	for(double e=0.1;e<0.5;e+=0.008) M_bins.push_back(e);
//	for(double e=0.5;e<=1;e+=0.05) M_bins.push_back(e);
	
	vector< vector< vector<int> > > data(Ebar_bins.size(), vector< vector<int> >(omega_bins.size(), vector<int>(M_bins.size(),0) ) );
	
	ofstream sout("spectra");
	ofstream lout("labels");
	
	//compute matrix elements
	vector<Eigen::MatrixXd> out(N);
	for(int i=0;i<N;i++) out[i]=eigvecs.adjoint()*szs[i]*eigvecs;
	
	//bin computed elements
	vector< vector<double> > vals(N),  unfolded_vals(N);
	vector<double>::iterator Ebar_it,omega_it,M_it;
	cout<<"limits: "<<eigvals(0)<<" "<<eigvals(eigvals.size()-2)<<endl;
	cout<<"sizes: "<<Ebar_bins.size()<<" "<<omega_bins.size()<<" "<<M_bins.size()<<endl;
	for(int m=0;m<rows-1;m++){
		for(int n=0;n<m;n++){ 
			Ebar_it=lower_bound(Ebar_bins.begin(),Ebar_bins.end(),0.5*(eigvals(n)+eigvals(m)));
			omega_it=lower_bound(omega_bins.begin(),omega_bins.end(),eigvals(m)-eigvals(n));
//			sout<<eigvals(m)<<" "<<eigvals(n)<<" ";			
			for(int i=0;i<N;i++){
				if( abs(out[i](m,n))>1e-6){ 
					M_it=lower_bound(M_bins.begin(),M_bins.end(),out[i](m,n));
//					cout<<"\t"<<eigvals(n)<<" "<<eigvals(m)<<endl;
//					cout<<Ebar_it-Ebar_bins.begin()<<" "<<omega_it-omega_bins.begin()<<" "<<M_it-M_bins.begin()<<endl;
					data[Ebar_it-Ebar_bins.begin()][omega_it-omega_bins.begin()][M_it-M_bins.begin()]++;
				}
//				sout<<out[i](m,n)<<" ";
				
			}
//			sout<<endl;
		}
	}
	
	//print binned elements
	lout<<Ebar_bins.size()<<" "<<omega_bins.size()<<endl;
	for(int i=0;i<(signed)Ebar_bins.size();i++){
		for(int j=0;j<(signed)omega_bins.size();j++){
			lout<<Ebar_bins[i]<<" "<<omega_bins[j]<<endl;
		}
	}

	for(int k=0;k<(signed)M_bins.size();k++){
		sout<<M_bins[k]<<" ";
		for(int i=0;i<(signed)Ebar_bins.size();i++){
			for(int j=0;j<(signed)omega_bins.size();j++){
				sout<<data[i][j][k]<<" ";
			}
		}
		sout<<endl;
	}
	sout.close();
	lout.close();
}	
int main(){
	Spectral s;
}		
//		Eigen::MatrixXd EH=(  (-1.*es.eigenvalues().array()).exp()  );
//		Eigen::DiagonalMatrix<double,-1,-1> EHD=Eigen::DiagonalMatrix<double,-1,-1> (EH);
//		Eigen::MatrixXd expH=eigvecs*EHD*eigvecs.adjoint();
	
//		cout<<"Z="<<EH.sum()<<endl;
//		for(int i=0;i<N;i++)
//			cout<<(szs[i]*expH).trace()<<" ";
//		cout<<endl;


//		Eigen::VectorXd eigvals;
//		for(int i=0;i<N;i++){
//			cout<<out[i]<<endl;
//			es.compute(out[i]);
//			eigvals=es.eigenvalues();
//			cout<<es.eigenvalues()<<endl;
//			Eigen_To_Std(eigvals,vals[i]);
//			sort(vals[i].begin(),vals[i].end());
//			unfolded_vals[i]=unfoldE(vals[i],50);
//		}
//		
//		sout.open("unfolded");
//		for(int i=0;i<unfolded_vals[0].size();i++){
//			for(int j=0;j<N;j++){
//				sout<<unfolded_vals[j][i]<<" ";
//			}sout<<endl;
//		} 
//		sout.close();

