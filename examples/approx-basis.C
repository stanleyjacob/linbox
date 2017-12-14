#include <cstdlib>
#include <string>
#include <iostream>
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h"
#include "linbox/algorithms/polynomial-matrix/order-basis.h"

#define TIMINGS_ON // printing timings
#define EXTRA_TIMINGS_ON // printing extra timings
#define VERBOSE_ON // printing more information (input specs, output degrees..)
//#define EXTRA_VERBOSE_ON // printing even more information (output matrices..)
//#define CHECK_ON // printing more information (input specs, output degrees..)

using namespace LinBox;
using namespace std;

template<typename Field, typename RandIter>
void bench_approximant_basis(
		const Field & GF,
		RandIter& rand,
		size_t m,
		size_t n,
		size_t d,
		size_t f,
		size_t t ) {
	typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field> MatrixP;

#ifdef EXTRA_VERBOSE_ON
	cout << sys << endl;
#endif // EXTRA_VERBOSE_ON
			
	// define the shift
	vector<int> shift_int(m,0);

	// compute order basis
	OrderBasis<Field> AppBas(GF);

	Timer chrono;

	if (f==0 || f==1){ // M-Basis
		MatrixP *sys = new MatrixP(GF, m, n, d);	
		// set the Serie at random
		for (size_t k=0;k<d;++k)
			for (size_t i=0;i<m;++i)
				for (size_t j=0;j<n;++j)
					rand.random(sys->ref(i,j,k));
		vector<size_t> shift(m,0);
		chrono.start(); // time the matrix creation
		MatrixP appbas(GF, m, m, d+1);
		chrono.stop();
#ifdef EXTRA_TIMINGS_ON
		cout << "TIME (M-Basis), creating basis: " << chrono.usertime() << " s" << endl;
#endif // EXTRA_TIMINGS_ON
		chrono.start(); // time the approximant basis computation
		AppBas.M_Basis(appbas, *sys, d, shift);
		chrono.stop();
#ifdef TIMINGS_ON
		cout << "TIME (M-Basis), computing basis: " << chrono.usertime() << " s" << endl;
#endif
		delete sys;
	}

	if (f==0 || f==2){ // PM-Basis
		MatrixP *sys = new MatrixP(GF, m, n, d);	
		// set the Serie at random
		for (size_t k=0;k<d;++k)
			for (size_t i=0;i<m;++i)
				for (size_t j=0;j<n;++j)
					rand.random(sys->ref(i,j,k));
		vector<size_t> shift(m,0);
		chrono.start(); // time the matrix creation
		MatrixP appbas(GF, m, m, d+1);
		chrono.stop();
#ifdef EXTRA_TIMINGS_ON
		cout << "TIME (PM-Basis), creating basis: " << chrono.usertime() << " s" << endl;
#endif // EXTRA_TIMINGS_ON
		chrono.start(); // time the approximant basis computation
		AppBas.PM_Basis(appbas, *sys, d, shift);
		chrono.stop();
#ifdef TIMINGS_ON
		cout << "TIME (PM-Basis), computing basis: " << chrono.usertime() << " s" << endl;
#endif
		delete sys;
	}


}

int main(int argc, char** argv){
	typedef Givaro::Modular<double> SmallField;
	typedef Givaro::Modular<RecInt::ruint128,RecInt::ruint256> LargeField;

	// default arguments:
	size_t m=64; // input matrix row dimension
	size_t n=32; // input matrix column dimension
	size_t d=32; // matrix degree
	size_t p=0; // size of the base field (default: random)
	size_t b=20; // entries bitsize
	size_t t=0; // threshold for pm-basis
	size_t f=0; // which function to benchmark: 0:mbasis; 3:newmbasis; 2:pmbasis; 3:newpmbasis.
	long seed = time(NULL);

	static Argument args[] = {
		{ 'm', "-m m", "Set row dimension of matrix series to m.", TYPE_INT, &m },
		{ 'n', "-n n", "Set column dimension of matrix series to n.", TYPE_INT, &n },
		{ 'd', "-d d", "Set degree of  matrix series to d.", TYPE_INT, &d },
		{ 'p', "-p p", "Set cardinality of prime field (if 0 --> will be random)", TYPE_INT, &p },
		{ 'b', "-b b", "Set bitsize of the matrix entries (if random prime)", TYPE_INT, &b },
		{ 'f', "-f f", "Set the targeted function to benchmark (0:mbasis; 3:newmbasis; 2:pmbasis; 3:newpmbasis.", TYPE_INT, &f },
		{ 't', "-t t", "Set degree threshold for pm-basis", TYPE_INT, &t },
		{ 's', "-s s", "Set the random seed to a specific value", TYPE_INT, &seed},
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);
	size_t logd=integer((uint64_t)d).bitsize();

#ifdef VERBOSE_ON
	switch (f) {
		case 0:
			cout << "~~~ Benchmarking algos MBASIS, PMBASIS, new MBASIS, new PMBASIS ~~~" << endl;
			break;
		case 1:
			cout << "~~~ Benchmarking algo MBASIS ~~~" << endl;
			break;
		case 2:
			cout << "~~~ Benchmarking algo PMBASIS ~~~" << endl;
			break;
		case 3:
			cout << "~~~ Benchmarking NEW algo MBASIS ~~~" << endl;
			break;
		case 4:
			cout << "~~~ Benchmarking NEW algo PMBASIS ~~~" << endl;
			break;
		default:
			cout << "!!! Invalid function provided through argument '-f' !!!" << endl;
			cout << "(should be among: 0:mbasis; 3:newmbasis; 2:pmbasis; 3:newpmbasis.)" << endl;
			cout << "Exiting..." << endl;
			return 0;
	}
	cout<<"--> input polynomial matrix: dimensions " << m << " x " << n << " and degree " << d << endl;
#endif

	if (p>0 or (b<26 && logd>b-4)) {
		if (p>0) {
#ifdef VERBOSE_ON
			cout << "--> using provided prime p = " << p << endl;
#endif
		}
		else {
			RandomPrimeIterator Rd(b,seed);	
			p = Rd.randomPrime();
#ifdef VERBOSE_ON
			cout << "--> using random prime with " << b << " bits, p = " << p << endl;
			cout << "    (not using FFT prime: degree too large for field bitsize)" << endl;
#endif
		}
		SmallField GF(p);
		typename SmallField::RandIter rand(GF,0,seed);
		bench_approximant_basis(GF,rand,m,n,d,f,t);
	}
	else if (b < 26){ // here, logd<=b-4
		RandomFFTPrime Rd(1<<b,seed);
		integer p = Rd.randomPrime(logd+1);
#ifdef VERBOSE_ON
		cout << "--> using random FFT prime with " << b << " bits, p = " << p << endl;
#endif
		SmallField GF(p);
		typename SmallField::RandIter rand(GF,0,seed);
		bench_approximant_basis(GF,rand,m,n,d,f,t);
	}
	else {
		RandomPrimeIterator Rd(b,seed);	
		integer pp = Rd.randomPrime();
#ifdef VERBOSE_ON
		cout << "--> using random prime with " << b << " bits, p = " << pp << endl;
		cout << "    (large bitsize: multi-precision FFT with Chinese remainderings" << endl;
#endif
		LargeField GF(pp);
		typename LargeField::RandIter G(GF,b,seed);
		bench_approximant_basis(GF,G,m,n,d,f,t);
	}
	
	return 0;
}
