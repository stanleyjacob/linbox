#include <cstdlib>
#include <string>
#include <iostream>
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h"
#include "linbox/algorithms/polynomial-matrix/order-basis.h"

#define TIMINGS_ON // printing timings
#define VERBOSE_ON // printing more information (input specs, output degrees..)
//#define EXTRA_VERBOSE_ON // printing even more information (output matrices..)
//#define EXTRA_TIMINGS_ON // printing extra timings
#define CHECK_ON // checking output

using namespace LinBox;
using namespace std;

/* ~~ SOME USEFUL FUNCTIONS ~~ */
template<typename PolMat>
void print_degree_matrix( const PolMat &pmat )
{
	const size_t d = pmat.degree();
	for ( size_t i=0; i<pmat.rowdim(); ++i ) {
		for ( size_t j=0; j<pmat.coldim(); ++j ) {
			int deg = d;
			while ( deg>=0 and pmat.get(i,j,deg) == 0 )
				--deg;
			cout << deg << "  ";
		}
		cout << endl;
	}
}

template <typename PolMat>
bool test_order( const PolMat &approx, const PolMat &series, const size_t order )
{
	PolynomialMatrixMulDomain<typename PolMat::Field> PMD(approx.field());
	PolMat prod( approx.field(), approx.rowdim(), series.coldim(), approx.size()+series.size() );
	PMD.mul( prod, approx, series );

	bool test = true;
	size_t d = 0;
	while ( test && d<order )
	{
		for ( size_t i=0; i<prod.rowdim(); ++i )
			for ( size_t j=0; j<prod.coldim(); ++j )
			{
				if ( prod.ref(i,j,d) != 0 )
				{
					test = false;
					//std::cout << d << "\t" << i << "\t" << j << std::endl;
				}
			}
		++d;
	}
	return test;
}

template <typename PolMat>
bool test_kernel( const PolMat &kerbas, const PolMat &pmat )
{
	PolynomialMatrixMulDomain<typename PolMat::Field> PMD(kerbas.field());
	PolMat prod( kerbas.field(), kerbas.rowdim(), pmat.coldim(), kerbas.size()+pmat.size()-1 );
	PMD.mul( prod, kerbas, pmat );

	bool test = true;
	size_t d = 0;
	while ( test && d<kerbas.size()+pmat.size()-1 )
	{
		for ( size_t i=0; i<prod.rowdim(); ++i )
			for ( size_t j=0; j<prod.coldim(); ++j )
			{
				if ( prod.ref(i,j,d) != 0 )
				{
					test = false;
					//cout << "degree " << d << endl;
				}
			}
		++d;
	}
	return test;
}

/* ~~ CREATE SHIFT OF GIVEN SHAPE ~~ */
vector<uint64_t> create_shift(
		size_t length,
		size_t shift_shape )
{
	// TODO build shift according to shift_shape
	vector<uint64_t> shift(length,0);
	return shift;
}

vector<int> create_shift(
		size_t length,
		size_t shift_shape,
		int min_value )
{
	// TODO build shift according to shift_shape
	vector<int> shift(length,min_value);
	return shift;
}

/* ~~ BENCHMARK APPROXIMANT BASIS ALGOS ~~ */
template<typename Field, typename RandIter>
void bench_approximant_basis(
		const Field & GF,
		RandIter& rand,
		size_t m,
		size_t n,
		size_t d,
		size_t func,
		size_t shift_shape,
		size_t threshold,
		bool resUpdate ) {

	typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field> MatrixP;
	typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field> PMatrix;

#ifdef VERBOSE_ON
	switch (func) {
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
		case 5:
			cout << "~~~ Benchmarking algo POPOV_PMBASIS ~~~" << endl;
			break;
		default:
			cout << "!!! Invalid function provided through argument '-f' !!!" << endl;
			cout << "(should be among: 0:all; 1:mbasis; 2:newmbasis; 3:pmbasis; 4:newpmbasis; 5:popov_pmbasis.)" << endl;
			cout << "Exiting..." << endl;
			return;
	}
	cout << "--> input polynomial matrix: dimensions " << m << " x " << n << " and degree " << d << endl;
#endif

	// compute order basis
	OrderBasis<Field> AppBas(GF);

	Timer chrono;

	if (func==0 || func==1){ // M-Basis
		vector<size_t> shift = create_shift(m,shift_shape);
		MatrixP *sys = new MatrixP(GF, m, n, d);	
		// set the Serie at random
		for (size_t k=0;k<d;++k)
			for (size_t i=0;i<m;++i)
				for (size_t j=0;j<n;++j)
					rand.random(sys->ref(i,j,k));
#ifdef VERBOSE_ON
		if ( m < 65 )
			cout << "--> input shift:\n" << shift << endl;
		else
			cout << "--> input shift: shape " << shift_shape << endl;
#endif // VERBOSE_ON
#ifdef EXTRA_VERBOSE_ON
		if ( m >= 65 )
			cout << "--> input shift:\n" << shift << endl;
		cout << "--> degrees in input matrix:" << endl;
		print_degree_matrix(*sys);
#endif // EXTRA_VERBOSE_ON
		chrono.start(); // time the matrix creation
		MatrixP appbas(GF, m, m, d+1);
		chrono.stop();
#ifdef EXTRA_TIMINGS_ON
		cout << "TIME (M-Basis), creating basis: " << chrono.usertime() << " s" << endl;
#endif // EXTRA_TIMINGS_ON
		chrono.clear(); chrono.start(); // time the approximant basis computation
		AppBas.M_Basis(appbas, *sys, d, shift);
		chrono.stop();
#ifdef TIMINGS_ON
		cout << "TIME (M-Basis), computing basis: " << chrono.usertime() << " s" << endl;
#endif
#ifdef EXTRA_VERBOSE_ON
		cout << "--> degrees in output basis:" << endl;
		print_degree_matrix(appbas);
#endif // EXTRA_VERBOSE_ON
#ifdef CHECK_ON
		cout << "--> output basis has correct order:" << test_order( appbas, *sys, d ) << endl;
#endif // CHECK_ON
		delete sys;
	}

	if (func==0 || func==2){ // PM-Basis
		vector<size_t> shift = create_shift(m,shift_shape);
		MatrixP *sys = new MatrixP(GF, m, n, d);	
		// set the Serie at random
		for (size_t k=0;k<d;++k)
			for (size_t i=0;i<m;++i)
				for (size_t j=0;j<n;++j)
					rand.random(sys->ref(i,j,k));
#ifdef VERBOSE_ON
		if ( m < 65 )
			cout << "--> input shift:\n" << shift << endl;
		else
			cout << "--> input shift: shape " << shift_shape << endl;
#endif // VERBOSE_ON
#ifdef EXTRA_VERBOSE_ON
		if ( m >= 65 )
			cout << "--> input shift:\n" << shift << endl;
		cout << "--> degrees in input matrix:" << endl;
		print_degree_matrix(*sys);
#endif // EXTRA_VERBOSE_ON
		chrono.clear(); chrono.start(); // time the matrix creation
		MatrixP appbas(GF, m, m, d+1);
		chrono.stop();
#ifdef EXTRA_TIMINGS_ON
		cout << "TIME (PM-Basis), creating basis: " << chrono.usertime() << " s" << endl;
#endif // EXTRA_TIMINGS_ON
		chrono.clear(); chrono.start(); // time the approximant basis computation
		AppBas.PM_Basis(appbas, *sys, d, shift);
		chrono.stop();
#ifdef TIMINGS_ON
		cout << "TIME (PM-Basis), computing basis: " << chrono.usertime() << " s" << endl;
#endif
#ifdef EXTRA_VERBOSE_ON
		cout << "--> degrees in output basis:" << endl;
		print_degree_matrix(appbas);
#endif // EXTRA_VERBOSE_ON
#ifdef CHECK_ON
		cout << "--> output basis has correct order:" << test_order( appbas, *sys, d ) << endl;
#endif // CHECK_ON
		delete sys;
	}

	if (func==0 || func==3){ // new M-Basis
		vector<int> shift = create_shift(m,shift_shape,0);
		PMatrix *sys = new PMatrix(GF, m, n, d);	
		// set the Serie at random
		for (size_t k=0;k<d;++k)
			for (size_t i=0;i<m;++i)
				for (size_t j=0;j<n;++j)
					rand.random(sys->ref(i,j,k));
#ifdef VERBOSE_ON
		if ( m < 65 )
			cout << "--> input shift:\n" << shift << endl;
		else
			cout << "--> input shift: shape " << shift_shape << endl;
#endif // VERBOSE_ON
#ifdef EXTRA_VERBOSE_ON
		if ( m >= 65 )
			cout << "--> input shift:\n" << shift << endl;
		cout << "--> degrees in input matrix:" << endl;
		print_degree_matrix(*sys);
#endif // EXTRA_VERBOSE_ON
		chrono.clear(); chrono.start(); // time the matrix creation
		PMatrix appbas(GF, m, m, d+1);
		chrono.stop();
#ifdef EXTRA_TIMINGS_ON
		cout << "TIME (mbasis), creating basis: " << chrono.usertime() << " s" << endl;
#endif // EXTRA_TIMINGS_ON
		chrono.clear(); chrono.start(); // time the approximant basis computation
		AppBas.mbasis(appbas, *sys, d, shift, resUpdate);
		chrono.stop();
#ifdef TIMINGS_ON
		cout << "TIME (mbasis), computing basis: " << chrono.usertime() << " s" << endl;
#endif
#ifdef EXTRA_VERBOSE_ON
		cout << "--> degrees in output basis:" << endl;
		print_degree_matrix(appbas);
#endif // EXTRA_VERBOSE_ON
#ifdef CHECK_ON
		cout << "--> output basis has correct order:" << test_order( appbas, *sys, d ) << endl;
#endif // CHECK_ON
		delete sys;
	}

	if (func==0 || func==4){ // new PM-Basis
		vector<int> shift = create_shift(m,shift_shape,0);
		PMatrix *sys = new PMatrix(GF, m, n, d);	
		// set the Serie at random
		for (size_t k=0;k<d;++k)
			for (size_t i=0;i<m;++i)
				for (size_t j=0;j<n;++j)
					rand.random(sys->ref(i,j,k));
#ifdef VERBOSE_ON
		if ( m < 65 )
			cout << "--> input shift:\n" << shift << endl;
		else
			cout << "--> input shift: shape " << shift_shape << endl;
#endif // VERBOSE_ON
#ifdef EXTRA_VERBOSE_ON
		if ( m >= 65 )
			cout << "--> input shift:\n" << shift << endl;
		cout << "--> degrees in input matrix:" << endl;
		print_degree_matrix(*sys);
#endif // EXTRA_VERBOSE_ON
		chrono.clear(); chrono.start(); // time the matrix creation
		PMatrix appbas(GF, m, m, d+1);
		chrono.stop();
#ifdef EXTRA_TIMINGS_ON
		cout << "TIME (pmbasis), creating basis: " << chrono.usertime() << " s" << endl;
#endif // EXTRA_TIMINGS_ON
		chrono.clear(); chrono.start(); // time the approximant basis computation
		AppBas.pmbasis(appbas, *sys, d, shift, threshold);
		chrono.stop();
#ifdef TIMINGS_ON
		cout << "TIME (pmbasis), computing basis: " << chrono.usertime() << " s" << endl;
#endif
#ifdef EXTRA_VERBOSE_ON
		cout << "--> degrees in output basis:" << endl;
		print_degree_matrix(appbas);
#endif // EXTRA_VERBOSE_ON
#ifdef CHECK_ON
		cout << "--> output basis has correct order:" << test_order( appbas, *sys, d ) << endl;
#endif // CHECK_ON
		delete sys;
	}

	if (func==0 || func==5){ // Popov PM-Basis
		vector<int> shift = create_shift(m,shift_shape,0);
		PMatrix *sys = new PMatrix(GF, m, n, d);	
		// set the Serie at random
		for (size_t k=0;k<d;++k)
			for (size_t i=0;i<m;++i)
				for (size_t j=0;j<n;++j)
					rand.random(sys->ref(i,j,k));
		chrono.clear(); chrono.start(); // time the matrix creation
#ifdef VERBOSE_ON
		if ( m < 65 )
			cout << "--> input shift:\n" << shift << endl;
		else
			cout << "--> input shift: shape " << shift_shape << endl;
#endif // VERBOSE_ON
#ifdef EXTRA_VERBOSE_ON
		if ( m >= 65 )
			cout << "--> input shift:\n" << shift << endl;
		cout << "--> degrees in input matrix:" << endl;
		print_degree_matrix(*sys);
#endif // EXTRA_VERBOSE_ON
		PMatrix appbas(GF, m, m, d+1);
		chrono.stop();
#ifdef EXTRA_TIMINGS_ON
		cout << "TIME (popov_pmbasis), creating basis: " << chrono.usertime() << " s" << endl;
#endif // EXTRA_TIMINGS_ON
		chrono.clear(); chrono.start(); // time the approximant basis computation
		AppBas.popov_pmbasis(appbas, *sys, d, shift, threshold);
		chrono.stop();
#ifdef TIMINGS_ON
		cout << "TIME (popov_pmbasis), computing basis: " << chrono.usertime() << " s" << endl;
#endif
#ifdef EXTRA_VERBOSE_ON
		cout << "--> degrees in output basis:" << endl;
		print_degree_matrix(appbas);
#endif // EXTRA_VERBOSE_ON
#ifdef CHECK_ON
		cout << "--> output basis has correct order:" << test_order( appbas, *sys, d ) << endl;
#endif // CHECK_ON
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
	size_t threshold=0; // threshold for pm-basis
	size_t func=0; // which function to benchmark: 0:all; 1:mbasis; 2:pmbasis; 3:newmbasis; 4:newpmbasis; 5:popov_pmbasis
	size_t shift_shape=0; // shape of the shift: 0:uniform; {10,11,...} same as {0,1,...} but with large minimum
	bool r=true; // whether to continuously update the residual, or recompute it
	long seed = time(NULL);

	static Argument args[] = {
		{ 'm', "-m m", "Set row dimension of matrix series to m.", TYPE_INT, &m },
		{ 'n', "-n n", "Set column dimension of matrix series to n.", TYPE_INT, &n },
		{ 'd', "-d d", "Set degree of  matrix series to d.", TYPE_INT, &d },
		{ 'p', "-p p", "Set cardinality of prime field (if 0 --> will be random)", TYPE_INT, &p },
		{ 'b', "-b b", "Set bitsize of the matrix entries (if random prime)", TYPE_INT, &b },
		{ 'f', "-f f", "Set the targeted function to benchmark (0:mbasis; 1:mbasis; 2:pmbasis; 3:newmbasis; 4:newpmbasis; 5:popov_pmbasis.", TYPE_INT, &func },
		{ 't', "-t t", "Set degree threshold for pm-basis", TYPE_INT, &threshold },
		{ 'r', "-r r", "Set resUpdate for mbasis", TYPE_BOOL, &r },
		{ 's', "-s s", "Set the shape of the shift ", TYPE_INT, &shift_shape},
		{ 'S', "-S S", "Set the random seed to a specific value", TYPE_INT, &seed},
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);
	size_t logd=integer((uint64_t)d).bitsize();

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
		bench_approximant_basis(GF,rand,m,n,d,func,shift_shape,threshold,r);
	}
	else if (b < 26){ // here, logd<=b-4
		RandomFFTPrime Rd(1<<b,seed);
		integer p = Rd.randomPrime(logd+1);
#ifdef VERBOSE_ON
		cout << "--> using random FFT prime with " << b << " bits, p = " << p << endl;
#endif
		SmallField GF(p);
		typename SmallField::RandIter rand(GF,0,seed);
		bench_approximant_basis(GF,rand,m,n,d,func,shift_shape,threshold,r);
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
		bench_approximant_basis(GF,G,m,n,d,func,shift_shape,threshold,r);
	}
	
	return 0;
}
