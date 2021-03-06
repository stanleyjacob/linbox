/* tests/test-qlup.C
 * Copyright (C) The LinBox group
 *
 * Time-stamp: <13 Nov 17 16:57:58 Jean-Guillaume.Dumas@imag.fr>
 * -----------------------------------------------------
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


/*! @file  tests/test-qlup.C
 * @ingroup tests
 * @brief  tests LQUP decomposition, solve, and nullspace of a random sparse matrice.
 * @test tests LQUP decomposition, solve, and nullspace of a random sparse matrice.
 */



#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>

#include <cstdio>

#include <linbox/matrix/sparse-matrix.h>
#include "linbox/algorithms/gauss.h"
#include "linbox/algorithms/gauss-gf2.h"
#include "linbox/blackbox/permutation.h"
#include "linbox/util/commentator.h"
#include <givaro/modular.h>
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/direct-sum.h"
#include "linbox/solutions/rank.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: LQUP decomposition of a random sparse matrix
 *
 * Constructs a random sparse matrix and computes its QLUP decomposition
 * using Sparse Gaussian elimination. Checks that the results match.
 */
template <class Field, class Blackbox, class RandStream >
bool testQLUP(const Field &F, size_t n, unsigned int iterations, int rseed, double sparsity = 0.05)
{
	bool res = true;

	commentator().start ("Testing Sparse elimination qlup", "testQLUP", iterations);

	size_t Ni = n;
	size_t Nj = n;
	integer card; F.cardinality(card);
	typename Field::RandIter generator (F,card,rseed);
	RandStream stream (F, generator, sparsity, n, n);

	for (size_t i = 0; i < iterations; ++i) {
		commentator().startIteration ((unsigned)i);


		stream.reset();

		Blackbox A (F, stream);

		std::ostream & report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

		F.write( report ) << endl;
		A.write( report,Tag::FileFormat::Maple ) << endl;

		DenseVector<Field> u(F,Nj), v(F,Ni), w1(F,Nj), w2(F,Ni), w3(F,Ni), w(F,Ni);
		for(auto it=u.begin();it!=u.end();++it)
			generator.random (*it);


		A.apply(v,u);


		unsigned long rank;

		Method::SparseElimination SE;
		SE.strategy(Specifier::PIVOT_LINEAR);
		GaussDomain<Field> GD ( F );
		typename Field::Element determinant;
		Blackbox L(F, A.rowdim(), A.coldim());
		Permutation<Field> Q(F,(int)A.rowdim());
		Permutation<Field> P(F,(int)A.coldim());

		GD.QLUPin(rank, determinant,
			  Q, L, A, P,
			  A.rowdim(), A.coldim() );

		Q.apply(w, L.apply(w3, A.apply(w2, P.apply(w1,u) ) ) );

		bool error = false;
		auto itv=v.begin();
		auto itw=w.begin();
		for( ; itw!=w.end();++itw,++itv) {
			if (! F.areEqual(*itw,*itv) ) {
				error = true;
			}
		}

		if (error) {
			res = false;

			report << "ERROR : matrix(" << u.size() << ",1,[";
			for(auto itu=u.begin(); itu!=u.end();++itu)
				report << *itu << ',';
			report << "]);\n[";
			for(auto itv2=v.begin(); itv2!=v.end();++itv2)
				report << *itv2 << ' ';
			report << "]  !=  [";
			for(auto itw2=w.begin(); itw2!=w.end();++itw2)
				report << *itw2 << ' ';
			report << "]" << std::endl;


			report << "w1: [";
			for(auto itw2=w1.begin(); itw2!=w1.end();++itw2)
				report << *itw2 << ' ';
			report << "]" << std::endl;
			report << "w2: [";
			for(auto itw2=w2.begin(); itw2!=w2.end();++itw2)
				report << *itw2 << ' ';
			report << "]" << std::endl;
			report << "w3: [";
			for(auto itw2=w3.begin(); itw2!=w3.end();++itw2)
				report << *itw2 << ' ';
			report << "]" << std::endl;
		}

		commentator().stop ("done");
		commentator().progress ();
	}

	commentator().stop (MSG_STATUS (res), (const char *) 0, "testQLUP");

	return res;
}

/* Test 2: LQUP solve of a random sparse matrix and a random dense vector
 *
 * Constructs a random sparse matrix and computes its QLUP decomposition
 * using Sparse Gaussian elimination.
 * Then solve using the decomposition and checks that the results match.
 */
template <class Field, class Blackbox, class RandStream>
bool testQLUPsolve(const Field &F, size_t n, unsigned int iterations, int rseed, double sparsity = 0.05)
{
	bool res = true;

	commentator().start ("Testing Sparse elimination qlup solve", "testQLUPsolve", iterations);

	size_t Ni = n;
	size_t Nj = n;
	integer card; F.cardinality(card);
	typename Field::RandIter generator (F,card,rseed);
	RandStream stream (F, generator, sparsity, n, n);

	GF2 F2;
	GF2::RandIter bitgenerator(F2,2,rseed);
	// GF2::Element randomsolve;

	for (size_t i = 0; i < iterations; ++i) {
		commentator().startIteration ((unsigned)i);

		stream.reset();
		Blackbox A (F, stream);

		std::ostream & report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

		F.write( report ) << endl;
		A.write( report, Tag::FileFormat::Maple ) << endl;

		DenseVector<Field> u(F,Nj), v(F,Ni), x(F,Nj), y(F,Ni);
		for(auto it=u.begin();it!=u.end();++it)
			generator.random (*it);

		A.apply(v,u);


		Method::SparseElimination SE;
		SE.strategy(Specifier::PIVOT_LINEAR);
		GaussDomain<Field> GD ( F );

		Blackbox CopyA ( A );

		GD.solvein(x, A, v /*, bitgenerator .random(randomsolve) */ );
		// report << "Random solving: " << randomsolve << std::endl;

		CopyA.apply(y, x);

		VectorDomain<Field> VD(F);


		if (! VD.areEqual(v,y)) {
			res=false;
			A.write( report, Tag::FileFormat::Maple ) << endl;

			report << "ERROR v: matrix(" << v.size() << ",1,[";
			for(auto itu=v.begin(); itu!=v.end();++itu)
				report << *itu << ',';
			report << "]);\n[";
			report << "ERROR y: matrix(" << y.size() << ",1,[";
			for(auto itu=y.begin(); itu!=y.end();++itu)
				report << *itu << ',';
			report << "]);\n[";
			for(auto itv=x.begin(); itv!=x.end();++itv)
				report << *itv << ' ';
			report << "]  !=  [";
			for(auto itw=y.begin(); itw!=y.end();++itw)
				report << *itw << ' ';
			report << "]" << std::endl;

		}

		commentator().stop ("done");
		commentator().progress ();
	}

	commentator().stop (MSG_STATUS (res), (const char *) 0, "testQLUPsolve");

	return res;
}

/* Test 2: LQUP nullspacebasis of a random sparse matrix
 *
 * Constructs a random sparse matrix and computes its QLUP decomposition
 * using Sparse Gaussian elimination (stores only U and P).
 * Then solve using the decomposition and checks that the results match.
 */
template <class Field, class Blackbox, class RandStream>
bool testQLUPnullspace(const Field &F, size_t n, unsigned int iterations, int rseed, double sparsity = 0.05)
{
	bool res = true;

	commentator().start ("Testing Sparse elimination qlup nullspacebasis", "testQLUPnullspace", iterations);

	size_t Ni = n;
	size_t Nj = n;
	integer card; F.cardinality(card);
	typename Field::RandIter generator (F,card,rseed);
	RandStream stream (F, generator, sparsity, n, n, rseed);

	for (size_t i = 0; i < iterations; ++i) {
		commentator().startIteration ((unsigned)i);

		stream.reset();
		Blackbox A (F, stream);

		std::ostream & report = commentator().report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

		F.write( report ) << endl;
		A.write( report, Tag::FileFormat::Maple ) << endl;


		Method::SparseElimination SE;
		SE.strategy(Specifier::PIVOT_LINEAR);
		GaussDomain<Field> GD ( F );

		Blackbox CopyA ( A );
		Blackbox X(F, A.coldim(), A.coldim() );

		GD.nullspacebasisin(X, CopyA );

		size_t nullity = X.coldim();

		DenseVector<Field> u(F,nullity);
		for(auto it=u.begin();it!=u.end();++it)
			generator.random (*it);
		DenseVector<Field> v(F,Nj);
		X.apply(v,u);
		report << "Random combination of the rows of the NullSpace basis" << std::endl;

		DenseVector<Field> w(F,Ni);
		A.apply(w, v);

		VectorDomain<Field> VD(F);

		if (! VD.isZero(w)) {
			res=false;
			A.write( report, Tag::FileFormat::Maple ) << endl;

			report << "ERROR u: matrix(" << u.size() << ",1,[";
			for(auto itu=u.begin(); itu!=u.end();++itu)
				report << *itu << ',';
			report << "]);\n[";
			report << "ERROR v: matrix(" << v.size() << ",1,[";
			for(auto itu=v.begin(); itu!=v.end();++itu)
				report << *itu << ',';
			report << "]);\n[";
			for(auto itv=w.begin(); itv!=w.end();++itv)
				report << *itv << ' ';
			report << "]  !=  0" << std::endl;

		}

		commentator().stop ("done");
		commentator().progress ();
	}

	commentator().stop (MSG_STATUS (res), (const char *) 0, "testQLUPnullspace");

	return res;
}

#define STOR_T SparseMatrixFormat::SparseSeq
// #define STOR_T Vector<Field>::SparseSeq
// #define STOR_T Sparse_Vector<Field::Element>
int main (int argc, char **argv)
{

	commentator().setMaxDepth (-1);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	// 	commentator().setMaxDetailLevel( 100000 );
	// 	commentator().setMaxDepth( 100000 );

	bool pass = true;

	static size_t n = 80;
	static integer q = 65519U;
	static integer bigQ("1234567890123456789012345678901234568123");
	//static integer q = 1000003U;
	static unsigned int iterations = 2;
	static double sparsity = 0.05;
	static int rseed = (int)time(NULL);

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 's', "-s S", "Sparse matrices with density S.", TYPE_DOUBLE,     &sparsity },
		{ 'r', "-r R", "Random generator seed.", TYPE_INT,     &rseed },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);
	srand ((unsigned int)rseed);

	commentator().start("QLUP  test suite", "qlup");
	commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
	<< "Seed: " << rseed << endl;

	{
		commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "over Givaro::Modular<uint32_t,uint64_t>" << endl;
		typedef Givaro::Modular<uint32_t,uint64_t> Field;
		Field F (q);
		typedef SparseMatrix<Field, STOR_T > Blackbox;
		typedef RandomSparseStream<Field, Blackbox::Row   > RandStream;
		if (!testQLUP<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity))
			pass = false;
		if (!testQLUPsolve<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity))
			pass = false;
		if (!testQLUPnullspace<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity))
			pass = false;
	}

	{
		commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "over Givaro::Modular<double>" << endl;
		typedef Givaro::Modular<double> Field;
		Field F (q);
		typedef SparseMatrix<Field, STOR_T > Blackbox;
		typedef RandomSparseStream<Field, Blackbox::Row   > RandStream;

		if (!testQLUP<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity))
			pass = false;
		if (!testQLUPsolve<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity))
			pass = false;
		if (!testQLUPnullspace<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity))
			pass = false;
	}

// 	{

// 		commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
// 		<< "over Givaro::Modular<Integer>" << endl;
// 		typedef Givaro::Modular<Integer> Field;
// 		Field F (bigQ);
// 		typedef SparseMatrix<Field, STOR_T > Blackbox;
// 		typedef RandomSparseStream<Field, Blackbox::Row   > RandStream;

// 		if (!testQLUP<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity))
// 			pass = false;
// 		if (!testQLUPsolve<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity))
// 			pass = false;
// 		if (!testQLUPnullspace<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity))
// 			pass = false;
// 	}

#if 1
	{
		commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
		<< "specialized over GF2>" << endl;
		typedef GF2 Field;
		Field F2;
		typedef LinBox::GaussDomain<LinBox::GF2>::Matrix Blackbox;
		typedef RandomSparseStreamGF2<Blackbox::Row_t> RandStream;
		if (!testQLUP<Field, Blackbox, RandStream> (F2, n, iterations, rseed, sparsity))
			pass = false;
		if (!testQLUPsolve<Field, Blackbox, RandStream> (F2, n, iterations, rseed, sparsity))
			pass = false;
	}
#endif

	commentator().stop(MSG_STATUS (pass),"QLUP test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
