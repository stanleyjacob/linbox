/*
 * examples/minpoly.C
 *
 * Copyright (C) 2005, 2010 D Saunders, J-G Dumas
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

/** \file examples/minpoly.C
 * @example  examples/minpoly.C
  \brief Minimal polynomial of a sparse matrix.
  \ingroup examples
  */
#include <linbox/linbox-config.h>
#include <iostream>
//! @bug this should be elsewhere
template <class Field, class Polynomial>
void printPolynomial (const Field &F, const Polynomial &v)
{
	for (int i = (int)v.size () ; i-- ; ) {
		F.write (std::cout, v[(size_t)i]);
		if (i > 0)
			std::cout << " x^" << i << " + ";
	}
	std::cout << std::endl;
}

#include <linbox/ring/modular.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/blackbox/compose.h>
#include <linbox/solutions/minpoly.h>

using namespace LinBox;
using namespace std;


int main (int argc, char **argv)
{
	commentator().setMaxDetailLevel (-1);
	commentator().setMaxDepth (-1);
	commentator().setReportStream (std::cerr);

	int a = argc;
	while(a--){
		cerr << "argv[" << a << "] = " << argv[a] << endl;
	}
	if (argc < 2) {
		cerr << "Usage: minpoly <matrix-file-in-SMS-format> [<p>]" << endl;
		return -1;
	}

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file " << argv[1] << endl; return -1; }

	if (argc != 3) {
		int process = 0;

		Method::Blackbox M;
#ifdef __LINBOX_HAVE_MPI
		Communicator C(&argc, &argv);
		process = C.rank();
		M.communicatorp(&C);
#endif

		Givaro::ZRing<Integer> ZZ;
		typedef SparseMatrix<Givaro::ZRing<Integer>> SpMat;
		/*
		SpMat A (ZZ);
		A.read (input);
		*/
		SpMat B (ZZ);
		B.read (input);
		typedef Transpose<SpMat> BB2;
		BB2 BT(B);
//		typedef DenseMatrix<Givaro::ZRing<Integer> > BB2;
//		BB2 BT(ZZ, B.coldim(), B.rowdim());
//		BT.random();
		Compose<SpMat, BB2 > A(B,BT);

		if(process == 0)
			cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;

		DenseVector<Givaro::ZRing<Integer> > m_A(ZZ);
		minpoly (m_A, A, M);

		if(process == 0){
			cout << "Minimal Polynomial is ";
			printPolynomial (ZZ, m_A);
		}

		//            minpoly (m_A, A, Method::DenseElimination() );

		//            if(process == 0){
		//                cout << "Minimal Polynomial is ";
		//                printPolynomial (ZZ, m_A);
		//            }

	}
	else{
		typedef Givaro::Modular<double> Field;
		double q = atof(argv[2]);
		Field F(q);
		SparseMatrix<Field> B (F);
		B.read (input);
		cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;

		DenseVector<Field> m_B(F);
		minpoly (m_B, B);

		cout << "Minimal Polynomial is ";
		printPolynomial (F, m_B);

		//                 minpoly (m_A, A, Method::DenseElimination() );
		//  		cout << "Minimal Polynomial is ";
		// 		printPolynomial (F, m_B);

	}

	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
