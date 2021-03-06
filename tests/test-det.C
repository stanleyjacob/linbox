/* tests/test-det.C
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 *
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
 *
 */


/*! @file  tests/test-det.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */



#include "linbox/linbox-config.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <givaro/givrational.h>
#include "linbox/util/commentator.h"
#include "givaro/modular.h"
#include "linbox/vector/blas-vector.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/solutions/det.h"
#include "linbox/solutions/methods.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: Determinant of nonsingular diagonal matrix with distinct entries
 *
 * Construct a random diagonal matrix with distinct entries and see that its
 * computed determinant equals the product of diagonal entries
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of iterations to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDiagonalDet1 (Field &F, size_t n, int iterations)
{
    typedef BlasVector<Field> Vector;
    typedef Diagonal <Field> Blackbox;

    commentator().start ("Testing nonsingular diagonal determinant (1)", "testDiagonalDet1", (unsigned int) iterations);

    bool ret = true;
    bool done;
    int i;
    size_t j, k;

    VectorDomain<Field> VD (F);

    Vector d(F,n);
    typename Field::Element pi, phi_wiedemann, phi_symm_wied, phi_blas_elimination, phi_sparseelim;
    typename Field::RandIter r (F);

    for (i = 0; i < iterations; i++) {
        commentator().startIteration ((unsigned int) i);

        F.assign(pi, F.one);

        for (j = 0; j < n; j++) {
            do {
                do r.random (d[j]); while (F.isZero (d[j]));

                done = true;
                for (k = 0; k < j; k++) {
                    if (F.areEqual (d[j], d[k])) {
                        done = false;
                        break;
                    }
                }
            } while (!done);

            F.mulin (pi, d[j]);
        }

        ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
        report << "Diagonal entries: ";
        VD.write (report, d);
        report << endl;

        report << "True determinant: ";
        F.write (report, pi);
        report << endl;

        Blackbox D (d);

        Method::Wiedemann WiedemannChoice;
        det (phi_wiedemann, D,  WiedemannChoice);
        F.write (report << "Computed determinant (Wiedemann) : ", phi_wiedemann) << endl;

        WiedemannChoice.symmetric(Specifier::SYMMETRIC);
        det (phi_symm_wied, D,  WiedemannChoice);
        F.write (report << "Computed determinant (Symmetric Wiedemann) : ", phi_symm_wied) << endl;

        det (phi_blas_elimination, D,  Method::DenseElimination ());
        F.write (report << "Computed determinant (DenseElimination) : ", phi_blas_elimination) << endl;

        det (phi_sparseelim, D,  Method::SparseElimination ());
        F.write (report << "Computed determinant (SparseElimination) : ", phi_sparseelim) << endl;

        if (!F.areEqual (pi, phi_wiedemann) || !F.areEqual (pi, phi_blas_elimination) || !F.areEqual(pi, phi_symm_wied)|| !F.areEqual(pi, phi_sparseelim)) {
            ret = false;
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: Computed determinant is incorrect" << endl;
        }

        commentator().stop ("done");
        commentator().progress ();
    }

    commentator().stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalDet1");

    return ret;
}

/* Test 2: Determinant of nonsingular diagonal matrix with nondistinct entries
 *
 * Construct a random diagonal matrix with nondistinct entries and see that its
 * computed determinant equals the product of diagonal entries
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of iterations to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDiagonalDet2 (Field &F, size_t n, int iterations)
{
    typedef BlasVector<Field> Vector;
    typedef Diagonal <Field> Blackbox;

    commentator().start ("Testing nonsingular diagonal determinant (2)", "testDiagonalDet2", (unsigned int) iterations);

    bool ret = true;
    int i, k;
    size_t j;

    Vector d(F,n);
    typename Field::Element pi, phi_wiedemann, phi_symm_wied, phi_blas_elimination, phi_sparseelim;
    typename Field::RandIter r (F);

    for (i = 0; i < iterations; i++) {
        commentator().startIteration ((unsigned int)i);

        F.assign(pi, F.one);

        size_t m = (n+1) / 2;
        for (j = 0; j < m; j++) {
            do r.random (d[j]); while (F.isZero (d[j]));
            F.mulin (pi, d[j]);
        }

        for (j = n / 2; j < n; j++) {
            k =int( (size_t)rand () % m );
            d[j] = d[(size_t)k];
            F.mulin (pi, d[(size_t)j]);
        }

            //ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
        ostream &report = commentator().report ();
        report << "Diagonal entries: ";
        printVector<Field> (F, report, d);

        report << "True determinant: ";
        F.write (report, pi);
        report << endl;

        Blackbox D (d);
        Method::Wiedemann WiedemannChoice;
        det (phi_wiedemann, D,  WiedemannChoice);
        WiedemannChoice.symmetric(Specifier::SYMMETRIC);
        report << "Computed determinant (Wiedemann) : ";
        F.write (report, phi_wiedemann);
        report << endl;

        det (phi_symm_wied, D,  WiedemannChoice);
        report << "Computed determinant (Symmetric Wiedemann) : ";
        F.write (report, phi_symm_wied);
        report << endl;

        det (phi_blas_elimination, D,  Method::DenseElimination ());
        report << "Computed determinant (DenseElimination) : ";
        F.write (report, phi_blas_elimination);
        report << endl;

        det (phi_sparseelim, D,  Method::SparseElimination ());
        report << "Computed determinant (SparseElimination) : ";
        F.write (report, phi_sparseelim);
        report << endl;

        if (!F.areEqual (pi, phi_wiedemann) || !F.areEqual (pi, phi_blas_elimination) || !F.areEqual(pi, phi_symm_wied) || !F.areEqual(pi, phi_sparseelim)) {
            ret = false;
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: Computed determinant is incorrect" << endl;
        }

        commentator().stop ("done");
        commentator().progress ();
    }

    commentator().stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalDet2");

    return ret;
}

/* Test 3: Determinant of singular diagonal matrix
 *
 * Construct a random diagonal matrix with one zero entry
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of iterations to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testSingularDiagonalDet (Field &F, size_t n, int iterations)
{
    typedef BlasVector <Field> Vector;
    typedef Diagonal <Field> Blackbox;

    commentator().start ("Testing singular diagonal determinant", "testSingularDiagonalDet",(size_t) iterations);

    bool ret = true;
    int i;
    size_t j;

    Vector d(F,n);
    typename Field::Element phi_wiedemann, phi_symm_wied, phi_blas_elimination, phi_sparseelim;
    typename Field::RandIter r (F);

    for (i = 0; i < iterations; i++) {
        commentator().startIteration ((unsigned int)i);

        for (j = 0; j < n; j++)
            r.random (d[j]);

        F.assign(d[rand () % n], F.zero);

        ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
        report << "Diagonal entries: ";
        printVector<Field> (F, report, d);

        Blackbox D (d);

        Method::Wiedemann WiedemannChoice;
        det (phi_wiedemann, D,  WiedemannChoice);
        report << "Computed determinant (Wiedemann) : ";
        F.write (report, phi_wiedemann);
        report << endl;

        WiedemannChoice.symmetric(Specifier::SYMMETRIC);
        det (phi_symm_wied, D,  WiedemannChoice);
        report << "Computed determinant (Symmetric Wiedemann) : ";
        F.write (report, phi_symm_wied);
        report << endl;

        det (phi_blas_elimination, D,  Method::DenseElimination ());
        report << "Computed determinant (DenseElimination) : ";
        F.write (report, phi_blas_elimination);
        report << endl;

        det (phi_sparseelim, D,  Method::SparseElimination ());
        report << "Computed determinant (SparseElimination) : ";
        F.write (report, phi_sparseelim);
        report << endl;

        if (!F.isZero (phi_wiedemann) || !F.isZero (phi_blas_elimination) || !F.isZero (phi_symm_wied)  || !F.isZero (phi_sparseelim) ) {
            ret = false;
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: Computed determinant is incorrect" << endl;
        }

        commentator().stop ("done");
        commentator().progress ();
    }

    commentator().stop (MSG_STATUS (ret), (const char *) 0, "testSingularDiagonalDet");

    return ret;
}

/* Test 4: Integer determinant
 *
 * Construct a random nonsingular diagonal sparse matrix and compute its
 * determinant over Z
 *
 * n - Dimension to which to make matrix
 * iterations - Number of iterations to run
 *
 * Returns true on success and false on failure
 */

bool testIntegerDet (size_t n, int iterations)
{
    commentator().start ("Testing integer determinant", "testIntegerDet", (unsigned int)iterations);

    bool ret = true;

    for (int i = 0; i < iterations; ++i) {
        commentator().startIteration ((unsigned int)i);
        Givaro::IntegerDom R;
        SparseMatrix<Givaro::IntegerDom> A (R, n, n);

        integer pi = 1;
        integer det_A_wiedemann, det_A_symm_wied, det_A_blas_elimination;

        for (unsigned int j = 0; j < n; ++j) {
            integer &tmp = A.refEntry (j, j);
            integer::nonzerorandom (tmp, 20*i + 1);
            integer::mulin (pi, tmp);
        }

        if (i % 2) {
            integer::negin(A.refEntry(1,1));
            integer::negin(pi);
        }

            //  if (i == iterations - 1) {
            //          A.setEntry(1, 2, A.getEntry(1,1));
            //          A.setEntry(2, 1, A.getEntry(2,2));
            //          pi = 0;
            //      }

            //GMP_Integers R;
            //A.write(cout,R);
        ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

        report << "True determinant: ";
        report << pi;
        report << endl;

        Method::Wiedemann WiedemannChoice;
        det (det_A_wiedemann, A, WiedemannChoice);
        report << "Computed integer determinant (Wiedemann): " << det_A_wiedemann << endl;


        WiedemannChoice.symmetric(Specifier::SYMMETRIC);
        det (det_A_symm_wied, A, WiedemannChoice);
        report << "Computed integer determinant (Symmetric Wiedemann): " << det_A_symm_wied << endl;

        det (det_A_blas_elimination, A, Method::DenseElimination());
        report << "Computed integer determinant (DenseElimination): " << det_A_blas_elimination << endl;


        if ((det_A_wiedemann != pi)||(det_A_blas_elimination != pi)||(det_A_symm_wied != pi))  {
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: Computed determinant is incorrect" << endl;
            ret = false;
        }

        commentator().stop ("done");
        commentator().progress ();
    }

    commentator().stop (MSG_STATUS (ret), (const char *) 0, "testIntegerDet");

    return ret;
}

/* Test 5: Integer determinant by generic methods
 *
 * Construct a random nonsingular diagonal sparse matrix and compute its
 * determinant over Z
 *
 * n - Dimension to which to make matrix
 * iterations - Number of iterations to run
 *
 * Returns true on success and false on failure
 */

bool testIntegerDetGen (size_t n, int iterations)
{
    commentator().start ("Testing integer determinant, generic methods", "testIntegerDeterminantGeneric", (unsigned int) iterations);

    bool ret = true;

    for (int i = 0; i < iterations; ++i) {
        commentator().startIteration ((unsigned int)i);
        Givaro::IntegerDom R;
        SparseMatrix<Givaro::IntegerDom> A (R, n, n);

        integer pi = 1;
        integer det_A, det_A_H, det_A_B, det_A_E;

        for (unsigned int j = 0; j < n; ++j) {
            integer &tmp = A.refEntry (j, j);
            integer::nonzerorandom (tmp, 20*i + 1);
            integer::mulin (pi, tmp);
        }

        if (i % 2) {
            integer::negin(A.refEntry(1,1));
            integer::negin(pi);
        }

        ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

        report << "True determinant: " << pi << endl;

        det (det_A, A);
        report << "Computed integer determinant (Default): " << det_A << endl;
        if (det_A != pi){
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: Default (integer dense) Computed determinant is incorrect" << endl;
            ret = false;
        }

        det (det_A_H, A, Method::Auto());
        report << "Computed integer determinant (Auto): " << det_A_H << endl;
        if (det_A_H != pi){
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: Auto Computed determinant is incorrect" << endl;
            ret = false;
        }

        det (det_A_B, A, Method::Blackbox());
        report << "Computed integer determinant (Blackbox): " << det_A_B << endl;
        if (det_A_B != pi){
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: Blackbox Computed determinant is incorrect" << endl;
            ret = false;
        }

        det (det_A_E, A, Method::Elimination());
        report << "Computed integer determinant (Elimination): " << det_A_E << endl;
        if (det_A_E != pi){
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: Elimination Computed determinant is incorrect" << endl;
            ret = false;
        }


        commentator().stop ("done");
        commentator().progress ();
            //commentator().progress (i, iterations);
    }

    commentator().stop (MSG_STATUS (ret), (const char *) 0, "testIntegerDeterminantGeneric");

    return ret;
}

/* Test 6: Rational determinant by generic methods
 *
 * Construct a random nonsingular diagonal sparse matrix and compute its
 * determinant over Z
 *
 * n - Dimension to which to make matrix
 * iterations - Number of iterations to run
 *
 * Returns true on success and false on failure
 */

bool testRationalDetGen (size_t n, int iterations)
{
    commentator().start ("Testing rational determinant, generic methods", "testRationalDeterminantGeneric", (unsigned int) iterations);

    bool ret = true;

    for (int i = 0; i < iterations; ++i) {
        commentator().startIteration ((unsigned int)i);
        Givaro::QField<Givaro::Rational> Q;
        SparseMatrix<Givaro::QField<Givaro::Rational> > A (Q, n, n);
        BlasMatrix <Givaro::QField<Givaro::Rational> > BB(Q, n, n);

        Givaro::QField<Givaro::Rational>::Element pi(1,1);
        Givaro::QField<Givaro::Rational>::Element det_A, det_B,det_A_H, det_B_H, det_A_B, det_B_B, det_A_E, det_B_E;

        for (unsigned int j = 0; j < n; ++j) {
            integer tmp_n;
            integer tmp_d;
            integer::nonzerorandom (tmp_n, 20*i + 1);
            integer::nonzerorandom (tmp_d, 20*i + 1);
            Givaro::QField<Givaro::Rational>::Element tmp;
            Q.init(tmp,tmp_n,tmp_d);
            A.setEntry(j,j,tmp);
            BB.setEntry(j,j,tmp);
            Q.mulin (pi, tmp);
        }

        if (i % 2) {
            Q.negin(A. refEntry(1,1));
            Q.negin(BB.refEntry(1,1));
            Q.negin(pi);
        }

        ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

        report << "True determinant: ";  Q.write(report,pi); report << endl;

        det (det_A, A);
        report << "Computed rational determinant (Default): "; Q.write(report, det_A); report << endl;
        if (!Q.areEqual(det_A ,pi)){
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: Default Computed determinant is incorrect" << endl;
            ret = false;
        }

        det (det_A_H, A, Method::Auto());
        report << "Computed rational determinant (Auto): "; Q.write(report, det_A_H); report << endl;
        if (!Q.areEqual(det_A_H ,pi)){
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: Auto Computed determinant is incorrect" << endl;
            ret = false;
        }

        det (det_A_B, A, Method::Blackbox());
        report << "Computed rational determinant (Blackbox): "; Q.write(report, det_A_B); report<< endl;
        if (!Q.areEqual(det_A_B , pi)){
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: Blackbox Computed determinant is incorrect" << endl;
            ret = false;
        }

        det (det_A_E, A, Method::Elimination());
        report << "Computed rational determinant (Elimination): "; Q.write(report, det_A_E); report << endl;
        if (!Q.areEqual(det_A_E ,pi)){
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: Elimination Computed determinant is incorrect" << endl;
            ret = false;
        }

        det (det_B, BB);
        report << "Computed rational determinant (Default): "; Q.write(report, det_A); report << endl;
        if (!Q.areEqual(det_B ,pi)){
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: (Dense) Default Computed determinant is incorrect" << endl;
            ret = false;
        }

        det (det_B_H, BB, Method::Auto());
        report << "Computed rational determinant (Auto): "; Q.write(report, det_A_H); report << endl;
        if (!Q.areEqual(det_B_H ,pi)){
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: (Dense) Auto Computed determinant is incorrect" << endl;
            ret = false;
        }
        det (det_B_B, BB, Method::Blackbox());
        report << "Computed rational determinant (Blackbox): "; Q.write(report, det_A_B); report<< endl;
        if (!Q.areEqual(det_B_B , pi)){
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: (Dense) Blackbox Computed determinant is incorrect" << endl;
            ret = false;
        }

        det (det_B_E, BB, Method::Elimination());
        report << "Computed rational determinant (DenseElimination): "; Q.write(report, det_A_E); report << endl;
        if (!Q.areEqual(det_B_E ,pi)){
            commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
                << "ERROR: (Dense) Elimination Computed determinant is incorrect" << endl;
            ret = false;
        }


        commentator().stop ("done");
        commentator().progress ();
            //commentator().progress (i, iterations);
    }

    commentator().stop (MSG_STATUS (ret), (const char *) 0, "testRationalDeterminantGeneric");

    return ret;
}


int main (int argc, char **argv)
{
    bool pass = true;

    static size_t n = 4;
    static integer q = 4093U;
    static int iterations = 1;

    static Argument args[] = {
        { 'n', "-n N", "Set dimension of test matrices to NxN", TYPE_INT,     &n },
        { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1]", TYPE_INTEGER, &q },
        { 'i', "-i I", "Perform each test for I iterations",    TYPE_INT,     &iterations },
        END_OF_ARGUMENTS
    };

    parseArguments (argc, argv, args);
    Givaro::Modular<int32_t> F (q);

    commentator().start("Determinant test suite", "det");

        // Make sure some more detailed messages get printed
    commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
    commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

    if (!testDiagonalDet1        (F, n, iterations)) pass = false;
    if (!testDiagonalDet2        (F, n, iterations)) pass = false;
    if (!testSingularDiagonalDet (F, n, iterations)) pass = false;
    if (!testIntegerDet          (n, iterations)) pass = false;
/*
  if (!testIntegerDetGen          (n, iterations)) pass = false;
  if (!testRationalDetGen          (n, iterations)) pass = false;
*/

    commentator().stop("determinant test suite");
    return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
