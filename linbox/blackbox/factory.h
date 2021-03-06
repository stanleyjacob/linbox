/* linbox/blackbox/factory.h
 * Copyright (C) 2002  LinBox
 *
 * Written by Bradford Hovinen <bghovine@math.uwaterloo.ca>
 *
 * ------------------------------------
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

#ifndef __LINBOX_blackbox_factory_H
#define __LINBOX_blackbox_factory_H

#include "linbox/util/error.h"
#include "linbox/vector/vector-traits.h"

namespace LinBox
{

	/** @brief A tool for computations with integer and rational matrices.
	 *
	 * The blackbox factory provides a facility for performing integer or rational
	 * computations by reducing modulo one or more primes and recovering the
	 * solution with Chinese Remaindering, lifting, or rational reconstruction. It
	 * is an interface that provides one method which, given a field, produces a
	 * black box representing a particular matrix over that field. The factory
	 * object may be passed to various procedures, such as rank, det, and solve,
	 * which will perform the required modular reductions to find integer or
	 * rational solutions.
	 *
	 * In the typical case, the user provides an object whose class inherits from
	 * BlackboxFactory and implements the method makeBlackbox. The object represents
	 * the original integer or rational version of the black box, whose data might
	 * require some modification (e.g. modular reduction) to produce a true black
	 * box. Alternatively, the resulting black box might merely be a
	 * reinterpretation of the data in the original object, as is the case where
	 * matrix entries are all nonnegative and smaller than the modulus.
	 */

	template <class Field, class Blackbox>
	class BlackboxFactory {
	public:

		/// Virtual destructor
		virtual ~BlackboxFactory () {}

		/** Given a field and vector type, construct a black box for the matrix
		 * over that field and using that vector type. This should be
		 * implemented by the user
		 * @param F Field over which to construct the black box
		 */
		virtual Blackbox *makeBlackbox (const Field &F) = 0;

		/** Compute and return the max-norm of the matrix.
		 *
		 * @param res Place to store result
		 */
		virtual integer &maxNorm (integer &res) = 0;

		/** Give the row dimension of the matrix
		*/
		virtual size_t rowdim () = 0;

		/** Give the column dimension of the matrix
		*/
		virtual size_t coldim () = 0;

	}; // BlackboxFactory

} // namespace LinBox

#endif // __LINBOX_blackbox_factory_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
