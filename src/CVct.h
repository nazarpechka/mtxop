#pragma once

typedef double FPTYPE;

#include "CMtx.h"
#include "CAlgError.h"

namespace MyAlgebra
{
	class CMtx;

	class CVct
	{
	private:
		FPTYPE  *vector;
		int     size;

	public:
		CVct(int size);

		FPTYPE & operator[](int ind) const;

		const CVct & operator=( const CVct & rhs );
		const CVct & operator=( FPTYPE val );

		CVct operator+( const CVct & rhs );
		CVct operator-( const CVct & rhs );

		CVct operator*( const CMtx & rhs );

		// Transpozycja - zamiana wektora wierszowego na kolumnowy i odwrotnie
		CVct operator~();
	};
}

