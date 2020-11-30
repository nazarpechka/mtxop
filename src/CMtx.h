#pragma once
#include "CVct.h"

namespace MyAlgebra
{
	class CVct;

	class CMtx
	{	
	private:
		FPTYPE  **row_ptr;
		int     row_cnt;
		int     col_cnt;

	public:
		static const FPTYPE ALG_PRECISION;

		// =========================================================================
		// KONSTRUKTORY:
		// =========================================================================

		// Tworzy macierz z mo¿liwoœci¹ losowej inicjalizacji
		CMtx( int row_cnt, int col_cnt, bool rand_init = false );

		// Tworzy kwadratow¹ macierz diagonaln¹
		CMtx( int row_cnt, FPTYPE diagonal );

		CMtx(const CMtx & rhs);

		// Jeœli potrzeba - nale¿y zadeklarowaæ i zaimplementowaæ inne konstruktory
		~CMtx();


		// =========================================================================
		// OPERATORY PRZYPISANIA:
		// =========================================================================

		const CMtx & operator=( const CMtx & rhs );

		// Zamiana macierzy na macierz diagonaln¹ 
		const CMtx & CMtx::operator=( const FPTYPE diagonal );

		// Operator pzzenosz¹cy
		const CMtx & operator=(CMtx && rhs);


		// =========================================================================
		// INDEKSOWANIE MACIERZY
		// =========================================================================

		FPTYPE * operator[](int row_ind);

		// =========================================================================
		// OPERACJE ALGEBRAICZNE
		// =========================================================================

		// Mnozenie macierzy przez wektor, rhs musi byæ wektorem kolumnowym
		CVct operator*( const CVct & rhs ) const;

		CMtx operator*( const CMtx & rhs ) const;

		// Mno¿enie macierzy przez sta³¹
		CMtx operator*( FPTYPE multiplier ) const;

		CMtx operator+( const CMtx & rhs ) const;
		CMtx operator-( const CMtx & rhs ) const;

		// Minus unarny - zmiana znaku wszystkich wspó³czynników macierzy
		CMtx operator-() const;

		// Transponowanie macierzy
		CMtx operator~() const;

		// Akceptuje tylko power >= -1:
		//    power = -1 - zwraca macierz odwrócon¹
		//    power = 0  - zwraca macierz jednostkow¹
		//    power = 1  - zwraca kopiê macierzy
		//    power > 1  - zwraca iloczyn macierzy 
		CMtx operator^( int power ) const;

		// Porównywanie macierzy z dok³adnoœci¹ do sta³ej ALG_PRECISION
		bool operator==(const CMtx && rhs) const;

		// Tylko do celów testowych - wypisuje macierz wierszami na stdout
		void display() const;

		// friend CMtx operator*( FPTYPE multiplier, const CMtx &rhs );
	};

	CMtx operator*(FPTYPE multiplier, const CMtx &rhs);
}


