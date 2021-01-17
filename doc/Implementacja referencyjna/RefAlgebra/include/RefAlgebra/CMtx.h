#pragma once
#include <string>

#include "types.h"


namespace RefAlgebra
{
	class CMtx
	{
	public:
		static const float ALG_PRECISION;
		// =========================================================================
		// KONSTRUKTORY:
		// =========================================================================

		// Tworzy macierz z mo�liwo�ci� losowej inicjalizacji
		CMtx(std::size_t row_cnt, std::size_t col_cnt, bool rand_init = false);

		// Tworzy kwadratow� macierz diagonaln�
		CMtx(std::size_t row_cnt, float diagonal);

		// Konstruktor kopiuj�cy
		CMtx(const CMtx& rhs);

		// Je�li potrzeba - nale�y zadeklarowa� i zaimplementowa� inne konstruktory
		~CMtx();


		// =========================================================================
		// OPERATORY PRZYPISANIA:
		// =========================================================================

		const CMtx& operator=(const CMtx& rhs);

		// Zamiana macierzy na macierz diagonaln� 
		const CMtx& operator=(float diagonal);

		// Pobranie wymiar�w macierzy (wykorzystywane przez wyj�tki)
		Shape shape() const;

		// =========================================================================
		// INDEKSOWANIE MACIERZY
		// =========================================================================

		float* operator[](std::size_t row_ind);
		const float* operator[](std::size_t row_ind) const;

		// =========================================================================
		// OPERACJE ALGEBRAICZNE
		// =========================================================================

		CMtx operator*(const CMtx& rhs) const;

		// Mno�enie macierzy przez sta��
		CMtx operator*(float multiplier) const;

		CMtx operator+(const CMtx& rhs) const;
		CMtx operator-(const CMtx& rhs) const;

		// Minus unarny - zmiana znaku wszystkich wsp�czynnik�w macierzy
		CMtx operator-() const;

		// Transponowanie macierzy
		CMtx operator~() const;

		// Akceptuje tylko power >= 0:
		//    power = 0  - zwraca macierz jednostkow�
		//    power = 1  - zwraca kopi� macierzy
		//    power > 1  - zwraca iloczyn macierzy 
		CMtx operator^(unsigned int power) const;

		// Por�wnywanie macierzy z dok�adno�ci� do sta�ej ALG_PRECISION
		bool operator==(const CMtx& rhs) const;

		// Tylko do cel�w testowych - wypisuje macierz wierszami na stdout
		void display() const;

		// friend CMtx operator*(float multiplier, const CMtx& rhs);
	private:
		float**	    row_ptr;
		std::size_t row_cnt;
		std::size_t col_cnt;

		void assignCopy(CMtx other);

		bool isSameShape(const CMtx& other) const;
		std::string className() const;
	};

	CMtx operator*(float multiplier, const CMtx& rhs);
}