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

		// Tworzy macierz z mo¿liwoœci¹ losowej inicjalizacji
		CMtx(std::size_t row_cnt, std::size_t col_cnt, bool rand_init = false);

		// Tworzy kwadratow¹ macierz diagonaln¹
		CMtx(std::size_t row_cnt, float diagonal);

		// Konstruktor kopiuj¹cy
		CMtx(const CMtx& rhs);

		// Jeœli potrzeba - nale¿y zadeklarowaæ i zaimplementowaæ inne konstruktory
		~CMtx();


		// =========================================================================
		// OPERATORY PRZYPISANIA:
		// =========================================================================

		const CMtx& operator=(const CMtx& rhs);

		// Zamiana macierzy na macierz diagonaln¹ 
		const CMtx& operator=(float diagonal);

		// Pobranie wymiarów macierzy (wykorzystywane przez wyj¹tki)
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

		// Mno¿enie macierzy przez sta³¹
		CMtx operator*(float multiplier) const;

		CMtx operator+(const CMtx& rhs) const;
		CMtx operator-(const CMtx& rhs) const;

		// Minus unarny - zmiana znaku wszystkich wspó³czynników macierzy
		CMtx operator-() const;

		// Transponowanie macierzy
		CMtx operator~() const;

		// Akceptuje tylko power >= 0:
		//    power = 0  - zwraca macierz jednostkow¹
		//    power = 1  - zwraca kopiê macierzy
		//    power > 1  - zwraca iloczyn macierzy 
		CMtx operator^(unsigned int power) const;

		// Porównywanie macierzy z dok³adnoœci¹ do sta³ej ALG_PRECISION
		bool operator==(const CMtx& rhs) const;

		// Tylko do celów testowych - wypisuje macierz wierszami na stdout
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