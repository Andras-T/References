#include <iostream>
#include "jarat.h"
#include "memtrace.h"


Jarat::Jarat(const std::string vonat, const Datum ind, const Datum erk, const std::string i_all, const std::string e_all)
	: vonat(vonat), i_allomas(i_all), e_allomas(e_all), indulas(ind), erkezes(erk)
{}

Jarat& Jarat::operator=(const Jarat& other)
{
	vonat = other.vonat;
	i_allomas = other.i_allomas;
	e_allomas = other.e_allomas;
	indulas = other.indulas;
	erkezes = other.erkezes;
	jegyek = other.jegyek;
	diakjegyek = other.diakjegyek;
	felnottjegyek = other.felnottjegyek;
	biciklijegyek = other.biciklijegyek;
	return *this;
}

void Jarat::foglal(const int k_osztaly, const int jegytipus, vectorClass <Vonat>& vonatok)
{
	Kocsi* aktualis;
	int hely= 0;
	int kocsiszam = 0;
	try
	{	//vonat megkeresese
	for (size_t i = 0; i < vonatok.size(); i++)
	{
		if (vonatok[i] == vonat)
		{
			aktualis = vonatok[i].getkocsi(k_osztaly);
			hely = aktualis->getszabad();
			aktualis->subszabad();
			kocsiszam = aktualis->getkocsiszam();
			break;
		}
	}	
		
	}
	catch (std::invalid_argument kivetel)
	{
		throw kivetel;
		return;
	}

	if (jegytipus == 0)	//alapjegy
	{
		Jegyek uj(hely, k_osztaly, kocsiszam);
		jegyek.push(uj);
		std::cout << "\n" << uj << " Jegyar: " << uj.getar() << std::endl;
	}
	else if (jegytipus == 1)	//felnottjegy
	{
		FelnottJegy uj(hely, k_osztaly, kocsiszam);
		felnottjegyek.push(uj);
		std::cout << "\n" << uj;
	}
	else if (jegytipus == 2)	//diakjegy
	{
		DiakJegy uj(hely, k_osztaly, kocsiszam);
		diakjegyek.push(uj);
		std::cout << "\n" << uj;
	}
	else if (jegytipus == 3)	//biciklijegy
	{
		Bicikli uj(hely, k_osztaly, kocsiszam);
		biciklijegyek.push(uj);
		std::cout << "\n" << uj;
	}
	else
		throw std::invalid_argument("Nincs ilyen jegytipus");

}

void Jarat::kiir(vectorClass <Vonat>* vonatok) const
{
	size_t index=0;
	for (size_t i = 0; i < vonatok->size(); i++)
	{
		if ((*vonatok)[i] == vonat)
		{
			index = i;
			break;
		}
	}
	std::cout << " " << vonat << ":\n Indulo allomas: " << i_allomas << " [" << indulas <<
		"]\n Erkezo allomas: " << e_allomas << " [" << erkezes << "]\n Elerheto osztalyok: -";
	for (size_t i = 0; i < (*vonatok)[index].getsize(); i++)
	{
		if (((*vonatok)[index])[i].getosztaly() == 1)
		{
			std::cout << "[1] ";
			break;
		}
	}
	for (size_t i = 0; i < (*vonatok)[index].getsize(); i++)
	{
		if (((*vonatok)[index])[i].getosztaly() == 2)
		{
			std::cout << "[2]";
			break;
		}
	}
	std::cout << std::endl;
}

void Jarat::jegykiir() const
{
	if (jegyek.size() > 0)
	{
		std::cout << vonat << ":\n Alapjegyek:\n";
		for (size_t i = 0; i < jegyek.size(); i++)
		{
			std::cout << "\n [" << i + 1 << "] " << jegyek[i];
		}
	}

	if (felnottjegyek.size() > 0)
	{
		std::cout << "\n Felnottjegyek:\n";
		for (size_t i = 0; i < felnottjegyek.size(); i++)
		{
			std::cout << "\n [" << i + 1 << "] " << felnottjegyek[i];
		}
	}
	if (diakjegyek.size() > 0)
	{
		std::cout << "\n Diakjegyek:\n";
		for (size_t i = 0; i < diakjegyek.size(); i++)
		{
			std::cout << "\n [" << i + 1 << "] " << diakjegyek[i];
		}
	}

	if (biciklijegyek.size())
	{
		std::cout << "\n Biciklijegyek:\n";
		for (size_t i = 0; i < biciklijegyek.size(); i++)
		{
			std::cout << "\n [" << i + 1 << "] " << biciklijegyek[i];
		}
	}
}