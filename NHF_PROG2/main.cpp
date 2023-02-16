#include <iostream>

#include "memtrace.h"


//osztalyok
#include "jarat.h"
#include "allomas.h"
#include "jegyek.h"
#include "vonat.h"
#include "kocsi.h"
#include "datum.h"
#include "vector.h"

//main-hez szukseges fuggvenyek
#include "allomaskezeles.h"
#include "vonatkezeles.h"
#include "jaratkezeles.h"


//megjeleniti a menut
void showMenu(std::string const menu[], const size_t menumeret)
{
	for (size_t i = 0; i < menumeret; i++)
	{
		std::cout << menu[i]  << std::endl;
	}
	std::cout << " ";
}

int main()
{
	/*-------------------------ADATOK BETOLTESE-------------------------*/
	Allomas* all;
	//betolti fajlbol az allomasokat
	all = new Allomas;	
	vectorClass <Vonat> vonatok;
	vectorClass <Jarat> jaratok;


	//betolti fajlbol a vonatokat
	Vonatkezeles::vonat_betolt(vonatok);
	Datum TESZT_IDO(11, 15);
	Jarat TESZT("Tenkes", TESZT_IDO, TESZT_IDO,"Pecs", "Vac");
	jaratok.push(TESZT);
	Jaratkezeles::megjelenit(jaratok, vonatok);
	jaratok[0].foglal(1, 1, vonatok);
	jaratok[0].foglal(1, 2, vonatok);
	jaratok[0].foglal(1, 3, vonatok);
	jaratok[0].kiir(&vonatok);
	//Menu merete
	size_t const menumeret = 7;
	//Menurendszer
	std::string const menu[menumeret] =
	{
		"\n [1]\tAllomasok kezelese\n",
		" [2]\tVonatok kezelese\n",
		" [3]\tKocsik kezelese\n",
		" [4]\tJaratok\n",
		" [5]\tFoglalas\n",
		" [6]\tJegyek\n",
		" [9]\tKilepes\n " };	

	std::string menupont, almenupont;

	/*-------------------------MENU-------------------------*/
	do
	{
	showMenu(menu, menumeret);
	std::cin >> menupont;
	

	/*-------------------------Allomasok kezelese-------------------------*/
	if (menupont == "1")
	{
		elvalaszto();	
		Allomaskezeles::allkezeles(*all);
		elvalaszto();
		std::cout << *all;
		elvalaszto();
	}
	/*-------------------------Vonatok kezelese-------------------------*/
	else if (menupont == "2")
	{
		Vonatkezeles::megjelenit(vonatok);
		std::cout << "\n [1]Vonat hozzaadasa \n [2] Utolso vonat eltavolitasa\n ";
		std::cin >> almenupont;
		if (almenupont == "1")
		{
			Vonatkezeles::addvonat(vonatok);
		}
		else if (almenupont == "2")
		{
			Vonatkezeles::delvonat(vonatok);
		}
		Vonatkezeles::megjelenit(vonatok);
	}
	/*-------------------------Kocsik kezelese-------------------------*/
	else if (menupont == "3")
	{
		Vonatkezeles::megjelenit(vonatok);
		std::cout << "\n [1] Uj kocsi vonathoz \n [2] kocsi eltavolitasa\n ";
		std::cin >> almenupont;
		if (almenupont == "1")
		{
			Vonatkezeles::kocsi_add(vonatok);
		}
		else if (almenupont == "2")
		{
			Vonatkezeles::kocsi_del(vonatok);
		}
		Vonatkezeles::megjelenit(vonatok);
	}
	/*-------------------------Jaratok kezelese-------------------------*/
	else if (menupont == "4")
	{
		elvalaszto();
		Jaratkezeles::megjelenit(jaratok, vonatok);
		std::cout << "\n [1] Uj jarat felvetele \n [2] Utolso jarat torlese \n ";
		std::cin >> almenupont;
		if (almenupont == "1")
		{
			Jaratkezeles::addjarat(vonatok, jaratok, *all);
		}
		else if (almenupont == "2")
		{
			Jaratkezeles::deljarat(jaratok);
		}
		Jaratkezeles::megjelenit(jaratok, vonatok);
		elvalaszto();
	}
	/*-------------------------Foglalasok kezelese-------------------------*/
	else if (menupont == "5")
	{
		elvalaszto();
		Jaratkezeles::megjelenit(jaratok, vonatok);
		if (jaratok.size() != 0)
		{
			std::cout << "\n Valassza ki a jarat szamat: ";
			std::cin >> almenupont;
			size_t sorszam = atoi(almenupont.c_str());
			try
			{
				Jaratkezeles::foglal(jaratok[sorszam - 1], vonatok);
			}
			catch (const std::invalid_argument)
			{
				std::cout << "\n Nincs ilyen jarat!\n";
			}	
		}	
		elvalaszto();
	}
	/*-------------------------Jegyek-------------------------*/
	else if (menupont == "6")
	{
		elvalaszto();
		if (jaratok.size() != 0)
		{
			Jaratkezeles::megjelenit(jaratok, vonatok);
			std::cout << "\n Valassza ki a jarat szamat: ";
			std::cin >> almenupont;
			size_t sorszam = atoi(almenupont.c_str());
			try
			{
				jaratok[sorszam - 1].jegykiir();
			}
			catch(std::invalid_argument)
			{
				std::cout << " !!Nincs ilyen!!\n";
			}
		}
		else
		{
			std::cout << "\n Meg nincsenek jegyek\n";
		}
		elvalaszto();
	}
	} while (menupont != "9");
	/*--------------------------------------------------*/
	delete all;
}