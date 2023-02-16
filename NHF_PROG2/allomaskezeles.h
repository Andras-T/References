#ifndef ALLOMASKEZELES_H
#define ALLOMASKEZELES_H

#include "allomas.h"
#include "vonatkezeles.h"

namespace Allomaskezeles
{
	/*!
	 *  uj allomas hozzaadasa
	 *
	 *      @param allomasok 
	 */
	void addAll(Allomas& allomasok)
	{
		std::string allomasnev;
		std::cout << "\n Allomasnev: ";
		std::cin >> allomasnev;

		try {
			allomasok.addAllomas(allomasnev);
		}
		catch (std::invalid_argument kivetel)
		{
			std::cout << "!! " << kivetel.what() << " !!" << std::endl;
		}
	}

	/*!
	 *  allomas kezeles (hozzaadas, torles)
	 *
	 *      @param allomasok 
	 */
	void allkezeles(Allomas& allomasok)
	{
		std::string choice;
		bool valasz = false;

		while (!valasz)
		{
			std::cout << std::endl << allomasok << std::endl;
			elvalaszto();
			std::cout << " 1 = Uj allomas hozzaadasa\n 2 = Utolso allomas eltavolitasa\n 3 = Kilepes\n ";
			std::cin >> choice;
			//Uj allomas
			if (choice == "1")
			{
				addAll(allomasok);
				valasz = true;
			}
			//Utolso allomas eltavolitasa
			else if (choice == "2")
			{
				try {
					allomasok.delAllomas();
				}
				catch (std::invalid_argument kivetel)
				{
					std::cout << " !!" << kivetel.what() << "!!\n";
				}
				valasz = true;
			}

			else if (choice == "3")
				valasz = true;
		}
	}
}

#endif // !ALLOMASKEZELES_H
