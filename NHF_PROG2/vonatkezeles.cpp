#include <iostream>
#include <string>
#include <fstream>
#include "vonatkezeles.h"
#include "memtrace.h"



void elvalaszto()
{
std::cout << "---------------------------------------------------------\n";
}

void Vonatkezeles::addvonat(vectorClass <Vonat>& vonatok, const char* filename)
{
	std::string vonatneve;
	std::cout << " Vonat neve: ";
	std::cin >> vonatneve;
	Vonat uj(vonatneve);
	for (size_t i = 0; i < vonatok.size(); i++)
	{
		if (vonatok[i]== vonatneve)
		{
			std::cout << "\n Mar letezik vonat ezzel a nevvel\n";
			return;
		}
	}
	std::cout << "\n addvonat \n";
	vonatok.push(uj);
	std::ofstream MyFile(filename);
	MyFile << vonatok;
	MyFile.close();
	std::cout << " Vonat hozzaadva\n ";
}

void Vonatkezeles::delvonat(vectorClass <Vonat>& vonatok, const char* filename)
	{
		try
		{
			vonatok.pop();
			std::ofstream MyFile(filename);
			MyFile << vonatok;
			MyFile.close();
		}
		catch (std::invalid_argument kivetel)
		{
			std::cout << " " << kivetel.what() << std::endl;
		}
	}

void Vonatkezeles::kocsi_add(vectorClass <Vonat>& vonatok, const char* filename)
	{
		std::string nev;
		std::cout << " Vonat neve: ";
		std::cin >> nev;
		//megkeresi  a megfelelo vonatot
		for (size_t i = 0; i < vonatok.size(); i++)	
		{
			Vonat* temp = &vonatok[i];
			if (temp->operator==(nev))
			{
				int kocsiszam, maxhely, osztaly;
				std::cout << "\n Kocsiszam: ";
				std::cin >> kocsiszam;
				std::cout << "\n Kocsiosztaly [1/2]: ";
				std::cin >> osztaly;
				//Nem megfelelo osztaly
				if (osztaly > 2 || osztaly < 0)
				{
					std::cout << " !!Nincs ilyen osztaly!!\n";
					return;
				}
				std::cout << "\n Max helyek szama: ";
				std::cin >> maxhely;
				//Kocsi hozzaadasa
				temp->addkocsi(maxhely, osztaly, kocsiszam);
				temp->kiir();
				std::ofstream MyFile(filename);
				//FILE-ba iras
				MyFile << vonatok;
				MyFile.close();
				return;
			}
		}
		std::cout << " !!Nincs ilyen vonat!!\n";
	}

void Vonatkezeles::kocsi_del(vectorClass <Vonat>& vonatok, const char* filename)
	{
		std::string nev;
		std::cout << "Vonat neve: ";
		std::cin >> nev;
		for (size_t i = 0; i < vonatok.size(); i++)	//Vonat megkeresese
		{
			Vonat* temp = &vonatok[i];
			if (temp->operator==(nev))
			{
				try
				{
					temp->delkocsi();
					std::ofstream MyFile(filename);
					MyFile << vonatok;
					MyFile.close();
				}
				catch (std::invalid_argument kivetel)	//nincs tobb kocsi
				{
					std::cout << kivetel.what();
				}

				return;
			}
		}
		std::cout << "Nincs ilyen vonat\n";
	}

void Vonatkezeles::szet_valogat(std::string buf, vectorClass <Vonat>& vonatok)
	{
		//elvalaszto
		std::string delimiter = ": ";
		//szo vegenek a pozicioja
		size_t pos = 0;
		size_t index = 0;

		Vonat aktualis;
		std::string vonatnev, kocsiszam, osztaly, maxhely;

		while ((pos = buf.find(delimiter)) != std::string::npos)
		{
			index++;

			switch (index)
			{
			case 1:	//vonatnev
				vonatnev = buf.substr(0, pos);
				aktualis.atnevez(vonatnev);
				vonatok.push(aktualis);
				delimiter = ",";
				break;
			case 2:	//kocsiszam
				kocsiszam = buf.substr(0, pos);
				break;
			case 3:	//kocsiosztaly
				osztaly = buf.substr(0, pos);
				delimiter = ";";
				break;
			case 4:	//maxhely
				maxhely = buf.substr(0, pos);
				vonatok[vonatok.size() - 1].addkocsi(atoi(maxhely.c_str()), atoi(osztaly.c_str()), atoi(kocsiszam.c_str()));
				index = 1;
				delimiter = ",";
				break;
			default:
				break;
			}
			//sztring roviditese
			buf.erase(0, pos + delimiter.length());
		}
	}

void Vonatkezeles::vonat_betolt(vectorClass <Vonat>& vonatok, const char* filename)
	{	
		// FILE megnyitása olvasashoz
		std::ifstream MyFile(filename);	
		//ideiglenes tarolo
		std::string buf;	
		//vonatok soronkenti beolvasasa (1 vonat = 1 sor)
		while (std::getline(MyFile, buf))	
		{
			Vonatkezeles::szet_valogat(buf, vonatok);
		}
		//FILE bezarasa
		MyFile.close();	
		//File iras modban
		std::ofstream MyReadFile(filename);
		//vonatok visszairasa fajlba
		MyReadFile << vonatok;	
		//FILE bezarasa
		MyReadFile.close();	
	}

void Vonatkezeles::megjelenit(vectorClass <Vonat>& vonatok)
	{
		elvalaszto();
		std::cout << "\n Vonatok:\n";
		for (size_t i = 0; i < vonatok.size(); i++)
		{
			std::cout << "\n [" << i + 1 << ".] ";
			vonatok[i].kiir();
		}
		elvalaszto();
	}
