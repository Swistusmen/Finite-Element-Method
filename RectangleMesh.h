#pragma once
#include <iostream>
#include "Element.h"

//domyslnie maksymalna liczba krawedzi dla kazdego wierzcholka- 4, komorki w srodku to prostokaty
//nie przechowuje elementow wprost w pamieci-bede je pozyskiwal przez obliczenia, zeby nie marnowac pamieci
//z racji tego ze element bedzie potrzebowal danych- pozniej pewnie zaimplemntuje kontener ze strukturami przechowujacymi dane

namespace data {
	class RectangleMesh {
	public:
		RectangleMesh(float H, float W, int nH, int nW, int* tab, int sizeofTab); //tab- tablica wartosci inicjalizujacych Elementy, element*2
		RectangleMesh(float H, float W, int nH, int nW, Element* elements);
		int maxIndexOfElement(); // jest tez liczba nodow
		bool isSuchAnElement(int ID); //narazie zwraca inta- czy istnieje, ale potem przekrztalce na zwracanie struktory posiadajacej info
		Element* getNodesID(int ID);
		Element& getElement(int ID);

	private:
		float h, w; //fizyczne parametry
		Element* elements;
		int nH; //wysokosc w nodach
		int nW;
	};
}