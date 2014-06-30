#include<string>

#include "Coord.h"

Coord::Coord(){}
Coord::Coord(int* coords, std::string chr) {
	this->coords[1]=coords[1];
	this->coords[2]=coords[2];
	this->coords[3]=coords[3];
	this->chr=chr;
}
int& Coord::operator[](const int& i) {
	return coords[i];
}

int Coord::get(const int& i) const {
	return coords[i];
}

bool Coord::operator<(const Coord& that) const { //should return 0 if overlap
	return chr < that.chr || (chr == that.chr && coords[2] <= that.coords[1]);
}

bool Coord::operator==(const Coord& that) const {
	return !(*this < that) && !(that < *this);
}
