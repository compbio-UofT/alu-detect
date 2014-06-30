#include<string>

class Coord {
	public:
		Coord();
		Coord(int* coords, std::string chr);
		int& operator[](const int& i);
		int get(const int& i) const;
		bool operator<(const Coord& that) const;
		bool operator==(const Coord& that) const;
		std::string chr;
	private:
		int coords[5];
};
