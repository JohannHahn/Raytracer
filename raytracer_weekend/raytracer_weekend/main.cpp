#include <iostream>

void printConsole(char* s)
{
	std::cout << s << std::endl;
}

int main()
{
	printConsole("test");
	std::cin.get();
}

