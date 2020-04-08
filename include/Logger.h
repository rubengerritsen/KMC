
#include <sstream>
#include <mutex>
#include <iostream>

class Logger: public std::stringstream
{
public:
    Logger() = default;

    ~Logger()
    {
        std::cout << this->str();
    }
};
