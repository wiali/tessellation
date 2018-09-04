#include "meshio.hpp"
#include <tessellation.h>
#include <QImage>

using namespace tessellation;

int main()
{
    tessellation::Tessellation::process("agate", "height", 2);

    return 0;
}

