#include "meshio.hpp"
#include <tessellation.h>
#include <QImage>

using namespace tessellation;

int main()
{
    tessellation::Tessellation::process("agate.obj", "height", 2);

    return 0;
}

