#ifndef NAVENGINE_INCLUDE_SATUTILS_CODE_GEN_HPP 
#define NAVENGINE_INCLUDE_SATUTILS_CODE_GEN_HPP 

#include <array>
#include <cstdint>

namespace satutils {

void CodeGenCA(std::array<bool,1023>& sequence, const uint8_t prn);
// void CodeGenL1C();
// void CodeGenL2CM();
// void CodeGenL2CL();
// void CodeGenL5I();
// void CodeGenL5Q();

} // namespace satutils
#endif
