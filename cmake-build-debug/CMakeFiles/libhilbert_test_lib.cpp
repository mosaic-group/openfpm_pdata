
extern "C"
{
#include "hilbertKey.h"
}
  int main()
{
  //An integer to handle errors
  int err;

  //Array to handle output
  uint64_t nextCoord[2];

  //Get the coordinates of the next cell
  getIntCoordFromHKey(nextCoord, 4, 2, 0, &err);

  return 0;
}
