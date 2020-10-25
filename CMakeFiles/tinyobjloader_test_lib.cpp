
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"


int main()
{
  std::string inputfile = "null.obj";
  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t> shapes;
  std::vector<tinyobj::material_t> materials;
  
  std::string err;
  bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, inputfile.c_str());

  return 0;
}
