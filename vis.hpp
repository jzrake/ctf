

#ifndef __VisualizationOpenGl_HEADER__
#define __VisualizationOpenGl_HEADER__

#include "hydro.hpp"

class VisualizationOpenGl
{
private:
  void load_texture();

public:
  VisualizationOpenGl();
  ~VisualizationOpenGl();
  void DrawScene();
} ;


#endif // __VisualizationOpenGl_HEADER__
