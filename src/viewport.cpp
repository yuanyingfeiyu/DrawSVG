#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox( float centerX, float centerY, float vspan ) {

  // Task 5 (part 2): 
  // Set svg coordinate to normalized device coordinate transformation. Your input
  // arguments are defined as normalized SVG canvas coordinates.
  this->centerX = centerX;
  this->centerY = centerY;
  this->vspan = vspan; 
  double tx[] = { 1.0, 0.0, -centerX + vspan,
				0.0, 1.0, -centerY + vspan,
				0.0, 0.0, 1.0 };

  double sc[] = { 1.0 /(2 * vspan), 0.0, 0.0,
				0.0, 1.0 / (2 * vspan), 0.0,
				0.0, 0.0, 1.0};

  // first translation and then scale and compose those two matrix and set it to svg_2_norm
  Matrix3x3 translation = Matrix3x3(tx);
  Matrix3x3 scale = Matrix3x3(sc);
  set_svg_2_norm(scale * translation);
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->centerX -= dx;
  this->centerY -= dy;
  this->vspan *= scale;
  set_viewbox( centerX, centerY, vspan );
}

} // namespace CMU462
