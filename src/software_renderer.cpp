#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CMU462 {


// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {
  sample_buffer = vector<unsigned char>(w * h * 4, 255);
  // set top level transformation
  transformation = svg_2_screen;

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y--;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();
}

void SoftwareRendererImp::fill_sample( int sx, int sy, const Color& color ) {
  if(sx < 0 || sx >= w) return;
  if(sy < 0 || sy >= h) return;

  sample_buffer[4 * (sx + sy * w)    ] = (uint8_t) (color.r * 255);
  sample_buffer[4 * (sx + sy * w) + 1] = (uint8_t) (color.g * 255);
  sample_buffer[4 * (sx + sy * w) + 2] = (uint8_t) (color.b * 255);
  sample_buffer[4 * (sx + sy * w) + 3] = (uint8_t) (color.a * 255);
};

void SoftwareRendererImp::fill_pixel( int x, int y, const Color& color ) {
  render_target[4 * (x + y * target_w)    ] = (uint8_t) (color.r * 255);
  render_target[4 * (x + y * target_w) + 1] = (uint8_t) (color.g * 255);
  render_target[4 * (x + y * target_w) + 2] = (uint8_t) (color.b * 255);
  render_target[4 * (x + y * target_w) + 3] = (uint8_t) (color.a * 255);
};


void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;
  this->w = sample_rate * target_w;
  this->h = sample_rate * target_h;
  sample_buffer = vector<unsigned char>(w * h * 4, 255);
}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;
  set_sample_rate(sample_rate);
}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

  // Task 5 (part 1):
  // Modify this to implement the transformation stack
  Matrix3x3 oldTransform = transformation;
  transformation = transformation * element->transform;

  switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
  }

  transformation = oldTransform;

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  // rasterize_point( p.x, p.y, point.style.fillColor );
  paint_point( p.x, p.y, point.style.fillColor);

}

void SoftwareRendererImp::paint_point( float x, float y, Color fillColor) {
  int sx = (int) floor(x);
  int sy = (int) floor(y);

  // fill sample - NOT doing alpha blending!
  // printf("%s \n", sample_buffer);
  // fill_sample(sx, sy, color);
  for (int i = 0; i < sample_rate; i++) {
    for (int j = 0; j < sample_rate; j++) {
      rasterize_point(sx * sample_rate + i, sy * sample_rate + j, fillColor);
    }
  }
}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int) floor(x);
  int sy = (int) floor(y);

  // check bounds
  if ( sx < 0 || sx >= w ) return;
  if ( sy < 0 || sy >= h ) return;

  // fill sample - NOT doing alpha blending!
  // printf("%s \n", sample_buffer);
  fill_sample(sx, sy, color);

}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {
  x0 = x0 * sample_rate;
  y0 = y0 * sample_rate;
  x1 = x1 * sample_rate;
  y1 = y1 * sample_rate;

  // Task 2: 
  // Implement line rasterization
  // implemented with Bresenham's algorithm 
  bool steep = abs(y1 - y0) > abs(x1 - x0);
  if(steep) {
    swap(x0, y0);
    swap(x1, y1);
  }

  if(x0 > x1) {
    swap(x0, x1);
    swap(y0, y1);
  }

  int dx = x1 - x0,
      dy = y1 - y0,
      y = y0,
      eps = 0;

  bool negative = dy * dx < 0;

  if (!negative) {
    for(int x = x0; x <= x1; x++ ) {
      if (steep) {
        rasterize_point(y, x, color);
      } else {
        rasterize_point(x, y, color);
      }
      
      eps += dy;
      if ((eps << 1) >= dx) {
        y++;
        eps -= dx;
      }
    }  
  } else {
    for(int x = x0; x <= x1; x++) {
      if (steep) {
        rasterize_point(y, x, color);
      } else {
        rasterize_point(x, y, color);
      }
      
      eps += dy;
      if ((eps << 1) <= - dx) {
        y--;
        eps += dx;
      }
    }
  }

  // task 2:
  // Implement line rasterization
  // to be continued: implemented with Xiaolin Wu's line algorithm

}

int sign(float value) {
  return value == 0 ? 0 : value > 0 ? 1 : -1;
}

bool in_triangle(Vector2D p1, Vector2D p2, Vector2D p3, Vector2D point) {
  Vector2D a = p2 - p1;
  Vector2D b = p3 - p2;
  Vector2D c = p1 - p3;

  Vector2D u1 = point - p1;
  Vector2D u2 = point - p2;
  Vector2D u3 = point - p3;

  int s1 = sign(cross(a, u1));
  double p = dot(a, u1) / a.norm2();
  if(s1 == 0 && p >= 0 && p <= 1) return true;

  int s2 = sign(cross(b, u2));
  p = dot(b, u2) / b.norm2();
  if(s2 == 0 && p >= 0 && p <= 1) return true;

  int s3 = sign(cross(c, u3));
  p = dot(c, u3) / c.norm2();
  if (s3 == 0 && p >= 0 && p <= 1) return true;

  return s1 == s2 && s2 == s3;
}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  x0 = x0 * sample_rate;
  y0 = y0 * sample_rate;
  x1 = x1 * sample_rate;
  y1 = y1 * sample_rate;
  x2 = x2 * sample_rate;
  y2 = y2 * sample_rate;

  // Task 3: 
  // Implement triangle rasterization
  float hmax = ceil(fmax(fmax(x0, x1), x2)), 
        hmin = floor(fmin(fmin(x0, x1), x2)), 
        vmax = ceil(fmax(fmax(y0, y1), y2)), 
        vmin = floor(fmin(fmin(y0, y1), y2));

  for(int i = hmin; i + 0.5 < hmax; i++) {
    for(int j = vmin; j + 0.5 < vmax; j++) {
      if(in_triangle(Vector2D(x0, y0), Vector2D(x1, y1), Vector2D(x2, y2), Vector2D(i + 0.5, j + 0.5))) {
        rasterize_point(i, j, color);
      }
    }
  }

}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 6: 
  // Implement image rasterization
  float dx = x1 - x0;
  float dy = y1 - y0;

  float period = 1.0f / sample_rate;
  float offset = period * 0.5f;
  float imgW = x1 - x0;
  float imgH = y1 - y0;

  for (float y = floor(y0); y <= floor(y1); y++) {
      for (float x = floor(x0); x <= floor(x1); x++) {
          for (int i = 0; i < sample_rate; i++) {
              for (int j = 0; j < sample_rate; j++) {
                  // sample locations (pixel coordinates)
                  float xs = x + j * period + offset;
                  float ys = y + i * period + offset;
                  
                  float u = (xs - x0) / imgW;
                  float v = (ys - y0) / imgH;
                  Color c = sampler->sample_trilinear(tex, u, v, imgW, imgH);

                  int sx = (int)floor(xs * sample_rate);
                  int sy = (int)floor(ys * sample_rate);
                  fill_sample(sx, sy, c);
              }
          }
      }
  }
}


void apply_filter(unsigned char* supersample_target, unsigned char* target_buffer, size_t sample_rate, size_t target_w, size_t target_h, int x, int y) {
  
  for(int bit = 0; bit < 4; bit++) {
    int sum = 0;
    for(int i = 0; i < sample_rate; i++) {
      for(int j = 0; j < sample_rate; j++) {
        sum += supersample_target[4*((x+i) + (y+j) * target_w * sample_rate) + bit];
      }
    }
    target_buffer[4*(x + y * target_w) + bit] = (uint8_t) sum / (sample_rate * sample_rate);
  }
}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 4: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 4".
  

  // first: apply unit-area box filter and set the render_target
  // printf("%d, %d", w, h);
  int filter_count = sample_rate * sample_rate;
  for(int x = 0; x < target_w; x++) {
    for(int y = 0; y < target_h; y++) {
        float r = 0, g = 0, b = 0, a = 0;
        for(int i = 0; i < sample_rate; i++) {
          for(int j = 0; j < sample_rate; j++) {
            int index = 4* (x * sample_rate + i + (y * sample_rate + j) * w);
            r += sample_buffer[index + 0] / 255.0f / filter_count;
            g += sample_buffer[index + 1] / 255.0f / filter_count;
            b += sample_buffer[index + 2] / 255.0f / filter_count;
            a += sample_buffer[index + 3] / 255.0f / filter_count;
          }
        }

        fill_pixel(x, y, Color(r, g, b, a));
    }
  }
  return;

}


} // namespace CMU462
