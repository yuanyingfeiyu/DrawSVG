#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Task 7: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

    for(int x = 0; x < level.width; x++) {
      for(int y = 0; y < level.height; y++) {

        int redSum = 0, greenSum = 0, blueSum = 0, alphaSum = 0;
        MipLevel& lastLevel = tex.mipmap[startLevel + i - 1];

                // four neighbors
        for(int j = 0; j < 2; j++) {
          for(int k = 0; k < 2; k++) {
            // RGBA
            int index = 4 * (x+j + (y+k) * lastLevel.width);
            redSum += lastLevel.texels[index];
            greenSum += lastLevel.texels[index+1];
            blueSum += lastLevel.texels[index+2];
            alphaSum += lastLevel.texels[index+3];
          }
        }

        level.texels[x + y * level.width] = redSum / 4;
        level.texels[x + y * level.width + 1] = greenSum / 4;
        level.texels[x + y * level.width + 2] = blueSum / 4;
        level.texels[x + y * level.width + 3] = alphaSum / 4;
      }
    }

  }

  // fill all 0 sub levels with interchanging colors (JUST AS A PLACEHOLDER)
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  if (level >= tex.mipmap.size())
    return Color(1,0,1,1);
  Color c = Color();  MipLevel& mip = tex.mipmap[level];
  vector<unsigned char>& texels = mip.texels;
  int x = (int) max(0.0f, floor(u * mip.width - 0.5f));
  int y = (int) max(0.0f, floor(v * mip.height - 0.5f));
  int index = 4 * (x + y * mip.width);
  uint8_to_float(&c.r, &texels[index]);
  return c;


}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task 6: Implement bilinear filtering
  if(level >= tex.mipmap.size()) 
    return Color(1,0,1,1);
  MipLevel& mip = tex.mipmap[level];
  vector<unsigned char>& texels = mip.texels;

  float x = u * mip.width, y = v * mip.height;

  int f00_x = (int) max(0.0f, floor(x - 0.5f));
  int f00_y = (int) max(0.0f, floor(y - 0.5f));

  int f01_x = f00_x;
  int f01_y = min((int) mip.height-1, f00_y+1);

  int f10_x = min((int) mip.width-1, f00_x+1);
  int f10_y = f00_y;

  int f11_x = min((int) mip.width-1, f00_x+1);
  int f11_y = min((int) mip.height-1, f00_y+1);

  float s = clamp(x - f00_x - 0.5f, 0.0f, 1.0f);
  float t = clamp(y - f00_y - 0.5f, 0.0f, 1.0f);

  Color c_f00 = Color(), c_f01 = Color(), c_f10 = Color(), c_f11 = Color();

  uint8_to_float(&c_f00.r, &texels[4 * (f00_x + f00_y * mip.width)]);
  uint8_to_float(&c_f01.r, &texels[4 * (f01_x + f01_y * mip.width)]);
  uint8_to_float(&c_f10.r, &texels[4 * (f10_x + f10_y * mip.width)]);
  uint8_to_float(&c_f11.r, &texels[4 * (f11_x + f11_y * mip.width)]);

  return (1.0f - t) * ((1.0f - s) * c_f00 + s * c_f10) + t * ((1.0f - s) * c_f01 + s * c_f11);
}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Task 7: Implement trilinear filtering

  // return magenta for invalid level
  float dudx = tex.width / u_scale;
  float dvdx = tex.height / u_scale;
  float dudy = tex.width / v_scale;
  float dvdy = tex.height / v_scale;

  float Lx2 = dudx * dudx + dvdx * dvdx;
  float Ly2 = dudy * dudy + dvdy * dvdy;

  float w = max(0.0f, log2f(sqrt(max(Lx2, Ly2)))); // IT IS VERY IMPORTANT TO CHECK BOUNDS
  int d = int(w);
  w = w - d;

  if (d >= tex.mipmap.size())
      return Color(1,0,1,1);
   
  Color a = this->sample_bilinear(tex, u, v, d);
  if ((d + 1) >= tex.mipmap.size())
      return a;

  Color b = this->sample_bilinear(tex, u, v, d + 1);
  return (1.0f - w) * a + (w) * b;

}

} // namespace CMU462
