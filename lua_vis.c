
#include "lualib.h"
#include "config.h"
#if (__MARA_USE_GLFW)


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "image.h"
#include "lauxlib.h"
#include "lunum.h"

#include "GL/glfw.h"


static int open_window(lua_State *L);
static int draw_texture(lua_State *L);
static int draw_lines3d(lua_State *L);



void lua_vis_load(lua_State *L)
{
  luaL_Reg vis_api[] = { { "open_window", open_window },
                         { "draw_texture", draw_texture },
                         { "draw_lines3d", draw_lines3d },
                         { NULL, NULL} };

  lua_newtable(L);
  luaL_setfuncs(L, vis_api, 0);
  lua_setglobal(L, "visual");
}



static int Autoplay     = 0;
static int WindowOpen   = 0;
static int WindowWidth  = 1024;
static int WindowHeight = 768;

static float xTranslate = 0.0;
static float yTranslate = 0.0;
static float zTranslate = 1.4;

static float RotationAngleX = 10;
static float RotationAngleY =  0;
static float ZoomFactor = 1.0;
static GLuint TextureMap;
static int ColormapIndex = 0;


static void KeyboardInput(int key, int state);
static void CharacterInput_draw_texture(int key, int state);
static void CharacterInput_draw_lines3d(int key, int state);
static void LoadTexture2d();
static void DrawBoundingBox(const double *bounds);


static lua_State *PresentLua = NULL;


int open_window(lua_State *L)
{
  double ClearColor[3] = { 0.1, 0.2, 0.2 };

  if (lua_gettop(L) > 0) {

    lua_getfield(L, 1, "window_size");
    if (!lua_isnil(L, -1)) {
      lua_rawgeti(L, -1, 1); WindowWidth  = lua_tonumber(L, -1); lua_pop(L, 1);
      lua_rawgeti(L, -1, 2); WindowHeight = lua_tonumber(L, -1); lua_pop(L, 1);
    }
    lua_pop(L, 1);

    lua_getfield(L, 1, "clear_color");
    if (!lua_isnil(L, -1)) {
      lua_rawgeti(L, -1, 1); ClearColor[0] = lua_tonumber(L, -1); lua_pop(L, 1);
      lua_rawgeti(L, -1, 2); ClearColor[1] = lua_tonumber(L, -1); lua_pop(L, 1);
      lua_rawgeti(L, -1, 3); ClearColor[2] = lua_tonumber(L, -1); lua_pop(L, 1);
    }
    lua_pop(L, 1);
  }
  printf("opening window (%d %d)\n", WindowWidth, WindowHeight);

  glfwInit();
  glfwOpenWindow(WindowWidth, WindowHeight, 0,0,0,0,0,0, GLFW_WINDOW);

  glClearColor(ClearColor[0], ClearColor[1], ClearColor[2], 0.0);
  glClearDepth(1.0);
  glDepthFunc(GL_LESS);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);
  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0, (float) WindowWidth / WindowHeight, 0.1, 100.0);
  glMatrixMode(GL_MODELVIEW);

  WindowOpen = 1;

  return 0;
}


int draw_texture(lua_State *L)
{
  if (!WindowOpen) {
    luaL_error(L, "there is no open window to draw in");
  }

  glfwSetKeyCallback(NULL);
  glfwSetCharCallback(CharacterInput_draw_texture);

  glfwDisable(GLFW_STICKY_KEYS);
  glfwDisable(GLFW_KEY_REPEAT);
  zTranslate = 1.4;

  PresentLua = L;
  LoadTexture2d();

  const double Lx0 = -0.5;
  const double Lx1 = +0.5;

  const double Ly0 = -0.5;
  const double Ly1 = +0.5;

  xTranslate = -0.5*(Lx0 + Lx1);
  yTranslate = -0.5*(Ly0 + Ly1);

  while (1) {

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    glTranslatef(-xTranslate, -yTranslate, -zTranslate);
    glScalef(ZoomFactor, ZoomFactor, ZoomFactor);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, TextureMap);

    glBegin(GL_QUADS);
    glNormal3f(0, 0, 1);
    glTexCoord2f(0, 0); glVertex3f(Lx0, Ly0, 0);
    glTexCoord2f(0, 1); glVertex3f(Lx0, Ly1, 0);
    glTexCoord2f(1, 1); glVertex3f(Lx1, Ly1, 0);
    glTexCoord2f(1, 0); glVertex3f(Lx1, Ly0, 0);
    glEnd();

    glFlush();
    glfwSwapBuffers();

    if (glfwGetKey(GLFW_KEY_ESC) || !glfwGetWindowParam(GLFW_OPENED)) {
      glfwCloseWindow();
      WindowOpen = 0;
      break;
    }

    if (Autoplay || glfwGetKey(GLFW_KEY_SPACE)) {
      break;
    }
  }
  return 0;
}
void CharacterInput_draw_texture(int key, int state)
{
  switch (key) {

  case 'z':
    ZoomFactor /= 1.1;
    break;

  case 'Z':
    ZoomFactor *= 1.1;
    break;

  case 'p':
    Autoplay ^= 1;
    break;

  case 'c':
    if (Mara_image_get_colormap(++ColormapIndex) == NULL) {
      ColormapIndex = 0;
    }
    LoadTexture2d();
    break;

  default:
    break;
  }
}

void _draw_single_line(double *points, int N)
{
  double smax=-1e16, smin=1e16;
  for (int m=0; m<N; ++m) {
    const double s = points[4*m + 3];
    if (s < smin) smin = s;
    if (s > smax) smax = s;
  }

  const float *cmap_data = Mara_image_get_colormap(ColormapIndex);

  glBegin(GL_LINE_STRIP);
  for (int m=0; m<N; ++m) {
    const double *x = &points[4*m];

    int cm = 255.0 * (x[3] - smin) / (smax - smin);
    double r = cmap_data[3*cm + 0];
    double g = cmap_data[3*cm + 1];
    double b = cmap_data[3*cm + 2];

    glColor3d(r, g, b);
    glVertex3d(x[0], x[1], x[2]);
  }
  glEnd();
}

int draw_lines3d(lua_State *L)
// -----------------------------------------------------------------------------
// @input : either a (N x 4) lunum array, or an integer-key table of such things
// @output: nil
// -----------------------------------------------------------------------------
{
  if (!WindowOpen) {
    luaL_error(L, "there is no open window to draw in");
  }

  struct Array **lines = NULL;
  int nlines;

  if (lunum_hasmetatable(L, 1, "array")) {
    nlines = 1;
    lines = (struct Array**) malloc(nlines*sizeof(struct Array*));
    lines[0] = lunum_checkarray1(L, 1);
  }
  else {
    nlines = lua_rawlen(L, 1);
    lines = (struct Array**) malloc(nlines*sizeof(struct Array*));

    for (int n=0; n<nlines; ++n) {
      lua_rawgeti(L, 1, n+1);
      struct Array *A = lunum_checkarray1(L, 2);
      lua_pop(L, 1);
      lines[n] = A;
    }
  }

  glfwSetKeyCallback(KeyboardInput);
  glfwSetCharCallback(CharacterInput_draw_lines3d);

  glfwEnable(GLFW_STICKY_KEYS);
  glfwEnable(GLFW_KEY_REPEAT);

  zTranslate = 1.8;
  PresentLua = L;

  xTranslate = 0.0;
  yTranslate = 0.0;

  const double bounds[] = { -0.5, +0.5,
                            -0.5, +0.5,
                            -0.5, +0.5 };
  while (1) {

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    glTranslatef(-xTranslate, -yTranslate, -zTranslate);
    glScalef(ZoomFactor, ZoomFactor, ZoomFactor);

    glRotatef(RotationAngleX, 1, 0, 0);
    glRotatef(RotationAngleY, 0, 1, 0);

    DrawBoundingBox(bounds);
    for (int n=0; n<nlines; ++n) {

      struct Array *A = lines[n];

      if (A->dtype != ARRAY_TYPE_DOUBLE || A->ndims != 2) {
	free(lines);
	luaL_error(L, "need a (N x 4) array of doubles");
	return 0;
      }
      if (A->shape[1] != 4) {
	free(lines);
	luaL_error(L, "need a (N x 4) array of doubles");
	return 0;
      }

      _draw_single_line((double*)A->data, A->shape[0]);
    }

    glFlush();
    glfwSwapBuffers();

    if (glfwGetKey(GLFW_KEY_ESC) || !glfwGetWindowParam(GLFW_OPENED)) {
      glfwCloseWindow();
      WindowOpen = 0;
      break;
    }

    if (Autoplay || glfwGetKey(GLFW_KEY_SPACE)) {
      break;
    }
  }

  free(lines);
  return 0;
}
void CharacterInput_draw_lines3d(int key, int state)
{
  switch (key) {

  case 'z':
    ZoomFactor /= 1.1;
    break;

  case 'Z':
    ZoomFactor *= 1.1;
    break;

  case 'p':
    Autoplay ^= 1;
    break;

  case 'c':
    if (Mara_image_get_colormap(++ColormapIndex) == NULL) {
      ColormapIndex = 0;
    }
    break;

  default:
    break;
  }
}

void LoadTexture2d()
{
  lua_State *L = PresentLua;

  if (lunum_upcast(L, 1, ARRAY_TYPE_DOUBLE, 0)) {
    lua_replace(L, -2);
  }
  struct Array *A = lunum_checkarray1(L, 1);

  const int Nx = A->shape[0];
  const int Ny = A->shape[1];

  const int sx = 1; // apply to graphics card data, not input data
  const int sy = Nx;

  double zmax=-1e16, zmin=1e16;

  GLfloat *TextureData = (GLfloat*) malloc(Nx*Ny*3*sizeof(GLfloat));
  double *zdata = (double*) malloc(Nx*Ny*sizeof(double));
  double *data = (double*) A->data;

  for (int i=0; i<Nx; ++i) {
    for (int j=0; j<Ny; ++j) {

      const double z = data[i*Ny + j];

      if (z < zmin) zmin = z;
      if (z > zmax) zmax = z;

      zdata[i*sx + j*sy] = z;
    }
  }

  const float *cmap_data = Mara_image_get_colormap(ColormapIndex);

  for (int i=0; i<Nx; ++i) {
    for (int j=0; j<Ny; ++j) {

      const int m = i*sx + j*sy;
      int cm = 255.0 * (zdata[m] - zmin) / (zmax - zmin);

      TextureData[3*m + 0] = cmap_data[3*cm + 0];
      TextureData[3*m + 1] = cmap_data[3*cm + 1];
      TextureData[3*m + 2] = cmap_data[3*cm + 2];
    }
  }

  glBindTexture(GL_TEXTURE_2D, TextureMap);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, Nx, Ny, 0, GL_RGB, GL_FLOAT, TextureData);

  free(TextureData);
  free(zdata);
}


void KeyboardInput(int key, int state)
{
  if (state != GLFW_PRESS) return;

  switch (key) {
  case GLFW_KEY_RIGHT : RotationAngleY += 3; break;
  case GLFW_KEY_LEFT  : RotationAngleY -= 3; break;

  case GLFW_KEY_DOWN  : RotationAngleX += 3; break;
  case GLFW_KEY_UP    : RotationAngleX -= 3; break;
  }
}



void DrawBoundingBox(const double *bounds)
{
  const double Lx0 = bounds[0];
  const double Lx1 = bounds[1];

  const double Ly0 = bounds[2];
  const double Ly1 = bounds[3];

  const double Lz0 = bounds[4];
  const double Lz1 = bounds[5];

  // Drawing a bonding box
  // ---------------------------------------------------------------------------
  glColor3d(0.6, 0.3, 0.3);
  glLineWidth(3.0);

  glBegin(GL_LINES);

  // x-edges
  // ---------------------------------------------------------------------------
  glVertex3f(Lx0, Ly0, Lz0);
  glVertex3f(Lx1, Ly0, Lz0);

  glVertex3f(Lx0, Ly0, Lz1);
  glVertex3f(Lx1, Ly0, Lz1);

  glVertex3f(Lx0, Ly1, Lz0);
  glVertex3f(Lx1, Ly1, Lz0);

  glVertex3f(Lx0, Ly1, Lz1);
  glVertex3f(Lx1, Ly1, Lz1);

  // y-edges
  // ---------------------------------------------------------------------------
  glVertex3f(Lx0, Ly0, Lz0);
  glVertex3f(Lx0, Ly1, Lz0);

  glVertex3f(Lx1, Ly0, Lz0);
  glVertex3f(Lx1, Ly1, Lz0);

  glVertex3f(Lx0, Ly0, Lz1);
  glVertex3f(Lx0, Ly1, Lz1);

  glVertex3f(Lx1, Ly0, Lz1);
  glVertex3f(Lx1, Ly1, Lz1);

  // z-edges
  // ---------------------------------------------------------------------------
  glVertex3f(Lx0, Ly0, Lz0);
  glVertex3f(Lx0, Ly0, Lz1);

  glVertex3f(Lx0, Ly1, Lz0);
  glVertex3f(Lx0, Ly1, Lz1);

  glVertex3f(Lx1, Ly0, Lz0);
  glVertex3f(Lx1, Ly0, Lz1);

  glVertex3f(Lx1, Ly1, Lz0);
  glVertex3f(Lx1, Ly1, Lz1);

  glEnd();
}


#else
void lua_vis_load(lua_State *L) { }
#endif // __MARA_USE_GLFW
