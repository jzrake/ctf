

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "lualib.h"
#include "lauxlib.h"
#include "lunum.h"

#include "GL/glfw.h"



static int open_window(lua_State *L);
static int draw_texture(lua_State *L);



int lua_vis_load(lua_State *L)
{
  luaL_Reg vis_api[] = { { "open_window", open_window },
			 { "draw_texture", draw_texture },
			 { NULL, NULL} };

  lua_newtable(L);
  luaL_setfuncs(L, vis_api, 0);
  lua_setglobal(L, "visual");

  return 0;
}


static int Autoplay     = 0;
static int WindowOpen   = 0;
static int WindowWidth  = 768;
static int WindowHeight = 768;

static float xTranslate = 0.0;
static float yTranslate = 0.0;
static float zTranslate = 1.4;

static float RotationAngleX = 220;
static float RotationAngleY =   0;
static float ZoomFactor = 1.0;
static GLuint TextureMap;


static void LoadTexture(lua_State *L);
static void KeyboardInput(int key, int state);
static void CharacterInput(int key, int state);




int open_window(lua_State *L)
{
  glfwInit();
  glfwOpenWindow(WindowWidth, WindowHeight, 0,0,0,0,0,0, GLFW_WINDOW);
  //  glfwEnable(GLFW_STICKY_KEYS);
  //  glfwEnable(GLFW_KEY_REPEAT);
  glfwSetKeyCallback(KeyboardInput);
  glfwSetCharCallback(CharacterInput);

  glClearColor(0.2, 0.1, 0.1, 0.0);
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

  LoadTexture(L);

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


void LoadTexture(lua_State *L)
{
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

  for (int i=0; i<Nx; ++i) {
    for (int j=0; j<Ny; ++j) {

      const int m = i*sx + j*sy;
      const float v = (zdata[m] - zmin) / (zmax-zmin);

      float ColorWidth = 0.05;
      float A = pow(0.5 / ColorWidth, 2);
      float r = exp(-A*pow(v-0.1,2)) + exp(-A*pow(v-0.4,2)) + exp(-A*pow(v-0.7,2));
      float g = exp(-A*pow(v-0.2,2)) + exp(-A*pow(v-0.5,2)) + exp(-A*pow(v-0.8,2));
      float b = exp(-A*pow(v-0.3,2)) + exp(-A*pow(v-0.6,2)) + exp(-A*pow(v-0.9,2));

      TextureData[3*m + 0] = r;
      TextureData[3*m + 1] = g;
      TextureData[3*m + 2] = b;
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

void CharacterInput(int key, int state)
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

  default:
    break;
  }
}
