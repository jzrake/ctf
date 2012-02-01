

#include "config.h"
#if (__MARA_USE_GLFW)
#include <stdlib.h>
#include <math.h>
#include "vis.hpp"
#include "valman.hpp"
extern "C" {
#include "lualib.h"
#include "lauxlib.h"
#include "GL/glfw.h"
}


static int WindowWidth  = 768;
static int WindowHeight = 768;

static float xTranslate = 0.0;
static float yTranslate = 0.0;
static float zTranslate = 1.4;

static float RotationAngleX = 220;
static float RotationAngleY =   0;
static float ZoomFactor = 1.0;
static GLuint TextureMap;

static int VariableIndex = 0;

static void KeyboardInput(int key, int state);
static void CharacterInput(int key, int state);
static void OpenWindow();


VisualizationOpenGl::VisualizationOpenGl()
{
  OpenWindow();
  glGenTextures(1, &TextureMap);
}

VisualizationOpenGl::~VisualizationOpenGl()
{
  glDeleteTextures(1, &TextureMap);
  glfwTerminate();
}

void VisualizationOpenGl::load_texture()
{
  static const int Nx = Mara->domain->aug_shape()[0];
  static const int Ny = Mara->domain->aug_shape()[1];

  GLfloat *TextureData = (GLfloat*) malloc(Nx*Ny*3*sizeof(GLfloat));

  const int sx = 1;
  const int sy = Nx;

  double zmax=-1e16, zmin=1e16;

  ValarrayManager M(Mara->domain->aug_shape(), Mara->domain->get_Nq());
  double *zdata = (double*) malloc(Nx*Ny*sizeof(double));

  for (int i=0; i<Nx; ++i) {
    for (int j=0; j<Ny; ++j) {

      std::valarray<double> P0 = Mara->PrimitiveArray[ M(i,j) ];
      const double z = P0[VariableIndex];

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


void VisualizationOpenGl::DrawScene()
{
  if (Mara->domain->get_Nd() != 2) return;

  const double Lx0 = Mara->domain->get_x0()[0];
  const double Lx1 = Mara->domain->get_x1()[0];

  const double Ly0 = Mara->domain->get_x0()[1];
  const double Ly1 = Mara->domain->get_x1()[1];

  xTranslate = -0.5*(Lx0 + Lx1);
  yTranslate = -0.5*(Ly0 + Ly1);

  this->load_texture();

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

    char keynum[] = "12345678";

    for (int i=0; i<Mara->domain->get_Nq(); ++i) {
      if (glfwGetKey(keynum[i])) {
	VariableIndex = i;
	this->load_texture();
	break;
      }
    }

    if (glfwGetKey(GLFW_KEY_SPACE)) {
      break;
    }
    if (glfwGetKey(GLFW_KEY_ESC) || !glfwGetWindowParam(GLFW_OPENED)) {
      glfwCloseWindow();

      Mara->visual = NULL;
      delete this;
      break;
    }
  }
}




void OpenWindow()
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

  default:
    break;
  }
}

#else

#include "vis.hpp"
VisualizationOpenGl::VisualizationOpenGl() { }
VisualizationOpenGl::~VisualizationOpenGl() { }
void VisualizationOpenGl::load_texture() { }
void VisualizationOpenGl::DrawScene() { }

#endif // __MARA_USE_GLFW
