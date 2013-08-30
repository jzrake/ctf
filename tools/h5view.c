
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <hdf5.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#define ESCAPE 27
#define FILENAME_LENGTH 1024
#define MAX_FILENAMES 4096

static char FileNameArr[MAX_FILENAMES][FILENAME_LENGTH];
static int  FileNameNum = 0;
static int  FileNameCur = 0;
static char *DataDirectory = NULL;
static char *data_set_arr[] = {"rho", "pre", "vx", "vy", "vz", "Bx", "By", "Bz"};
static int   data_set_num = 0;

static GLuint  WindowWidth  = 512;
static GLuint  WindowHeight = 768;
static GLuint  WindowId = 0;
static GLuint  texture_id = 0;

static double  cut_increment = 1e-2;
static double  data_cut[2];
static double *raw_data = NULL;
static int     raw_data_dims[2] = {0,0};
static int     colorbar_id = 3;
static int     log_scale = 0;
static int     user_cuts = 0;
static int     screenshot_lock = 0;
static int     take_screenshot_on_draw = 0;
static int     auto_play = 0;
static int     show_filename = 1;

static void take_screenshot();
static void auto_scale();
static void get_rgb(double val, GLfloat rgb[3]);
static double scale(double z);


void RefreshFileList()
{
  char cmd[1024];
  sprintf(cmd, "ls %s/*", DataDirectory);
  FILE *ls = popen(cmd, "r");
  char c[FILENAME_LENGTH];

  FileNameNum = 0;

  while (1) {
    fgets(c, FILENAME_LENGTH, ls);
    if (feof(ls)) {
      //      fclose(ls);
      break;
    }
    if (strstr(c, ".h5")) {
      if (c[strlen(c)-1] == '\n') {
	c[strlen(c)-1] = '\0';
      }
      strncpy(FileNameArr[FileNameNum], c, FILENAME_LENGTH);
      printf("%s\n", FileNameArr[FileNameNum]);
      ++FileNameNum;
    }
  }
  if (FileNameNum == 0) {
    fprintf(stderr, "ERROR: there are no files to look at\n");
    glutDestroyWindow(WindowId);
    exit(2);
  }
}

void LoadFile()
{
  int i,j;
  const char *fname = FileNameArr[FileNameCur];
  const char *gname = "prim";
  const char *dname = data_set_arr[data_set_num];

  printf("loading %s/prim/%s\n", fname, dname);

  hid_t h5f = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t grp = H5Gopen(h5f, gname, H5P_DEFAULT);
  hid_t set = H5Dopen(grp, dname, H5P_DEFAULT);
  hid_t spc = H5Dget_space(set);

  hsize_t dims[2];
  H5Sget_simple_extent_dims(spc, dims, NULL);
  printf("data set has size [%lld %lld]\n", dims[0], dims[1]);

  raw_data_dims[0] = dims[0];
  raw_data_dims[1] = dims[1];

  double *raw_dataT = (double*) malloc(dims[0] * dims[1] * sizeof(double));
  raw_data = (double*) realloc(raw_data, dims[0] * dims[1] * sizeof(double));
  H5Dread(set, H5T_NATIVE_DOUBLE, spc, spc, H5P_DEFAULT, raw_dataT);

  int Nx = dims[0];
  int Ny = dims[1];

  /* transpose the data */
  for (i=0; i<Nx; ++i) {
    for (j=0; j<Ny; ++j) {
      raw_data[i + j*Nx] = raw_dataT[i*Ny + j];
    }
  }

  free(raw_dataT);
  H5Sclose(spc);
  H5Dclose(set);
  H5Gclose(grp);
  H5Fclose(h5f);
}

void LoadTexture()
{
  int Nx = raw_data_dims[0];
  int Ny = raw_data_dims[1];
  int i,j,m;
  double min, max, z;
  GLfloat *texture_data = (GLfloat*) malloc(Nx * Ny * 3 * sizeof(GLfloat));

  if (!user_cuts) {
    auto_scale();
  }

  min = data_cut[0];
  max = data_cut[1];

  printf("data cuts [%4.3f, %4.3f] %s-scaled\n", min, max,
	 log_scale ? "log" : "linearly");

  for (i=0; i<Nx; ++i) {
    for (j=0; j<Ny; ++j) {
      m = i*Ny + j;
      z = (scale(raw_data[m]) - min) / (max - min);
      get_rgb(z, &texture_data[3*m]);
    }
  }

  glBindTexture(GL_TEXTURE_2D, texture_id);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
  glTexImage2D(GL_TEXTURE_2D, 0, 3, Nx, Ny, 0, GL_RGB, GL_FLOAT, texture_data);
  glBindTexture(GL_TEXTURE_2D, 0);
  free(texture_data);
}

void glutPrint(float x, float y, const char *text)
{
  int blending = glIsEnabled(GL_BLEND);
  if(!text || !strlen(text)) return;
  glEnable(GL_BLEND);
  glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
  glRasterPos3f(x, y, 1e-3);
  while (*text) {
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, *text);
    text++;
  }
  if (!blending) glDisable(GL_BLEND);
}

void InitGL(int Width, int Height)
{
  glClearColor(0.9f, 0.9f, 0.9f, 0.0f);
  glClearDepth(1.0);
  glDepthFunc(GL_LESS);
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);
  glMatrixMode(GL_MODELVIEW);
}

void ReSizeGLScene(int Width, int Height)
{
  if (Height==0) Height=1;
  glViewport(0, 0, Width, Height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);
  glMatrixMode(GL_MODELVIEW);
}

void DrawGLScene()
{
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
  glLoadIdentity();

  glTranslated(0, 0, -5.0);

  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, texture_id);
  glBegin(GL_QUADS);
  glTexCoord2d(0, 0); glVertex3d(-1.0, 2.0, 0.0);
  glTexCoord2d(1, 0); glVertex3d( 1.0, 2.0, 0.0);
  glTexCoord2d(1, 1); glVertex3d( 1.0,-2.0, 0.0);
  glTexCoord2d(0, 1); glVertex3d(-1.0,-2.0, 0.0);
  glEnd();
  glBindTexture(GL_TEXTURE_2D, 0);
  glDisable(GL_TEXTURE_2D);

  if (show_filename) {
    glutPrint(-1.0 + 1e-2, -2.0 + 1e-2, FileNameArr[FileNameCur]);
  }

  if (take_screenshot_on_draw) {
    take_screenshot_on_draw = 0;
    take_screenshot();
  }
  if (auto_play) {
    if (FileNameCur < FileNameNum - 1) {
      ++FileNameCur;
      LoadFile();
      LoadTexture();
      if (screenshot_lock) {
	take_screenshot_on_draw = 1;
      }
    }
  }

  glutSwapBuffers();
}

void keyPressed(unsigned char key, int x, int y)
{
  switch (key) {
  case ESCAPE:
    glutDestroyWindow(WindowId);
    exit(0);
  case 'r':
    RefreshFileList();
    break;
  case ']':
    if (FileNameCur < FileNameNum - 1) {
      ++FileNameCur;
      LoadFile();
      LoadTexture();
      if (screenshot_lock) {
	take_screenshot_on_draw = 1;
      }
    }
    break;
  case '}':
    auto_play ^= 1;
    break;
  case '1': data_set_num = 0; LoadFile(); LoadTexture(); break;
  case '2': data_set_num = 1; LoadFile(); LoadTexture(); break;
  case '3': data_set_num = 2; LoadFile(); LoadTexture(); break;
  case '4': data_set_num = 3; LoadFile(); LoadTexture(); break;
  case '5': data_set_num = 4; LoadFile(); LoadTexture(); break;
  case '6': data_set_num = 5; LoadFile(); LoadTexture(); break;
  case '7': data_set_num = 6; LoadFile(); LoadTexture(); break;
  case '8': data_set_num = 7; LoadFile(); LoadTexture(); break;
  case '[':
    if (FileNameCur > 0) {
      --FileNameCur;
      LoadFile();
      LoadTexture();
      if (screenshot_lock) {
	take_screenshot_on_draw = 1;
      }
    }
    break;
  case 'c':
    ++colorbar_id;
    if (colorbar_id == 5) colorbar_id = 0;
    printf("switch to color bar %d\n", colorbar_id);
    LoadTexture();
    break;
  case 'l':
    log_scale = !log_scale;
    LoadTexture();
    break;
  case 'a':
    user_cuts = 0; // auto-scale data by removing user cuts
    LoadTexture();
    break;
  case ',':
    user_cuts = 1;
    data_cut[0] -= cut_increment * (data_cut[1] - data_cut[0]);
    LoadTexture();
    break;
  case '.':
    user_cuts = 1;
    data_cut[0] += cut_increment * (data_cut[1] - data_cut[0]);
    LoadTexture();
    break;
  case '<':
    user_cuts = 1;
    data_cut[1] -= cut_increment * (data_cut[1] - data_cut[0]);
    LoadTexture();
    break;
  case '>':
    user_cuts = 1;
    data_cut[1] += cut_increment * (data_cut[1] - data_cut[0]);
    LoadTexture();
    break;
  case 's':
    take_screenshot();
    break;
  case 'S':
    screenshot_lock ^= 1;
    printf("%s screenshot lock\n", screenshot_lock ? "enabling" : "disabling");
    break;
  case 'f':
    show_filename ^= 1;
    break;
  }
}

int main(int argc, char **argv)
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
  glutInitWindowSize(WindowWidth, WindowHeight);
  glutInitWindowPosition(512, 0);
  WindowId = glutCreateWindow("Flying Grid");
  glutDisplayFunc(&DrawGLScene);
  //  glutFullScreen();
  glutIdleFunc(&DrawGLScene);
  glutReshapeFunc(&ReSizeGLScene);
  glutKeyboardFunc(&keyPressed);
  InitGL(WindowWidth, WindowHeight);

  if (argc == 1) {
    DataDirectory = "data";
  }
  else {
    DataDirectory = argv[1];
  }

  glGenTextures(1, &texture_id);

  RefreshFileList();
  LoadFile();
  LoadTexture();

  glutMainLoop();

  return 0;
}

void auto_scale()
{
  int Nx = raw_data_dims[0];
  int Ny = raw_data_dims[1];
  int m,i,j;
  double z;

  data_cut[0] = scale(raw_data[0]);
  data_cut[1] = scale(raw_data[0]);

  for (i=0; i<Nx; ++i) {
    for (j=0; j<Ny; ++j) {
      m = i*Ny + j;
      z = scale(raw_data[m]);
      if (z < data_cut[0]) data_cut[0] = z;
      if (z > data_cut[1]) data_cut[1] = z;
    }
  }
}

double scale(double z)
{
  if (log_scale) {
    return log10(z);
  }
  else {
    return z;
  }
}

void get_rgb(double val, GLfloat rgb[3])
{
  double rrr, ggg, bbb;

  if( colorbar_id == 0 ){
    double nexp = 8.0;
    rrr = exp(-nexp*pow(val-5./6.,2.0)) + .25*exp(-nexp*pow(val+1./6.,2.0));
    ggg = exp(-nexp*pow(val-3./6.,2.0));
    bbb = exp(-nexp*pow(val-1./6.,2.0)) + .25*exp(-nexp*pow(val-7./6.,2.0));
  }else if(colorbar_id == 1){
    if( val < .1 ){
      bbb = 4.*(val+.15);
      ggg = 0.0;
      rrr = 0.0;
    }else if( val < .35){
      bbb = 1.0;
      ggg = 4.*(val-.1);
      rrr = 0.0;
    }else if( val < .6 ){
      bbb = 4.*(.6-val);
      ggg = 1.;
      rrr = 4.*(val-.35);
    }else if( val < .85){
      bbb = 0.0;
      ggg = 4.*(.85-val);
      rrr = 1.;
    }else{
      bbb = 0.0;
      ggg = 0.0;
      rrr = 4.*(1.1-val);
    }
  }else if(colorbar_id == 2){
    rrr = 2.*val;
    ggg = 1.2*val;
    bbb = .8*val;
  }else if(colorbar_id == 3){
    double gam = .8;
    double Amp;
    double r0,g0,b0;
    double hi,lo,x1,x2,x3,x4;
    hi = .8;
    lo = .1;
    if( val > hi ) Amp = .3 + .7*(1.-val)/(1.-hi);
    else if( val < lo ) Amp = .3 + .7*(val)/(lo);
    else Amp = 1.0;

    x1 = .5;
    x2 = .325;
    x3 = .15;
    x4 = 0.;

    if( val > x1 )      r0 = 1.;
    else if( val > x2 ) r0 = (val-x2)/(x1-x2);
    else if( val > x3 ) r0 = 0.;
    else if( val > x4 ) r0 = (val-x3)/(x4-x3);
    else                r0 = 1.;

    x1 = .6625;
    x2 = .5;
    x3 = .275;
    x4 = .15;

    if( val > x1 )      g0 = 0.;
    else if( val > x2 ) g0 = (val-x1)/(x2-x1);
    else if( val > x3 ) g0 = 1.;
    else if( val > x4 ) g0 = (val-x4)/(x3-x4);
    else                g0 = 0.;

    x1 = .325;
    x2 = .275;

    if( val > x1 )      b0 = 0.;
    else if( val > x2 ) b0 = (val-x1)/(x2-x1);
    else                b0 = 1.;

    rrr = pow(Amp*r0,gam);
    ggg = pow(Amp*g0,gam);
    bbb = pow(Amp*b0,gam);
  }else if(colorbar_id == 4){
    if( val < .1 ){
      bbb = 4.*(val+.125);
      ggg = 0.0;
      rrr = 0.0;
    }else if( val < .375){
      bbb = 1.0;
      ggg = 4.*(val-.125);
      rrr = 0.0;
    }else if( val < .625 ){
      bbb = 4.*(.625-val);
      rrr = 4.*(val-.375);
      ggg = bbb;
      if( rrr > bbb ) ggg = rrr;
    }else if( val < .875){
      bbb = 0.0;
      ggg = 4.*(.875-val);
      rrr = 1.;
    }else{
      bbb = 0.0;
      ggg = 0.0;
      rrr = 4.*(1.125-val);
    }
  }else if(colorbar_id == 5){
    rrr = val;
    ggg = val;
    bbb = val;
  }else{
    rrr = 1.0;
    ggg = 1.0;
    bbb = 1.0;
  }

  rgb[0] = rrr;
  rgb[1] = ggg;
  rgb[2] = bbb;
}

void take_screenshot()
{
  static int num = 0;
  char fname[1024];
  char convert[1024];
  num += 1;

  sprintf(fname, "%04d.ppm", num);
  sprintf(convert, "convert -flip %04d.ppm %04d.png; rm %04d.ppm", num, num, num);

  printf("writing a screenshot to %s\n", fname);
 
  int dimx = glutGet(GLUT_WINDOW_WIDTH);
  int dimy = glutGet(GLUT_WINDOW_HEIGHT);
 
  size_t imsize = 3*dimx*dimy;
  char *pixels = (char*) malloc(imsize*sizeof(char));
  glReadPixels(0, 0, dimx, dimy, GL_RGB, GL_UNSIGNED_BYTE, pixels);
 
  FILE *fp = fopen(fname, "wb");
  fprintf(fp, "P6\n%d %d\n255\n", dimx, dimy);
  fwrite(pixels, sizeof(char), imsize, fp);
  fclose(fp);
 
  free(pixels);
  system(convert);
}
