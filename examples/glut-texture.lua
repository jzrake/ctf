--------------------------------------------------------------------------------
--
-- Demonstrates how to display a simple 2d texture as an image using the GLUT
-- bindings.
--
--------------------------------------------------------------------------------

require 'luagl'
require 'luaglut'
local array = require 'array'

local ESCAPE = 27
local window = 0
local texture = array.vector(1, 'int')

local function InitGL(Width, Height)
   glClearColor(0.0, 0.0, 0.0, 0.0)
   glClearDepth(1.0)
   glDepthFunc(GL_LESS)
   glEnable(GL_DEPTH_TEST)
   glEnable(GL_NORMALIZE)
   glShadeModel(GL_SMOOTH)
   glMatrixMode(GL_PROJECTION)
   glLoadIdentity()
   gluPerspective(45.0, Width/Height, 0.1, 100.0)
   glMatrixMode(GL_MODELVIEW)
end

local function ReSizeGLScene(Width, Height)
   if Height == 0 then Height = 1 end
   glViewport(0, 0, Width, Height)
   glMatrixMode(GL_PROJECTION)
   glLoadIdentity()
   gluPerspective(45.0, Width/Height, 0.1, 100.0)
   glMatrixMode(GL_MODELVIEW)
end

local function DrawGLScene()
   glClear(bit32.bor(GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT))
   glLoadIdentity()

   glTranslated(0.0, 0.0, -3.0)
   glEnable(GL_TEXTURE_2D)
   glBindTexture(GL_TEXTURE_2D, texture[0])
   glBegin(GL_QUADS)
   glTexCoord2d(0, 0); glVertex3d(-1.0, 1.0, 0.0)
   glTexCoord2d(1, 0); glVertex3d( 1.0, 1.0, 0.0)
   glTexCoord2d(1, 1); glVertex3d( 1.0,-1.0, 0.0)
   glTexCoord2d(0, 1); glVertex3d(-1.0,-1.0, 0.0)
   glEnd()
   glBindTexture(GL_TEXTURE_2D, 0)
   glDisable(GL_TEXTURE_2D)

   glutSwapBuffers()
end

local function keyPressed(key, x, y)
   if key == ESCAPE then
      glutDestroyWindow(window)
      os.exit(0)
   end
end

local function LoadTextures()
   local data = array.array({32,32,3}, 'float')
   local datavec = data:vector()

   for i=0,#datavec/3-1 do
      datavec[3*i + 0] = math.random()
      datavec[3*i + 1] = math.random()
      datavec[3*i + 2] = math.random()
   end

   glGenTextures(1, texture:pointer())
   glBindTexture(GL_TEXTURE_2D, texture[0])
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
   glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL)
   glTexImage2D(GL_TEXTURE_2D, 0, 3, 32, 32, 0, GL_RGB, GL_FLOAT, data:pointer())
end

local function main()
   glutInit(arg)
   glutInitDisplayMode(bit32.bor(GLUT_RGBA, GLUT_DOUBLE, GLUT_ALPHA, GLUT_DEPTH))
   glutInitWindowSize(640, 480)
   glutInitWindowPosition(0, 0)
   window = glutCreateWindow("texture example")
   glutDisplayFunc(DrawGLScene)
   glutIdleFunc(DrawGLScene)
   glutReshapeFunc(ReSizeGLScene)
   glutKeyboardFunc(keyPressed)
   InitGL(640, 480)
   LoadTextures()
   glutMainLoop()
   return 1
end

main()
