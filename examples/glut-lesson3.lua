--------------------------------------------------------------------------------
--
-- Lesson 3 of NEHE's OpenGL tutorials, implemented in Lua. Before running this
-- example, type `make lua-glut` to build the GL extension.
--
-- http://nehe.gamedev.net/tutorial/adding_colour/13003
--
--------------------------------------------------------------------------------

require 'luagl'
require 'luaglut'

local ESCAPE = 27
local window = 0


local function InitGL(Width, Height)
   glClearColor(0.0, 0.0, 0.0, 0.0)
   glClearDepth(1.0)
   glDepthFunc(GL_LESS)
   glEnable(GL_DEPTH_TEST)
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
   glClear(GL_COLOR_BUFFER_BIT + GL_DEPTH_BUFFER_BIT)
   glLoadIdentity()
   glTranslated(-1.5, 0.0, -6.0)
   glBegin(GL_POLYGON)
   glColor3d(  1.0, 0.0, 0.0)
   glVertex3d( 0.0, 1.0, 0.0)
   glColor3d(  0.0, 1.0, 0.0)
   glVertex3d( 1.0,-1.0, 0.0)
   glColor3d(  0.0, 0.0, 1.0)
   glVertex3d(-1.0,-1.0, 0.0)
   glEnd()

   glTranslated(3.0, 0.0, 0.0)

   glColor3d(0.5, 0.5, 1.0)
   glBegin(GL_QUADS)
   glVertex3d(-1.0, 1.0, 0.0)
   glVertex3d( 1.0, 1.0, 0.0)
   glVertex3d( 1.0,-1.0, 0.0)
   glVertex3d(-1.0,-1.0, 0.0)
   glEnd()

   glutSwapBuffers()
end

local function keyPressed(key, x, y)
   if key == ESCAPE then
      glutDestroyWindow(window)
      os.exit(0)
   end
end

local function main()
   glutInit(arg)
   glutInitDisplayMode(GLUT_RGBA + GLUT_DOUBLE + GLUT_ALPHA + GLUT_DEPTH)
   glutInitWindowSize(640, 480)
   glutInitWindowPosition(0, 0)
   window = glutCreateWindow("Jeff Molofee's GL Code Tutorial ... NeHe '99")
   glutDisplayFunc(DrawGLScene)
   glutFullScreen()
   glutIdleFunc(DrawGLScene)
   glutReshapeFunc(ReSizeGLScene)
   glutKeyboardFunc(keyPressed)
   InitGL(640, 480)
   glutMainLoop()
   return 1
end

main()
